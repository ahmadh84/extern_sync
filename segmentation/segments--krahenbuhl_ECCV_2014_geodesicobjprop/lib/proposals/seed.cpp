/*
    Copyright (c) 2014, Philipp Krähenbühl
    All rights reserved.
	
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the Stanford University nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.
	
    THIS SOFTWARE IS PROVIDED BY Philipp Krähenbühl ''AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL Philipp Krähenbühl BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
	 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
	 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
	 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include "seed.h"
#include "util.h"
#include "saliency.h"
#include "util/algorithm.h"
#include "util/optimization.h"
#include <queue>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <random>

SeedFunction::~SeedFunction() {
}
void SeedFunction::save(std::ostream &s) const {
}
void SeedFunction::load(std::istream &s) {
}
VectorXi ImageSeedFunction::compute(const OverSegmentation &os, int M) const {
	try {
		return computeImageOverSegmentation( dynamic_cast<const ImageOverSegmentation&>( os ), M );
	}
	catch (std::bad_cast e) {
		throw std::invalid_argument("Image seed function not supported for ImageOverSegmentation.");
	}
}

VectorXi RegularSeed::computeImageOverSegmentation(  const ImageOverSegmentation & ios, int M ) const {
	// Pick the seeds on a regular grid (only supported on images)
	VectorXi r( M );
	const RMatrixXs & s = ios.s();
	int ny = std::max( (int)(sqrt(1.0*M*s.rows()/s.cols())+0.5), 1 );
	for( int j=0,k=0; j<ny; j++ ) {
		int i0 = j*M/ny, i1 = (j+1)*M/ny;
		for( int i=i0; i<i1; i++,k++ )
			r[k] = s( (int)((j+0.5)/ny*s.rows()), (int)((i-i0+0.5)/(i1-i0)*s.cols()) );
	}
	return r;
}
GeodesicSeed::GeodesicSeed( float pow, float const_w, float min_d ) : pow_(pow), const_w_(const_w), min_d_( min_d ) {
}
VectorXi GeodesicSeed::compute( const OverSegmentation & os, int M ) const {
	VectorXi seeds( M );
	// Place the first seed in the geodesic center
	const Edges & edges = os.edges();
	VectorXf ew = os.edgeWeights().array().pow( pow_ ) + const_w_;
	seeds[0] = geodesicCenter( edges, ew );
	GeodesicDistance gdist( edges, ew );
	gdist.update( seeds[0] );
	// Place all other seeds at the max geodesic distance
	for( int it=1; it<M; it++ ) {
		float d = gdist.d().maxCoeff( &seeds[it] );
		if( d <= min_d_ ) {
			// Stop early
			return seeds.head( it );
		}
		if( it+1<M )
			gdist.update( seeds[it] );
	}
	return seeds;
}
void GeodesicSeed::save(std::ostream &s) const {
	float param[3] = {pow_, const_w_, min_d_};
	s.write( (const char*)param, sizeof(param) );
}
void GeodesicSeed::load(std::istream &s) {
	float param[3];
	s.read( (char*)param, sizeof(param) );
	pow_ = param[0];
	const_w_ = param[1];
	min_d_ = param[2];
}
VectorXi RandomSeed::compute( const OverSegmentation & os, int M ) const {
	// Pick the seeds randomly
	return randomChoose( os.Ns(), M );
}
VectorXi SegmentationSeed::compute( const OverSegmentation & os, int M ) const {
	const float h_weight = 0.2;
	const Edges & edges = os.edges();
	VectorXf ew = os.edgeWeights();
	//// Segment the image a bit ////
	// Run a simple hierachical clustering algorithm to generate the seed set
	std::priority_queue< Node, std::vector<Node>, std::greater<Node> > q;
	for (int i=0; i<edges.size(); i++ )
		q.push(Node( i, ew[i]));
	
	UnionFindSet ufset( os.Ns() );
	std::vector<int> height(  os.Ns(), 0 );
	for( int m=M; m< os.Ns() && !q.empty(); ) {
		Node e = q.top();
		q.pop();
		int a = ufset.find( edges[e.to].a ), b = ufset.find( edges[e.to].b );
		if( a != b ) {
			float merge_cost = ew[e.to] + h_weight*std::max(height[a],height[b]);
			if( e.w < merge_cost )
				q.push(Node( e.to, merge_cost ) );
			else {
				ufset.merge(a,b);
				height[a] = height[b] = std::max( height[a], height[b])+1;
				m++;
			}
		}
	}
	VectorXi seg = -VectorXi::Ones( os.Ns());
	int nseg=0;
	for( int i=0; i< os.Ns(); i++ ) {
		if( seg[ufset.find(i)] == -1 )
			seg[i] = seg[ufset.find(i)] = nseg++;
		else
			seg[i] = seg[ufset.find(i)];
	}
	
	//// Segment the image a bit ////
	
	// Compute the graph of each segment
	std::vector< Edges > graphs( nseg );
	std::vector< std::vector<float> > weights( nseg );
	std::vector< std::vector<int> > ids( nseg );
	for( int i=0; i<seg.size(); i++ )
		ids[ seg[i] ].push_back( i );
	for( int i=0; i<edges.size(); i++ )
		if( seg[edges[i].a] == seg[edges[i].b] ) {
			graphs [seg[edges[i].a]].push_back( edges[i] );
			weights[seg[edges[i].a]].push_back( ew[i] );
		}
	// Find the geodesic centers of the segments
	VectorXi r( nseg );
	for( int s=0; s<nseg; s++ ) {
		// Create the subgraph of segment s
		std::unordered_map<int,int> id_map;
		int nid=0;
		for( int i: ids[s] )
			id_map[i] = nid++;
		
		// Compress the graph
		Edges mappend_graph;
		for( int i=0; i<graphs[s].size(); i++ )
			mappend_graph.push_back( Edge( id_map[ graphs[s][i].a ], id_map[ graphs[s][i].b ] ) );
		
		// Find the geodesic center
		int c = geodesicCenter( mappend_graph, VectorXf::Map( weights[s].data(), weights[s].size() ) );
		r[s] = ids[s][c];
	}
	return r;
}

LearnedSeed::LearnedSeed() {
	// Basic features
	features_.addPosition();
	features_.addColor();
	// Weighted geodesics
	features_.addGeodesic(0,1,1);
	features_.addGeodesic(1,1,2e-3);
// 	features_.addGeodesic(1,2,2e-3);
	features_.addGeodesic(1,3,2e-3);
	// Geodesics to bnd
	features_.addGeodesicBnd(0,1,1);
	features_.addGeodesicBnd(1,1,2e-3);
	features_.addGeodesicBnd(1,3,2e-3);
	
}
void LearnedSeed::train( const std::vector< ImageOverSegmentation > &ios, const std::vector<VectorXs> & lbl, const std::vector<VectorXf> & weight, int max_feature ) {
	std::vector< SeedFeatureVector > f;
	for( const auto & g: ios )
		f.push_back( features_.create( g ) );
	train( f, lbl, weight, max_feature );
}
void LearnedSeed::train( const std::vector< ImageOverSegmentation* > &ios, const std::vector<VectorXs> & lbl, const std::vector<VectorXf> & weight, int max_feature ) {
	std::vector< SeedFeatureVector > f;
	for( const auto * g: ios )
		f.push_back( features_.create( *g ) );
	train( f, lbl, weight, max_feature );
}
void LearnedSeed::train( const std::vector< ImageOverSegmentation > &ios, const std::vector<VectorXs> & lbl, int max_feature ) {
	std::vector<VectorXf> weight( lbl.size() );
	for( int i=0; i<lbl.size(); i++ )
		weight[i] = VectorXf::Ones( lbl[i].size() );
	train( ios, lbl, weight, max_feature );
}
void LearnedSeed::train( const std::vector< ImageOverSegmentation* > &ios, const std::vector<VectorXs> & lbl, int max_feature ) {
	std::vector<VectorXf> weight( lbl.size() );
	for( int i=0; i<lbl.size(); i++ )
		weight[i] = VectorXf::Ones( lbl[i].size() );
	train( ios, lbl, weight, max_feature );
}

class LogisticSeedTrainer: public EnergyFunction {
protected:
	std::vector<int> ids_;
	const std::vector<SeedFeatureVector> & f_;
	const std::vector<VectorXb> & l_;
public:
	LogisticSeedTrainer( const std::vector<SeedFeatureVector> & f, const std::vector<VectorXb> & l, const std::vector<bool> & active ): f_(f), l_(l) {
		for( int i=0; i<(int)active.size(); i++ )
			if( active[i] )
				ids_.push_back( i );
	}
	virtual VectorXf initialGuess() const {
		const int N = f_[0].dim();
		return VectorXf::Zero( N );
	}
	virtual VectorXf gradient(const VectorXf &x, float &e) const {
		const float L2_norm = 1e-5;
		e = 0;
		VectorXf g = 0*x;
#pragma omp parallel for
		for( int i=0; i<(int)ids_.size(); i++ ) {
			int id = ids_[i];
			RMatrixXf fm = f_[id];
			// Subtract the colwise mean (to make things more numerically stable) [might need to do it a few times]
			fm.rowwise() -= fm.colwise().mean();
			fm.rowwise() -= fm.colwise().mean();
			fm.rowwise() -= fm.colwise().mean();
			
			// Compute the response
			VectorXf fx = fm * x;
			float mx = fx.maxCoeff(), spos=0;
			
			// Compute the prob dist
			VectorXf efx = (fx.array()-mx).exp();
			efx.array() /= efx.array().sum()+1e-10;
			for( int j=0; j<efx.size(); j++ )
				if( l_[id][j] )
					spos += efx[j];
			
			// Update the energy
#pragma omp atomic
			e -= spos;
			
			// Update the gradient
			VectorXf gg = ((l_[id].cast<float>().array()-spos)*efx.array()).matrix().transpose()*fm;
#pragma omp critical
			g -= gg;
		}
		e += 0.5*L2_norm*x.squaredNorm();
		g += L2_norm*x;
		return g;
	}
};
class LogLogisticSeedTrainer: public EnergyFunction {
protected:
	std::vector<int> ids_;
	const std::vector<SeedFeatureVector> & f_;
	const std::vector<VectorXb> & l_;
	const std::vector<VectorXf> & w_;
public:
	LogLogisticSeedTrainer( const std::vector<SeedFeatureVector> & f, const std::vector<VectorXb> & l, const std::vector<VectorXf> & w, const std::vector<bool> & active ): f_(f), l_(l), w_(w) {
		for( int i=0; i<(int)active.size(); i++ )
			if( active[i] )
				ids_.push_back( i );
	}
	virtual VectorXf initialGuess() const {
		const int N = f_[0].dim();
		return VectorXf::Zero( N );
	}
	virtual VectorXf gradient(const VectorXf &x, float &e) const {
		const float L2_norm = 0, EPS = 1e-4;
		e = 0;
		double se = 0;
		VectorXf g = 0*x;
#pragma omp parallel for
		for( int i=0; i<(int)ids_.size(); i++ ) {
			int id = ids_[i];
			RMatrixXf fm = f_[id];
			// Subtract the colwise mean (to make things more numerically stable) [might need to do it a few times]
			fm.rowwise() -= fm.colwise().mean();
			fm.rowwise() -= fm.colwise().mean();
			fm.rowwise() -= fm.colwise().mean();
			
			// Compute the response
			VectorXf fx = fm * x;
			
			// Compute the prob dist
			// TODO: Make sure the mx and does not land in w_[id]=0
			float mx = fx.maxCoeff();
			VectorXf efx = w_[id].array()*(fx.array()-mx).exp();
			
			// Compute the positive prob dist
			int n_pos = (l_[id].array() && (w_[id].array()>EPS)).cast<int>().sum();
			VectorXf fx_pos( n_pos ), w_pos( n_pos );
			for( int j=0,k=0; j<efx.size(); j++ )
				if( l_[id][j] && w_[id][j] > EPS ) {
					w_pos[k] = w_[id][j];
					fx_pos[k] = fx[j];
					k++;
				}
			float mx_pos = fx_pos.maxCoeff();
			VectorXf efx_pos = w_pos.array()*(fx_pos.array()-mx_pos).exp();
			// Update the energy
#pragma omp atomic
			se -= log(efx_pos.array().sum())+mx_pos - (log(efx.array().sum())+mx);
			
			// Update the gradient
			efx.array() /= efx.array().sum();
			efx_pos.array() /= efx_pos.array().sum();
			VectorXf gg = VectorXf::Zero( x.size() );
			for( int j=0,k=0; j<efx.size(); j++ ) {
				if( l_[id][j] )
					gg += efx_pos[k++] * fm.row(j);
				gg -= efx[j] * fm.row(j);
			}
			
#pragma omp critical
			g -= gg;
		}
		e = se + 0.5*L2_norm*x.squaredNorm();
		g += L2_norm*x;
		return g;
	}
};
void LearnedSeed::train(std::vector< SeedFeatureVector > &f, const std::vector<VectorXs> & lbl, const std::vector<VectorXf> & weight, int max_seed) {
	// Initialize the labels
	std::vector< VectorXb > l( lbl.size() );
	std::vector< VectorXf > w = weight;
	std::vector< bool > active( lbl.size() );
	std::vector< VectorXb > hit( lbl.size() );
	int n_obj = 0, n_got_tot = 0;
	for( int i=0; i<lbl.size(); i++ ) {
		n_obj += lbl[i].maxCoeff()+1;
		hit[i] = VectorXb::Zero( lbl[i].maxCoeff()+1 );
		l[i] = lbl[i].array()>=0;
		active[i] = l[i].any();
	}
	for( int it=0; it<max_seed; it++ ) {
		LogLogisticSeedTrainer trainer( f, l, w, active );
		// TODO: Maybe a few random initiaziations would work better
		float e;
		VectorXf new_w = minimizeLBFGS( trainer, e, false );
		
		int n_got = 0, n_active = 0;
		// Update the seeds and features
		for( int i=0; i<lbl.size(); i++ ) 
			if( active[i] ){
				RMatrixXf fm = f[i];
				// Subtract the colwise mean (to make things more numerically stable)
				fm.rowwise() -= fm.colwise().mean();
				fm.rowwise() -= fm.colwise().mean();
				fm.rowwise() -= fm.colwise().mean();
				
				// Compute the response
				VectorXf fw = fm*new_w;
				int mx=0;
				fw.maxCoeff( &mx );
				
				// Update the feature
				f[i].update( mx );
				
				n_active += active[i];
				// Remove the current segment (or part of it)
				if( l[i][mx] ) {
// 					float ww = w[i][mx]+0.25;
					for( int k=0; k<l[i].size(); k++ )
						if(l[i][k] && lbl[i][k] == lbl[i][mx]) {
							// Update the weight and label
// 							w[i][k] -= ww;
// 							if( w[i][k] < 1e-4 ) {
								l[i][k] = 0;
								w[i][k] = 1;
// 							}
						}
					active[i] = l[i].any();
					if(!hit[i][ lbl[i][mx] ])
						n_got++;
					hit[i][ lbl[i][mx] ] = 1;
				}
		}
		
		// Add the learned feature
		w_.push_back( new_w );
		n_got_tot += n_got;
		printf("[%3d]  Total objects %5.1f%%       Got %4d / %4d objects\r", it, 100.*n_got_tot/n_obj, n_got,n_active );
		fflush( stdout );
	}
	printf("Training got %d / %d objects: %0.1f%%                            \n", n_got_tot, n_obj, 100.*n_got_tot/n_obj );
}
VectorXi LearnedSeed::computeImageOverSegmentation( const ImageOverSegmentation & ios, int M ) const {
	SeedFeatureVector f = features_.create( ios );
	
	VectorXb used = VectorXb::Zero( ios.Ns() );
	std::vector<int> r;
	for( int i=0; i<M && i<w_.size(); i++ ) {
		RMatrixXf fm = f;
		fm.rowwise() -= fm.colwise().mean();
		fm.rowwise() -= fm.colwise().mean();
		fm.rowwise() -= fm.colwise().mean();
		
		int mx=0;
		VectorXf fw = fm*w_[i];
		fw.maxCoeff( &mx );
		if( !used[mx] ) {
			used[mx] = 1;
			r.push_back( mx );
		}
		
		f.update( mx );
	}
	return VectorXi::Map( r.data(), r.size() );
}
void LearnedSeed::load(std::istream &s) {
	int M;
	s.read((char*)&M,sizeof(M));
	w_.resize(M);
	for( int i=0; i<M; i++ )
		loadMatrixX( s, w_[i] );
}
void LearnedSeed::save(std::ostream &s) const {
	int M = w_.size();
	s.write((const char*)&M,sizeof(M));
	for( int i=0; i<M; i++ )
		saveMatrixX( s, w_[i] );
}
VectorXi SaliencySeed::computeImageOverSegmentation( const ImageOverSegmentation & ios, int M) const {
	std::mt19937 rand(0);
	VectorXf sal = Saliency().saliency( ios.image(), ios.s() );
	float ns = sal.array().sum();
	VectorXi r( M );
	for( int i=0; i<M; i++ ) {
		float v = rand()*ns / RAND_MAX;
		int j=0;
		for(j=0;j<sal.size() && v>sal[j]; j++ )
			v -= sal[j];
		r[i] = j;
	}
	return r;
}
void LearnedSeed::load(const std::string &s) {
	std::ifstream is(s);
	if(!is.is_open())
		throw std::invalid_argument( "Could not open file '"+s+"'!" );
	load(is);
}
void LearnedSeed::save(const std::string &s) const {
	std::ofstream os(s);
	if(!os.is_open())
		throw std::invalid_argument( "Could not write file '"+s+"'!" );
	save(os);
}
