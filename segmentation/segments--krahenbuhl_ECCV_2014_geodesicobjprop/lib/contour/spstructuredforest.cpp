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
#include "spstructuredforest.h"
#include "util/util.h"
#include <unordered_map>
#include <iostream>

#define FULL_PATCH
static RMatrixXu8 toPatch( const VectorXus & bnd_data, int PS ) {
	RMatrixXu8 bnd = RMatrixXu8::Zero(PS,PS);
	for( int i: bnd_data )
		bnd.data()[i] = 1;
	RMatrixXu8 r = RMatrixXu8::Zero(PS,PS);
	int id = 1;
	for( int j=0; j<PS; j++ )
		for( int i=0; i<PS; i++ )
			if (!bnd(j,i) && !r(j,i)){
				// Run DFS
				std::vector<Vector2i> q;
				q.push_back({i,j});
				while(!q.empty()) {
					VectorXi p = q.back();
					q.pop_back();
#ifdef FULL_PATCH
					if( bnd(p[1], p[0]) && !r(p[1], p[0]) )
						r(p[1], p[0]) = id;
#endif
					if( bnd(p[1], p[0]) || r(p[1], p[0]) )
						continue;
					r(p[1], p[0]) = id;
					if( 0<p[0]    && r(p[1],p[0]-1)==0 ) q.push_back({p[0]-1,p[1]});
					if( 0<p[1]    && r(p[1]-1,p[0])==0 ) q.push_back({p[0],p[1]-1});
					if( p[0]+1<PS && r(p[1],p[0]+1)==0 ) q.push_back({p[0]+1,p[1]});
					if( p[1]+1<PS && r(p[1]+1,p[0])==0 ) q.push_back({p[0],p[1]+1});
				}
				id++;
			}
			
#ifdef FULL_PATCH
	for( int j=0; j<PS; j++ )
		for( int i=0; i<PS; i++ ) {
			if( i && r(j,i) == 0 )
				r(j,i) = r(j,i-1);
			else if( j && r(j,i) == 0 )
				r(j,i) = r(j-1,i);
		}
	for( int j=PS-1; j>=0; j-- )
		for( int i=PS-1; i>=0; i-- ) {
			if( i && r(j,i-1) == 0 )
				r(j,i-1) = r(j,i);
			else if( j && r(j-1,i) == 0 )
				r(j-1,i) = r(j,i);
		}
	const int n_bnd = (r.array() == 0).cast<int>().sum();
	if( n_bnd > 0 )
		printf("#0 = %d\n", n_bnd );
#endif
	return r;
}

SPStructuredForest::SPStructuredForest(const StructuredForest &o):settings_(o.settings_) {
	const int PS = settings_.out_patch_size;
	std::unordered_map< std::string, int > phash;
	std::vector<RMatrixXu8> patches;
	patches.push_back( RMatrixXu8::Ones( PS, PS ) );
	for ( auto t: o.forest_.trees() ) {
		LabelTree lt;
		lt.set( t );
		
		// Translate the labels (and compress)
		for( int i=0; i<t.data().size(); i++ ) {
			RangeData d = t.data()[i];
			if( d.begin < d.end ) {
				VectorXus p = o.patch_ids_.segment( d.begin, d.end-d.begin );
				std::sort( p.data(), p.data()+p.size() );
				std::string hp( (char*)p.data(), sizeof(p[0])*p.size() );
				if( phash.count( hp ) )
					lt.data()[i] = phash[ hp ];
				else {
					lt.data()[i] = phash[ hp ] = patches.size();
					patches.push_back( toPatch( p, PS ) );
				}
			}
			else
				lt.data()[i] = 0;
		}
		
		// Add the new tree
		forest_.addTree( lt );
	}
	// Store the data
	patches_ = RMatrixXu8::Zero( patches.size(), PS*PS );
	n_patch_lbl_ = VectorXu8::Zero( patches.size() );
	for( int i=0; i<(int)patches.size(); i++ ) {
		std::copy( patches[i].data(), patches[i].data()+PS*PS, patches_.data()+i*patches_.cols() );
		n_patch_lbl_[i] = patches[i].maxCoeff();
	}
}
void SPStructuredForest::load(const std::string &fn) {
	std::ifstream s( fn );
	if(!s.is_open())
		throw std::invalid_argument( "Could not open file '"+fn+"'!" );
	forest_.load( s );
	// Load patch_ids
	int sz[2];
	s.read( (char*)sz, sizeof(sz) );
	patches_ = RMatrixXu8(sz[0],sz[1]);
	s.read( (char*)patches_.data(), sizeof(unsigned short)*sz[0]*sz[1] );
	s.read( (char*)&settings_, sizeof(settings_) );
}
void SPStructuredForest::save(const std::string &fn) const {
	std::ofstream s( fn );
	forest_.save( s );
	// Save patch_ids
	int sz[2] = {(int)patches_.rows(),(int)patches_.cols()};
	s.write( (const char*)sz, sizeof(sz) );
	s.write( (const char*)patches_.data(), sizeof(unsigned char)*sz[0]*sz[1] );
	s.write( (const char*)&settings_, sizeof(settings_) );
}
RMatrixXf SPStructuredForest::computeHist(const SFFeatures &f, const RMatrixXs &s, const Edges &edges, int N_EHIST_BIN, bool per_tree ) const {
	const int Ns = getN( edges ), W = s.cols(), H = s.rows(), PS = settings_.out_patch_size;
	const int Rp = PS/2, pW = W+PS, pH = H+PS;
	eassert( PS == 2*Rp );
	
	std::vector<int> os;
	for( int j=0; j<PS; j++ )
		for( int i=0; i<PS; i++ )
			os.push_back( i+j*pW );
	
	std::unordered_map<Edge, int> edge_map;
	for( int i=0; i<(int)edges.size(); i++ )
		edge_map[ edges[i] ] = i;
	
	// Pad the segmentation
	RMatrixXs pad_s = -RMatrixXs::Ones( pH, pW );
	pad_s.block( Rp, Rp, H, W ) = s;
	
	// Upsample all edges
	RMatrixXf ehist = RMatrixXf::Zero( edges.size(), (per_tree?forest_.nTrees()/2:1)*N_EHIST_BIN );
	
	// Offset for the segment lookup
	const int ND=13;
	const int OX[] = { 0, 1,-1, 0, 0, 1,-1, 1,-1, 2,-2, 0, 0, 2,-2, 2,-2};
	const int OY[] = { 0, 0, 0, 1,-1, 1, 1,-1,-1, 0, 0, 2,-2, 2, 2,-2,-2};
	int c[ND] = {};
	bool d[ND] = {};
	
	const int stride = f.x_[1] - f.x_[0], max_lbl = n_patch_lbl_.maxCoeff(), max_seg = ND+1;
	std::vector<int> id_remap( Ns, -1 );
	std::vector<float> hist( max_seg*max_lbl, 0 ), scnt( max_seg, 0 );
	std::vector<int> segs( max_seg );
	for( int t=0; t<forest_.nTrees(); t++ ) {
		const int bin_o = (per_tree?t/2:0)*N_EHIST_BIN;
		for( int i=0; i<f.id_.size(); i++ ) {
			const int x = f.x_[i], y = f.y_[i];
			if( (((x+y)/stride)&1) == (t&1) ) {
				const int cx = x+Rp, cy = y+Rp;
				// Are we close to an edge
				bool any = false;
				for( int k=0; k<ND; k++ ) {
					c[k] = pad_s(cy+OY[k],cx+OX[k]);
					d[k] = (c[k]>=0 && (!k || c[k]!=c[0]));
					if(k)
						any = any || d[k];
				}
				if( c[0]>=0 && any ) {
					// Evaluate the forest
					int id = forest_.tree(t).predictData( f, i );
					
					// Find a (unique) list of all adjacent segments and build the id_map
					int nid = 0;
					for( int k=0; k<ND; k++ )
						if( d[k] && id_remap[ c[k] ]==-1 ) {
							segs[nid] = c[k];
							id_remap[ c[k] ] = nid++;
						}
					// Project the patch onto the superpixels
					const short * ps = pad_s.data() + (x + y*pW);
					const uint8_t * ppatch = patches_.data()+id*patches_.cols();
					const int n_lbl = n_patch_lbl_[id];
					
					// NOTE: Looping over all PS*PS is a bit slow [~60-80ms on peppers.png]
					for( int k=0; k<PS*PS; k++ ) {
						const int sg = ps[os[k]], pid = ppatch[k]-1;
						if( sg>=0 && id_remap[sg]>=0 && pid>=0 ) {
							// Remap the id
							int rid = id_remap[ sg ];
							// Add an element to the histogram for segment ps[os[k]] and edge ppatch[k]
							hist[ rid*n_lbl + pid ]++;
							scnt[ rid ]++;
						}
					}
					for( int k1=0; k1+1<nid; k1++ )
						for( int k2=k1+1; k2<nid; k2++ ) {
							const int s1 = segs[k1], s2 = segs[k2];
							// For each edge [prject onto it]
							const int a = id_remap[ s1 ], b = id_remap[ s2 ];
							if( scnt[ a ] > 0 && scnt[ b ] > 0 ) {
								auto it = edge_map.find( Edge( s1, s2 ) );
								if( it != edge_map.end() ) {
									int eid = it->second;
									// Edge_strength = 1 - sum_l h1(l) h2(l) [or hist intersection kernel]
									float es = 1;
									for( int k=0; k<n_lbl; k++ )
// 										es -= hist[ a*n_lbl + k ] / scnt[ a ] * hist[ b*n_lbl + k ] / scnt[ b ];
										es -= std::min( hist[ a*n_lbl + k ] / scnt[ a ], hist[ b*n_lbl + k ] / scnt[ b ] );
									int bin = (int)(std::min(0.99999f, std::max(0.0f, es))*N_EHIST_BIN);
									ehist( eid, bin_o + bin ) += 1;
								}
							}
							else
								printf("Error a histgram count is 0: %d %d\n",(int)scnt[a], (int)scnt[b]);
						}
					// Clear the idmap and histogram
					memset( hist.data(), 0, nid*n_lbl*sizeof(float) );
					memset( scnt.data(), 0, nid*sizeof(float) );
					for( int k=0; k<nid; k++ )
						id_remap[ segs[k] ] = -1;
				}
			}
		}
	}
	return ehist;
}
VectorXf SPStructuredForest::detect(const Image8u &rgb_im, const RMatrixXs &s, const Edges &edges) const {
	const int N_EHIST_BIN=21;
	if( rgb_im.W() != s.cols() || rgb_im.H() != s.rows() )
		throw std::invalid_argument("Image and segmentation size do no match!");
	SFFeatures f( rgb_im );
	RMatrixXf ehist = computeHist( f, s, edges, N_EHIST_BIN );
	// Evaluate the histogram
	VectorXf r = VectorXf::Zero( edges.size() );
	// NOTE: This can be a bit slow (10ms)
	for( int i=0; i<edges.size(); i++ ) {
// 		// Pick the 75th percentile
// 		float p = 0.75;
// 		r[i] = 0;
// 		VectorXf h = ehist.row(i);
// 		h.array() /= h.array().sum()+1e-10;
// 		for( int k=0; k<=N_EHIST_BIN; k++ ) {
// 			if( p <= 0 ) {
// 				r[i] = 1.0*k/N_EHIST_BIN;
// 				break;
// 			}
// 			if( k < N_EHIST_BIN )
// 				p -= h[k];
// 		}
		// Pick the mean
		r[i] = 0;
		VectorXf h = 1*ehist.row(i);
		h.array() /= h.array().sum()+1e-10;
		for( int k=0; k<N_EHIST_BIN; k++ ) {
			r[i] += h[k]*k/(N_EHIST_BIN-1);
		}
	}
	return r.array();
}
RMatrixXf SPStructuredForest::computeHistPerTree(const Image8u &im, const RMatrixXs &s, const Edges &edges, int n_bin, bool normalize ) const {
	if( im.W() != s.cols() || im.H() != s.rows() )
		throw std::invalid_argument("Image and segmentation size do no match!");
	RMatrixXf r( edges.size(), n_bin*forest_.nTrees()/2 );
	SFFeatures f( im );
	RMatrixXf h = computeHist( f, s, edges, n_bin, true );
	if( normalize ) {
		// Normalize each block of data (belonging to a tree)
		for( int i=0; i<forest_.nTrees()/2; i++ )
			h.array().block(0,i*n_bin,h.rows(),n_bin).colwise() /= h.array().block(0,i*n_bin,h.rows(),n_bin).rowwise().sum();
	}
	return h;
}
