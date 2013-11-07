// A simple mex wrapper around  the graph-cut code of Yuri Boykov (yuri@csd.uwo.ca) and Vladimir Kolmogorov (v.kolmogorov@cs.ucl.ac.uk)
//
// by Ian Endres (iendres2@uiuc.edu)
#include <stdio.h>
#include <stdint.h>
#include "mex.h"
#include "../graph.h"
#include "../graph.cpp"
#include "../maxflow.cpp"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{
   //mex_maxflow(unary_cost, pw_cost, pw_edge)

   if(nrhs<3) {
      mexErrMsgTxt("mex_maxflow(unary_cost, pw_cost, pw_edge)");
   }

   double *unary = (double *)mxGetData(prhs[0]);
   double *pw = (double *)mxGetData(prhs[1]);
   double *edges = (double *)mxGetData(prhs[2]);

   const int *unary_dim = mxGetDimensions(prhs[0]);
   const int num_nodes = unary_dim[0];
   const int *pw_dim = mxGetDimensions(prhs[1]);
   const int num_edges = pw_dim[0];


//   printf("nodes: %d, edges: %d\n", num_nodes, num_edges);



   typedef Graph<int64_t,int64_t,int64_t> LongLongGraph;
   LongLongGraph g(num_nodes, num_edges); 

   int *node_list = new int[num_nodes];
   

   for(int i=0; i<num_nodes; i++) {
//      printf("Adding node %d\n", i);
      node_list[i] = g.add_node();
      g.add_tweights(node_list[i], (int64_t)unary[i], (int64_t)unary[num_nodes+i]);
   }


   for(int i=0; i<num_edges; i++) {
//      printf("Adding edge %d, %d, %d\n", i, (int)edges[i]-1, node_list[(int)edges[i]-1]);
      g.add_edge(node_list[(int)edges[i]-1], node_list[(int)edges[i+num_edges]-1], (int64_t)pw[i], (int64_t)pw[i]);
   }
//   printf("Computing flow\n");
   double flow_val = (double)g.maxflow();

   if(nlhs>0) {
      plhs[0] = mxCreateNumericMatrix(num_nodes, 1, mxDOUBLE_CLASS, mxREAL);
      double *cut_assignment = (double *)mxGetPr(plhs[0]);

      for(int i=0; i<num_nodes; i++) {
         cut_assignment[i] = (double)g.what_segment(node_list[i]);
      }
   }

   if(nlhs>1) {
      plhs[1] = mxCreateDoubleScalar((double)flow_val);
   }

   delete [] node_list;
   node_list = 0;
}
