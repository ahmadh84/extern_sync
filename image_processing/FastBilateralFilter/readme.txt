This is the MATLAB implementation of the fast O(1) bilateral filter described in the following papers:

[1] K.N. Chaudhury, D. Sage, and M. Unser, "Fast O(1) bilateral filtering using trigonometric range kernels," IEEE Trans. Image 
Processing, vol. 20, no. 11, 2011.

[2] K.N. Chaudhury, "Acceleration of the shiftable O(1) algorithm for bilateral filtering and non-local means," 
arXiv:1203.5128v1. 

Also included is the MATLAB implementation of joint bilateral filtering of multiband images.

Author     : Kunal N. Chaudhury (kchaudhu@math.princeton.edu)
Date       : March 2012

To run the software use:
=======================


[ outImg , param ]  =  shiftableBF(inImg, sigma1, sigma2, w, tol)


INPUT
=====

inImg      : grayscale image
sigma1     : width of spatial Gaussian
sigma2     : width of range Gaussian
[-w, w]^2  : domain of spatial Gaussian (note: w must be odd)
tol        : truncation parameter 


OUTPUT
======

outImg     : filtered image
param      : lsit of parameters


 
