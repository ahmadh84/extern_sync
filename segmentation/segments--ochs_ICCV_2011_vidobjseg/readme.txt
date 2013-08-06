**************************************************************************************
*                                                                                    *
* Linux 64bit binaries for turning a sparse labeling into dense regions              *
*                                                                                    *
* Copyright (c) 2011 Peter Ochs                                                      *
*                                                                                    *
* ------------                                                                       *
* Terms of use                                                                       *
* ------------                                                                       *
* This program is provided for research purposes only. Any commercial                *
* use is prohibited. If you are interested in a commercial use, please               *
* contact the copyright holder.                                                      * 
*                                                                                    *
* If you used this program in your research work, you should cite the                *
* following publication:                                                             *
*                                                                                    *
* P. Ochs and T. Brox.                                                               *
* Object Segmentation in Video: A Hierarchical Variational Approach                  *
*   for Turning Point Trajectories into Dense Regions,                               *
* IEEE International Conference on Computer Vision (ICCV), Barcelona, November 2011. *
*                                                                                    *
* This program is distributed WITHOUT ANY WARRANTY.                                  *
*                                                                                    *
**************************************************************************************

-----------------
Table of Contents
-----------------
1) How to use the binary
2) How to compute the superpixels
3) Converter: track file -> sparse labels
4) Evaluation
5) References
6) Bug Report (version 09.Jan 2012)

------------------------
1) How to use the binary
------------------------
The binary file "dens100" provides the algorithm from [1] to turn a sparse 
segmentation into dense regions. It is executed by the command 
  "./dens100 <arg1.ppm> <arg2> <arg3.ppm> <arg4*> <arg5*>"
The arguments are explained below:
<arg1.ppm> is the original input image. [e.g. marple1_001.ppm]
<arg2>     is the (non-negative integer) number of hierarchy levels (0,1,2,...,<arg2>-1), 
           ordered from finer to coarser superpixel levels.
           (0 -> single level model, else -> multi-level model)
           The superpixel levels are expected to be given as piecewise constant images
           named as follows:
           "<arg1>_seg0.ppm", "<arg1>_seg1.ppm", ... .
           A piecewise constant region on a coarser superpixel level 'i' must be the union of 
           piecewise constant regions on a finer superpixel level 'j', 'j'<'i'. 
           (See (2) on how to get the superpixel).
<arg3.ppm> is the image with sparse labels. Each color identifies a different label. 
           White pixels symbolize no label information.
<arg4*>    (optional) is the regularization weight between data term and regularization term, 
           which is denoted by 'alpha' in [1]. It is a real number in ]0,1[. The smaller 
           'alpha' the more preference to smoothing. [default: alpha=0.3]
<arg5*>    (optional) is the output style. [default = 0]
           (0: image with dense regions,
            1: additional image with dense regions overlayed to the original image,
            2: as '1' with initial labels)

---------------------------------
2) How to compute the superpixels
---------------------------------
We used the source code (for Linux/Mac, 32/64 bits) from 

  http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html#algorithms 

which implements [3] to compute the superpixels. In order to create piecewise constant images, 
we used the file "createSuperpixel.m" to execute their superpixel algorithm, which creates 3 
superpixel levels thresholding the boundary probability set at 0.05, 0.1, and 0.15. 
It takes 2 input arguments:
<arg1>     image folder
<arg2>     image name (".ppm" is added automatically).


-------------------
3) Additional tools
-------------------
The provided binary "conv_Track2Img" converts a track file as obtained from [2] to 
a sparse label image which can be used as input (<arg3.ppm>) for "dens100". 
It takes 3 arguments:
<arg1.ppm> is an arbitrary original image from the desired video sequence. 
           [e.g. marple1_001.ppm] It is used only to determine the image size.
<arg2.dat> is a track file as obtained from [2].
<arg3>     is a (non-negative integer) number giving the frame number in the video 
           sequence. [e.g. 0 to obtain the sparse labels for the first image 
           in the video sequence]


-------------
4) Evaluation 
-------------
In order to evaluate your results convert your result images into a track 
file for each sequence and use the evaluation tool provided in [2].


-------------
5) References 
-------------
[1] P. Ochs and T. Brox
Object Segmentation in Video: A Hierarchical Variational Approach 
  for Turning Point Trajectories into Dense Regions,
International Conference on Computer Vision (ICCV), Barcelona, Spain, Nov. 2011.
http://lmb.informatik.uni-freiburg.de/people/ochs/publication_all.html

[2] T. Brox, J. Malik
Object segmentation by long term analysis of point trajectories,
European Conference on Computer Vision (ECCV), Crete, Greece, Springer, LNCS, Sept. 2010.
http://lmb.informatik.uni-freiburg.de/people/brox/publication_selected.html
Download Benchmark from: http://lmb.informatik.uni-freiburg.de/resources/datasets/

[3] P. Arbelaez, M. Maire, C. Fowlkes, and J. Malik. 
Contour detection and hierarchical image segmentation. 
IEEE Transactions on Pattern Analysis and Machine Intelligence (TPAMI), 2011.
http://www.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/index2.html


-------------
6) Bug report 
-------------
- 09. Feb 2012:
	* Executables for evaluation removed. For the intermediate way of evaluation see 4). 
		New evaluation executable which converts the images to a track file will come very 
		soon.
- 09. Jan 2012:
  * Implementation of projection onto l1-unit ball changed. (Faster now).
  * Non-local neighbors implementation added. (Missing before).
  * Added possibility to start the programm with given track-file from [2].
- 20. Dez 2011:
  * Output label color adapted to input label color.
