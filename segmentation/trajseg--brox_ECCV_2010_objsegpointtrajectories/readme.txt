This program is for scientific use only. Any commercial use of this
package or parts of it is prohibited.
________________________________________________________________

Motion segmentation binary for 64 bit Linux
________________________________________________________________

(c) Thomas Brox 2010

If you use this program, you should cite the following paper:

T. Brox, J. Malik: Object segmentation by long term analysis of 
point trajectories, Proc. European Conference on Computer Vision, 
2010. 

---------------

Usage: 

./motionsegBM bmfFile startFrame numberOfFrames sampling *numberOfObjects 
              *affineMergingFlag

bmfFile is a text file with a very short header, comprising the 
number of images in the sequence and 1. After the header all 
image files of the sequences are listed separated by line breaks. 
See marple2.bmf for an example. All input files must be in the 
PPM format (P6). 

startFrame is the frame where the computation is started. 
Usually this is frame 0. If you want to start later in the 
sequence you may specify another value.

numberOfFrames is the number of frames for which you want to 
run the computation. Make sure that the value is not larger than 
the total number of frames of your sequence.

sampling specifies the subsampling parameter. If you specify 8 (a 
good default value), only every 8th pixel in x and y direction is 
taken into account. If you specify 1, the sampling will be dense 
(be careful, memory consumption and computation time will be very 
large in this setting). 

*numberOfObjects is an optinal parameter that specifies the number of 
objects assumed for the segmentation. If you do not specify this number, 
the number of objects will be determined automatically, as originally 
proposed in the paper. The original option has the value -1.
Note that specifying the number of objects can still result in fewer
clusters, due to empty clusters during the optimization.

*affineMergingFlag is a flag that allows to turn off the affine merging 
post-processing step. Enter 0 to turn this merging off. Default value 
is 1 as originally proposed in the paper. This parameter must be given 
as the 6th argument; specify -1 for *numberOfObjects to use this flag in
conjunction with the automatic model selection.

The output can be found in a subdirectory of the directory where 
the input sequence is stored. It comprises a text file 
TracksNumberOfFrames.dat with all the tracks and their labels. 
For further details how to interpret the text file have a look at
readWriteTracks.cpp. Additionally, the directory contains the 
computed optical flow fields used for tracking (in the Middlebury 
flo format), as well as a visualization of the tracked points 
Tracking???.ppm and their labels Segments???.ppm. 

You can ignore error messages reporting a singular matrix. This
is the normal outcome of a testing procedure. 

Linking problems: 

motionsegBM requires a dynamic library called libg2c.so.0 that
is no longer installed by default on all systems. The current 
package includes libg2c.so.0, so in case of trouble running 
motionsegBM, your system administrator can install it for you 
or you can set the path to this file via 
export LD_LIBRARY_PATH=<Directory>:$LD_LIBRARY_PATH 

_________________________________________________________________

Example
__________________________________________________________________

Example 1:
./conv2ppm
./motionsegBM marple2/marple2.bmf 0 50 8

This converts all JPEG images of the example sequence to PPMs
and runs the motion segmentation program on the first 50 images
of this sequence. 

Example 2:
./conv2ppm
./motionsegBM marple2/marple2.bmf 0 50 8 -1 0

This converts all JPEG images of the example sequence to PPMs
and runs the motion segmentation program on the first 50 images
of this sequence with the proposed automatic model selection, but
without the affine merging step.

_________________________________________________________________

Bugs
__________________________________________________________________

10.2.2012: There has been a typo in the file marple2Def.dat, which
           has been corrected. This does not affect the errors 
           reported by the evaluation tool.

8.2.2011: There has been a minor mistake in the file car4Def.dat, 
          which has been corrected. This affects the errors reported
          by the evaluation tool. The corrected version usually 
          reports smaller errors.   


Please report bugs to brox@informatik.uni-freiburg.de

