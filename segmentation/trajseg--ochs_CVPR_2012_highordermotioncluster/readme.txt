This program is for scientific use only. Any commercial use of this
package or parts of it is prohibited.
________________________________________________________________

Higher Order Motion Segmentation binary for 64 bit Linux
________________________________________________________________

(c) Peter Ochs 2012

If you use this program, you should cite the following paper:

P. Ochs, T. Brox: Higher Order Motion Models and Spectral Clustering,
IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2012.

---------------

Usage: 

./motionsegOB bmfFile startFrame numberOfFrames sampling

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


./conv2ppm
./motionsegOB marple2/marple2.bmf 0 50 8

This converts all JPEG images of the example sequence to PPMs
and runs the motion segmentation program on the first 50 images
of this sequence. 

_________________________________________________________________

Bugs
__________________________________________________________________

Please report bugs to brox@informatik.uni-freiburg.de

