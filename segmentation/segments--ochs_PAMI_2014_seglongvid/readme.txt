Linux 64bit binaries for running the motion segmentation as in the paper


==================================================================
||  Segmentation of moving objects by long term video analysis  ||
==================================================================


Copyright (c) 2013 Peter Ochs

__________________________________________________________________

Terms of use
__________________________________________________________________

This program is provided for research purposes only. Any commercial
use is prohibited. If you are interested in a commercial use, please 
contact the copyright holder. 

If you used this program in your research work, you should cite the 
following publication:

Peter Ochs, Jitendra Malik, Thomas Brox 
Segmentation of moving objects by long term video analysis
IEEE Transactions on Pattern Analysis and Machine Intelligence, preprint, 2013. 

This program is distributed WITHOUT ANY WARRANTY.

__________________________________________________________________

Usage MoSeg binary
__________________________________________________________________

Usage: 
------

./MoSeg filestructure.cfg sequenceName startFrame numberOfFrames 
        sampling regWeightSegmentation

filestructure.cfg   is a document containing the information about
input and output directories. See the example files for its usage.

sequenceName   is only the name of the sequence. This name is taken 
as foldername relative to the input path specified in the 
filestructure. It assumes that the required bmfFile is located in 
the directory with the same name. 
(bmfFile is a text file with a very short header, comprising the 
number of images in the sequence and 1. After the header all 
image files of the sequences are listed separated by line breaks. 
See marple8.bmf for an example. All input files must be in the 
PPM format (P6).) 

startFrame   is the frame where the computation is started. 
Usually this is frame 0. If you want to start later in the 
sequence you may specify another value.

numberOfFrames   is the number of frames for which you want to 
run the computation. Make sure that the value is not larger than 
the total number of frames of your sequence.

sampling   specifies the subsampling parameter. If you specify 8 (a 
good default value), only every 8th pixel in x and y direction is 
taken into account. If you specify 1, the sampling will be dense 
(be careful, memory consumption and computation time will be very 
large in this setting). 

regWeightSegmentation   is the weight used in the segmentation of the 
trajectories with spatial regularization. Please refer to the paper for 
details. The default value is 60.

The output comprises a text file TracksNumberOfFrames.dat with all 
the tracks and their labels. For further details how to interpret 
the text file have a look at readWriteTracks.cpp. Additionally, 
a visualization of the tracked points Tracking???.ppm and their labels 
Segments???.ppm is given. The computed optical flow fields used for 
tracking (in the Middlebury flo format) are stored in a subfolder 
of the Data-folder.

You can ignore error messages reporting a singular matrix. This
is the normal outcome of a testing procedure. 


Linking problems: 
-----------------

MoSeg requires a dynamic library called libg2c.so.0 that
is no longer installed by default on all systems. The current 
package includes libg2c.so.0, so in case of trouble running 
MoSeg, your system administrator can install it for you 
or you can set the path to this file via 
export LD_LIBRARY_PATH=<Directory>:$LD_LIBRARY_PATH 


Example:
--------

./prepareDemoSequence.sh
./MoSeg filestructureTrainingSet.cfg marple8 0 10 8 

This converts all JPEG images of the example sequence to PPMs
and runs the motion segmentation program on the first 10 images
of this sequence. 

Run on all sequences:
---------------------

Use prepareTrainingSequences.sh or prepareTestSequences.sh to 
convert all jpg-files to PPM format (P6).
Use exec_training_all or exec_test_all to run the provided 
motion segmentation binary on all sequences of the training set 
or the test set. Please refer to these files for the required arguments.

__________________________________________________________________

Usage MoSeg-densify binary (requires CUDA 5.5, compute capability 2.0)
__________________________________________________________________

Usage: 
------

./dens100gpu filestructure.cfg imageName trackFile frameIndex resultDir

(Required arguments are shown if less than 4 arguments are passed.)

filestructure.cfg   is a document containing the information about
input and output directories and a few parameters. See the example 
files for its usage.

imageName   is the name of the image relative to the data directory 
specified in the filestructure. If a folder is to be processed 
(indicated by frameIndex=-1) 'relative_path/image.ppm', where 
'image' should NOT be replaced by the imagename, should be given.

trackFile   is the track file relative to the tracks directory specified
in the filestructure.

frameIndex   is index of the considered frame (starting with 0) as 
written in the track file. Giving -1 processes the complete folder 
where imageName refers to. 

resultDir   is a path relative to the resultDir specified in the 
filestructure where the results will be stored.

The output comprises a track-file for the dense results, the densified 
result as a segmentation image, and the segmentation overlayed to the 
original image.


Example:
--------

./dens100gpu filestructureTrainingSetDensify.cfg marple8/image.ppm 
             OchsBroxMalik8_0000_0010_0000060.00/marple8/Tracks10.dat -1 
             OchsBroxMalik8_0000_0010_0000060.00/marple8/DenseSegmentation

This runs the densify-code on the first 10 frames of the demo sequence
marple8. Before you can use this, the motion segmentation code must 
have been run. 

Run on all sequences:
---------------------

Use exec_densify_training_all or exec_densify_test_all to run the provided 
densify binary on all sequences of the training set or the test set. 
Please refer to these files for the required arguments.

__________________________________________________________________

Usage Evaluation binaries for the 
Freiburg Berkeley Motion Segmentation Dataset (FBMS)
__________________________________________________________________

YOU MUST RUN YOUR METHOD ON ALL SEQUENCES USING EXACTLY THE SAME
CODE AND PARAMETERS. 

Usage: 
------

./MoSegEvalAll allShots.txt all allTracksFileList.txt threshold


allShots.txt 
is a txt file containing the number of sequences 
to be evaluated and the references to the files defining the 
evaluation for a sequence. 
See allShotsTrainingSet.txt or allShotsTestSet.txt for details.

all 
specifies that the complete sequences are considered.

allTracksFileList.txt
is a txt file containing the references to the track-files 
generated by your method. 
See eval_test_all or eval_training_all for an automatic 
generation and evaluation for the training and test set.

threshold (default: 0.75)
specifies a value for an object specific F-measure to be 
accepted and counted as an object, i.e., objects with 
F-measure higher than threshold increase the number of object counter.
For your final evaluation you should use the default value 0.75.


The evaluation code provides a result file for each sequence and
evaluation length in the directory of your trajectory file. 
Additionally it provides a summary for all sequences in the main
directory. Refer to the paper to learn how to interpret the 
numbers. The evaluation source code is also provided as a 
reference.


Evaluation on all sequences:
----------------------------
You can use eval_training_all or eval_test_all to evaluate your 
method on the full training or test set. Please refer to these 
files for the required arguments.

_________________________________________________________________

Links
__________________________________________________________________

The full benchmark dataset can be downloaded from: 

http://lmb.informatik.uni-freiburg.de/resources/datasets/FBMS_Trainingset.tar.gz 
http://lmb.informatik.uni-freiburg.de/resources/datasets/FBMS_Testset.tar.gz


The evaluation code for the evaluation binaries provided in this package 
can be downloaded from here:

http://lmb.informatik.uni-freiburg.de/resources/datasets/fbms/pami2013_MoSeg_eval.tar.gz


The source code for computing the flow-variation used in the 
definition of the affinities can be found here:

http://lmb.informatik.uni-freiburg.de/resources/binaries/pami2013_FlowVar.tar.gz

__________________________________________________________________

Bugs
__________________________________________________________________

Please report any bugs to Peter Ochs:    ochs@cs.uni-freiburg.de

