CC=g++_48

OS_NAME=$(shell uname -s)
ifeq ($(OS_NAME),Linux)
  #LAPACKLDFLAGS=/usr/lib64/atlas/libsatlas.so   # single-threaded blas
  LAPACKLDFLAGS=/usr/lib64/atlas/libf77blas.so.3  # multi-threaded blas
endif
ifeq ($(OS_NAME),Darwin)  # Mac OS X
  LAPACKLDFLAGS=-framework Accelerate # for OS X
endif
LAPACKCFLAGS=-Dinteger=int
STATICLAPACKLDFLAGS=-fPIC -Wall -g -fopenmp -static -static-libstdc++ -I/usr/local/include/opencv -I/usr/local/include  /usr/local/lib/libopencv_contrib.a /usr/local/lib/libopencv_stitching.a /usr/local/lib/libopencv_nonfree.a /usr/local/lib/libopencv_superres.a /usr/local/lib/libopencv_ocl.a /usr/local/lib/libopencv_ts.a /usr/local/lib/libopencv_videostab.a /usr/local/lib/libopencv_gpu.a /usr/local/lib/libopencv_photo.a /usr/local/lib/libopencv_objdetect.a /usr/local/lib/libopencv_legacy.a /usr/local/lib/libopencv_video.a /usr/local/lib/libopencv_ml.a /usr/local/lib/libopencv_calib3d.a /usr/local/lib/libopencv_features2d.a /usr/local/lib/libopencv_highgui.a /usr/local/share/OpenCV/3rdparty/lib/libIlmImf.a /usr/local/share/OpenCV/3rdparty/lib/liblibjasper.a /usr/local/lib/libopencv_imgproc.a /usr/local/lib/libopencv_flann.a /usr/local/lib/libopencv_core.a /usr/lib64/libtiff.a /usr/lib64/libpng12.a /opt/libjpeg-turbo/lib64/libjpeg.a /usr/lib64/libz.a -lswscale -lavutil -lavformat -lavcodec -lrt -lpthread -lm -ldl -lstdc++ # statically linked version

CFLAGS= -fPIC -Wall -g -std=c++11 $(LAPACKCFLAGS) -fopenmp -DUSE_OPENMP -O3 $(shell pkg-config --static --cflags --libs opencv)
LDFLAGS=-fPIC -Wall -g -ljpeg -lpng -fopenmp $(shell pkg-config --static --cflags --libs opencv) 

SOURCES := $(shell find . -name '*.cpp')
OBJ := $(SOURCES:%.cpp=%.o)
HEADERS := $(shell find . -name '*.h')

all: deepmatching

.cpp.o:  %.cpp %.h
	$(CC) -o $@ $(CFLAGS) -c $+

deepmatching: $(HEADERS) $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS) $(LAPACKLDFLAGS)

deepmatching-static: $(HEADERS) $(OBJ)
	$(CC) -o $@ $^ $(STATICLAPACKLDFLAGS)


clean:
	rm -f $(OBJ) deepmatching *~ *.pyc .gdb_history
