Full compilation on my 64bit machine:
make MATLABDIR=/opt/matlab LDIRS="-L /opt/matlab/bin/glnxa64/" CFLAGS="-O3 -march=nocona -ffast-math -fomit-frame-pointer -fPIC"

libchi2.so compilation (the only thing needed for python):
make libchi2.so
python chi2.py





