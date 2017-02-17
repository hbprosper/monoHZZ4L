clang++ -c analyzeHZZ4L.cc -I../include
clang++ analyzeHZZ4L.o \
	-L../lib -lmonoHZZ4L \
	-L$DELPHES -lDelphes -o ../bin/analyzeHZZ4L
rm *.o


