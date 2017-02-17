clang++ -c mkNtuple.cc -I../include
clang++ mkNtuple.o \
	-L../lib -lmonoHZZ4L \
	-L$DELPHES -lDelphes -o mkNtuple
rm *.o


