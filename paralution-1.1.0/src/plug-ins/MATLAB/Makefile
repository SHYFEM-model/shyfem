# Compiler
OCTAVE=mkoctfile
MATLAB=mex

# Flags
FLAGS= -lparalution

# Include
INC=-I../../../build/inc
LIB=-L../../../build/lib

default: octave

octave: octave_paralution_pcg

matlab: matlab_paralution_pcg

octave_paralution_pcg: octave_paralution_pcg.o
	$(OCTAVE) --mex octave_paralution_pcg.o $(LIB) $(FLAGS)
	mv octave_paralution_pcg.mex paralution_pcg.mex

octave_paralution_pcg.o: paralution_pcg.cpp
	$(OCTAVE) --mex -c paralution_pcg.cpp -o octave_paralution_pcg.o $(INC)

matlab_paralution_pcg: paralution_pcg.cpp
	$(MATLAB)  -largeArrayDims paralution_pcg.cpp $(INC) $(LIB) $(FLAGS)


clean:
	rm -rf *.o *.mex *.mexa64
