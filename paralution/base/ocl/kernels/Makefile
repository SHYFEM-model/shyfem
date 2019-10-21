# Compiler
CC=g++

# Flags
FLAGS=-O3

# Include
INC=-I../../../../build/inc
LIB=../../../../build/lib/libparalution.so

default: cl_to_hpp hpp_to_cl

cl_to_hpp: cl_to_hpp.o
	$(CC) -o cl_to_hpp cl_to_hpp.o $(LIB) $(FLAGS) -Wall

cl_to_hpp.o: cl_to_hpp.cpp
	$(CC) $(INC) -c cl_to_hpp.cpp

hpp_to_cl: hpp_to_cl.o
	$(CC) -o hpp_to_cl hpp_to_cl.o $(LIB) $(FLAGS) -Wall

hpp_to_cl.o: hpp_to_cl.cpp
	$(CC) $(INC) -c hpp_to_cl.cpp

clean:
	rm -rf *.o cl_to_hpp hpp_to_cl
