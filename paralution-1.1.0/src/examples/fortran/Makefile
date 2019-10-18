# *************************************************************************
#
#    PARALUTION   www.paralution.com
#
#    Copyright (C) 2015  PARALUTION Labs UG (haftungsbeschränkt) & Co. KG
#                        Am Hasensprung 6, 76571 Gaggenau
#                        Handelsregister: Amtsgericht Mannheim, HRA 706051
#                        Vertreten durch:
#                        PARALUTION Labs Verwaltungs UG (haftungsbeschränkt)
#                        Am Hasensprung 6, 76571 Gaggenau
#                        Handelsregister: Amtsgericht Mannheim, HRB 721277
#                        Geschäftsführer: Dimitar Lukarski, Nico Trost
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# *************************************************************************



# PARALUTION version 1.1.0 


CC=g++
FC=gfortran
FLAGS=-O3 -lstdc++ -fopenmp -lcudart -lcusparse -lcublas -lOpenCL -L/usr/local/cuda/lib64
INC=-I../../../build/inc
LIB=../../../build/lib/libparalution.so
OBJ=mmio.o paralution_fortran.o paralution_solver.o

default: paralution_solver

paralution_solver: $(OBJ)
	$(FC) -o paralution_solver $(OBJ) $(LIB) $(FLAGS) -Wall

paralution_solver.o: paralution_solver.f90
	$(FC) $(FLAGS) -c paralution_solver.f90

mmio.o: mmio.f
	$(FC) $(FLAGS) -c mmio.f

paralution_fortran.o: ../../plug-ins/paralution_fortran.cpp
	$(CC) $(FLAGS) $(INC) -c ../../plug-ins/paralution_fortran.cpp

clean:
	rm -rf *.o paralution_solver
