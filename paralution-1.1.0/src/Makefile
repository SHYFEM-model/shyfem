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

FEMDIR = ../..
FEMBIN = $(FEMDIR)/fembin

include $(FEMDIR)/Rules.make

#CFLAGS =-O3 -fopenmp -fPIC
# For debug mode
#CFLAGS =-ggdb -Wall -fPIC
#IC_MKL =-mkl

# For systems with cross 32-64bit lib support
# (this is mostly for MacOS)
# please add -m64 to CFLAGS and CFLAGS

OMPFLAGS = -fopenmp

# CUDA GPU support 
ifeq ($(GPU),CUDA)
  SUPPORT_CUDA =-DSUPPORT_CUDA
  CUDA_OBJ = base/gpu/backend_gpu.o base/gpu/gpu_allocate_free.o base/gpu/gpu_matrix_bcsr.o base/gpu/gpu_matrix_coo.o base/gpu/gpu_matrix_csr.o base/gpu/gpu_matrix_dense.o base/gpu/gpu_matrix_dia.o base/gpu/gpu_matrix_ell.o base/gpu/gpu_matrix_hyb.o base/gpu/gpu_matrix_mcsr.o base/gpu/gpu_vector.o
endif

# OpenCL support
ifeq ($(GPU),OpenCL)
  SUPPORT_OCL =-DSUPPORT_OCL
  LOCL=-lOpenCL -L/usr/lib/x86_64-linux-gnu
  OCL_OBJ =  base/ocl/backend_ocl.o base/ocl/ocl_allocate_free.o base/ocl/ocl_matrix_bcsr.o base/ocl/ocl_matrix_coo.o  base/ocl/ocl_matrix_csr.o base/ocl/ocl_matrix_dense.o base/ocl/ocl_matrix_dia.o base/ocl/ocl_matrix_ell.o base/ocl/ocl_matrix_hyb.o base/ocl/ocl_matrix_mcsr.o base/ocl/ocl_vector.o
endif

# MIC/Intel Xeon Phi support
ifeq ($(GPU),MIC)
  SUPPORT_MIC =-DSUPPORT_MIC
  MIC_CC = icc
  #MIC_FLAGS = -O3 -Wall -openmp -vec-report3 -offload-option,mic,compiler,"-z defs"
  MIC_FLAGS = -O3 -Wall -openmp -offload-option,mic,compiler,"-z defs" -fPIC
  LMIC= -offload-option,mic,compiler,"-z defs"
  MIC_OBJ = base/mic/mic_allocate_free.o  base/mic/backend_mic.o base/mic/mic_vector.o  base/mic/mic_matrix_csr.o base/mic/mic_matrix_coo.o base/mic/mic_matrix_dense.o base/mic/mic_matrix_mcsr.o base/mic/mic_matrix_bcsr.o base/mic/mic_matrix_dia.o base/mic/mic_matrix_ell.o base/mic/mic_matrix_hyb.o
endif


#  ------------------------------------------------------
# Special Libraries
#  ------------------------------------------------------

# Uncomment for MKL support
# SUPPORT_MKL =-DSUPPORT_MKL $(IC_MKL) 



# ----------------------------------------------------------------------------------------

SUPPORT_LIB =$(SUPPORT_MKL) $(SUPPORT_CUDA) $(SUPPORT_OCL) $(SUPPORT_MIC)
LFLAGS=$(OPENMPFLAG) $(LIBCUDA) $(IC_MKL) $(LOCL) $(LMIC)

#PARALUTION_LIB =libparalution.a
PARALUTION_LIB =../../femlib/libparalution.a
PARALUTION_INC =-I.


HOST_OBJ = base/host/host_affinity.o base/host/host_conversion.o base/host/host_io.o  base/host/host_matrix_bcsr.o base/host/host_matrix_coo.o base/host/host_matrix_csr.o base/host/host_matrix_dense.o base/host/host_matrix_dia.o base/host/host_matrix_ell.o base/host/host_matrix_hyb.o base/host/host_matrix_mcsr.o base/host/host_stencil_laplace2d.o base/host/host_vector.o


OBJ = $(OCL_OBJ) $(CUDA_OBJ) $(HOST_OBJ) $(MIC_OBJ) base/backend_manager.o base/base_matrix.o base/base_paralution.o base/base_stencil.o base/base_vector.o base/local_matrix.o base/local_stencil.o base/local_vector.o base/operator.o base/vector.o utils/allocate_free.o utils/log.o utils/math_functions.o utils/time_functions.o solvers/chebyshev.o solvers/iter_ctrl.o  solvers/mixed_precision.o solvers/solver.o solvers/deflation/dpcg.o solvers/direct/inversion.o solvers/direct/lu.o solvers/direct/qr.o solvers/eigenvalue/ampe_sira.o solvers/eigenvalue/cg_hn.o solvers/krylov/bicgstab.o  solvers/krylov/cg.o solvers/krylov/cr.o solvers/krylov/fgmres.o solvers/krylov/gmres.o solvers/krylov/idr.o solvers/multigrid/amg.o  solvers/multigrid/base_amg.o solvers/multigrid/base_multigrid.o solvers/multigrid/multigrid.o  solvers/preconditioners/preconditioner_ai.o solvers/preconditioners/preconditioner_as.o solvers/preconditioners/preconditioner_blockprecond.o solvers/preconditioners/preconditioner.o solvers/preconditioners/preconditioner_multicolored.o solvers/preconditioners/preconditioner_multicolored_gs.o solvers/preconditioners/preconditioner_multicolored_ilu.o solvers/preconditioners/preconditioner_multielimination.o solvers/preconditioners/preconditioner_saddlepoint.o

default: lib

install: all
	rm -rf ../build/
	mkdir ../build
	mkdir ../build/bin
	mkdir ../build/inc
	mkdir ../build/lib
	cp -rf examples/amg examples/as-precond examples/assembly2 examples/assembly examples/async examples/benchmark examples/bicgstab examples/block-precond examples/cg-amg examples/cg examples/cmk examples/complex examples/direct examples/dpcg examples/fixed-point examples/gmres examples/idr examples/me-preconditioner examples/mixed-precision examples/power-method examples/simple-spmv examples/sira examples/sp-precond examples/stencil examples/tns ../build/bin/
	cp $(PARALUTION_LIB) ../build/lib/
	find . -name '*.hpp' | cpio -pavd ../build/inc/
	cp -rf plug-ins/OpenFOAM ../build/inc/plug-ins/

all: 	lib amg as-precond assembly2 assembly async benchmark bicgstab block-precond cg-amg cg cmk complex direct dpcg fixed-point gmres idr me-preconditioner mixed-precision power-method simple-spmv sira sp-precond stencil tns


fem: lib

lib: 	$(OBJ)
# For "ar"
# For MIC use "xair"; for everything else use "ar"
#	xiar -qoffload-build rcs $(PARALUTION_LIB) $(OBJ)
#	ar rcs $(PARALUTION_LIB) $(OBJ)
ifeq ($(GPU),MIC)
	xiar -qoffload-build rcs $(PARALUTION_LIB) $(OBJ)
else
	ar rcs $(PARALUTION_LIB) $(OBJ)
endif
#
# For ".so"
#	$(CCP) -shared -Wl,-soname,$(PARALUTION_LIB) -o  $(PARALUTION_LIB) $(OBJ) -lc


clean: 
	rm -rf *.a *.o examples/*.o examples/amg examples/as-precond examples/assembly2 examples/assembly examples/async examples/benchmark examples/bicgstab examples/block-precond examples/cg-amg examples/cg examples/cmk examples/complex examples/direct examples/dpcg examples/fixed-point examples/gmres examples/idr examples/me-preconditioner examples/mixed-precision examples/power-method examples/simple-spmv examples/sira examples/sp-precond examples/stencil examples/tns tmp/*.o base/*.o base/host/*.o base/gpu/*.o base/ocl/*.o base/mic/*.o solvers/direct/*.o solvers/krylov/*.o solvers/multigrid/*.o solvers/*.o solvers/eigenvalue/*.o solvers/preconditioners/*.o solvers/deflation/*.o utils/*.o *.mtx *.vec *.dat examples/*.mtx examples/*.vec examples/*.dat ../thirdparty/matrix-market/*.o ../build/* libparalution.a $(PARALUTION_LIB)

cleanall: clean

#
# general
#
base/base_paralution.o: base/base_paralution.cpp base/base_paralution.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/base_paralution.cpp -o base/base_paralution.o

base/operator.o: base/operator.cpp base/operator.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/operator.cpp -o base/operator.o

base/vector.o: base/vector.cpp base/vector.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/vector.cpp -o base/vector.o


base/local_matrix.o: base/local_matrix.cpp base/local_matrix.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/local_matrix.cpp -o base/local_matrix.o

base/global_matrix.o: base/global_matrix.cpp base/global_matrix.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/global_matrix.cpp -o base/global_matrix.o


base/local_vector.o: base/local_vector.cpp base/local_vector.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/local_vector.cpp -o base/local_vector.o

base/global_vector.o: base/global_vector.cpp base/global_vector.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/global_vector.cpp -o base/global_vector.o

base/base_matrix.o: base/base_matrix.cpp base/base_matrix.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/base_matrix.cpp -o base/base_matrix.o


# 
# host 
#
base/host/host_matrix_csr.o: base/host/host_matrix_csr.cpp base/host/host_matrix_csr.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/host/host_matrix_csr.cpp -o base/host/host_matrix_csr.o

base/host/host_matrix_mcsr.o: base/host/host_matrix_mcsr.cpp base/host/host_matrix_mcsr.hpp
	$(CCP) $(CPFLAGS) $(SUPPORT_LIB) -c base/host/host_matrix_mcsr.cpp -o base/host/host_matrix_mcsr.o 

base/host/host_matrix_bcsr.o: base/host/host_matrix_bcsr.cpp base/host/host_matrix_bcsr.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/host/host_matrix_bcsr.cpp -o base/host/host_matrix_bcsr.o 

base/host/host_matrix_coo.o: base/host/host_matrix_coo.cpp base/host/host_matrix_coo.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/host/host_matrix_coo.cpp -o base/host/host_matrix_coo.o 

base/host/host_matrix_dia.o: base/host/host_matrix_dia.cpp base/host/host_matrix_dia.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/host/host_matrix_dia.cpp -o base/host/host_matrix_dia.o 

base/host/host_matrix_ell.o: base/host/host_matrix_ell.cpp base/host/host_matrix_ell.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/host/host_matrix_ell.cpp -o base/host/host_matrix_ell.o 

base/host/host_matrix_hyb.o: base/host/host_matrix_hyb.cpp base/host/host_matrix_hyb.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/host/host_matrix_hyb.cpp -o base/host/host_matrix_hyb.o 

base/host/host_matrix_dense.o: base/host/host_matrix_dense.cpp base/host/host_matrix_dense.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/host/host_matrix_dense.cpp -o base/host/host_matrix_dense.o 

base/host/host_vector.o: base/host/host_vector.cpp base/host/host_vector.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/host/host_vector.cpp -o base/host/host_vector.o

base/host/host_conversion.o: base/host/host_conversion.cpp base/host/host_conversion.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/host/host_conversion.cpp -o base/host/host_conversion.o

base/host/host_affinity.o: base/host/host_affinity.cpp base/host/host_affinity.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/host/host_affinity.cpp -o base/host/host_affinity.o

base/host/host_io.o: base/host/host_io.cpp base/host/host_io.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/host/host_io.cpp -o base/host/host_io.o

base/host/host_stencil_laplace2d.o: base/host/host_stencil_laplace2d.cpp base/host/host_stencil_laplace2d.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/host/host_stencil_laplace2d.cpp -o base/host/host_stencil_laplace2d.o

base/base_vector.o: base/base_vector.cpp base/base_vector.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/base_vector.cpp -o base/base_vector.o

base/backend_manager.o: base/backend_manager.cpp base/backend_manager.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/backend_manager.cpp -o base/backend_manager.o

base/base_stencil.o: base/base_stencil.cpp base/base_stencil.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/base_stencil.cpp -o base/base_stencil.o

base/local_stencil.o: base/local_stencil.cpp base/local_stencil.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/local_stencil.cpp -o base/local_stencil.o

#
# third party
#
../thirdparty/matrix-market/mmio.o: ../thirdparty/matrix-market/mmio.c ../thirdparty/matrix-market/mmio.h	
	$(CC) $(CPFLAGS) $(OMPFLAGS) -c ../thirdparty/matrix-market/mmio.c -o ../thirdparty/matrix-market/mmio.o

#
# solver
#
solvers/krylov/cg.o: solvers/krylov/cg.cpp solvers/krylov/cg.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/krylov/cg.cpp -o solvers/krylov/cg.o

solvers/krylov/cr.o: solvers/krylov/cr.cpp solvers/krylov/cr.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/krylov/cr.cpp -o solvers/krylov/cr.o

solvers/krylov/gmres.o: solvers/krylov/gmres.cpp solvers/krylov/gmres.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/krylov/gmres.cpp -o solvers/krylov/gmres.o

solvers/krylov/fgmres.o: solvers/krylov/fgmres.cpp solvers/krylov/fgmres.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/krylov/fgmres.cpp -o solvers/krylov/fgmres.o

solvers/krylov/bicgstab.o: solvers/krylov/bicgstab.cpp solvers/krylov/bicgstab.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/krylov/bicgstab.cpp -o solvers/krylov/bicgstab.o

solvers/krylov/idr.o: solvers/krylov/idr.cpp solvers/krylov/idr.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/krylov/idr.cpp -o solvers/krylov/idr.o

solvers/deflation/dpcg.o: solvers/deflation/dpcg.cpp solvers/deflation/dpcg.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/deflation/dpcg.cpp -o solvers/deflation/dpcg.o

solvers/eigenvalue/ampe_sira.o: solvers/eigenvalue/ampe_sira.cpp solvers/eigenvalue/ampe_sira.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/eigenvalue/ampe_sira.cpp -o solvers/eigenvalue/ampe_sira.o

solvers/eigenvalue/cg_hn.o: solvers/eigenvalue/cg_hn.cpp solvers/eigenvalue/cg_hn.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/eigenvalue/cg_hn.cpp -o solvers/eigenvalue/cg_hn.o

solvers/solver.o: solvers/solver.cpp solvers/solver.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/solver.cpp -o solvers/solver.o

solvers/multigrid/multigrid.o: solvers/multigrid/multigrid.cpp solvers/multigrid/multigrid.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/multigrid/multigrid.cpp -o solvers/multigrid/multigrid.o

solvers/multigrid/base_multigrid.o: solvers/multigrid/base_multigrid.cpp solvers/multigrid/base_multigrid.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/multigrid/base_multigrid.cpp -o solvers/multigrid/base_multigrid.o

solvers/multigrid/amg.o: solvers/multigrid/amg.cpp solvers/multigrid/amg.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/multigrid/amg.cpp -o solvers/multigrid/amg.o

solvers/multigrid/base_amg.o: solvers/multigrid/base_amg.cpp solvers/multigrid/base_amg.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/multigrid/base_amg.cpp -o solvers/multigrid/base_amg.o

solvers/chebyshev.o: solvers/chebyshev.cpp solvers/chebyshev.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/chebyshev.cpp -o solvers/chebyshev.o

solvers/preconditioners/preconditioner.o: solvers/preconditioners/preconditioner.cpp solvers/preconditioners/preconditioner.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/preconditioners/preconditioner.cpp -o solvers/preconditioners/preconditioner.o

solvers/preconditioners/preconditioner_multielimination.o: solvers/preconditioners/preconditioner_multielimination.cpp solvers/preconditioners/preconditioner_multielimination.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/preconditioners/preconditioner_multielimination.cpp -o solvers/preconditioners/preconditioner_multielimination.o

solvers/preconditioners/preconditioner_multicolored.o: solvers/preconditioners/preconditioner_multicolored.cpp solvers/preconditioners/preconditioner_multicolored.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/preconditioners/preconditioner_multicolored.cpp -o solvers/preconditioners/preconditioner_multicolored.o

solvers/preconditioners/preconditioner_multicolored_gs.o: solvers/preconditioners/preconditioner_multicolored_gs.cpp solvers/preconditioners/preconditioner_multicolored_gs.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/preconditioners/preconditioner_multicolored_gs.cpp -o solvers/preconditioners/preconditioner_multicolored_gs.o

solvers/preconditioners/preconditioner_multicolored_ilu.o: solvers/preconditioners/preconditioner_multicolored_ilu.cpp solvers/preconditioners/preconditioner_multicolored_ilu.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/preconditioners/preconditioner_multicolored_ilu.cpp -o solvers/preconditioners/preconditioner_multicolored_ilu.o

solvers/preconditioners/preconditioner_saddlepoint.o: solvers/preconditioners/preconditioner_saddlepoint.cpp solvers/preconditioners/preconditioner_saddlepoint.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/preconditioners/preconditioner_saddlepoint.cpp -o solvers/preconditioners/preconditioner_saddlepoint.o

solvers/preconditioners/preconditioner_blockprecond.o: solvers/preconditioners/preconditioner_blockprecond.cpp solvers/preconditioners/preconditioner_blockprecond.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/preconditioners/preconditioner_blockprecond.cpp -o solvers/preconditioners/preconditioner_blockprecond.o

solvers/preconditioners/preconditioner_as.o: solvers/preconditioners/preconditioner_as.cpp solvers/preconditioners/preconditioner_as.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/preconditioners/preconditioner_as.cpp -o solvers/preconditioners/preconditioner_as.o

solvers/preconditioners/preconditioner_ai.o: solvers/preconditioners/preconditioner_ai.cpp solvers/preconditioners/preconditioner_ai.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/preconditioners/preconditioner_ai.cpp -o solvers/preconditioners/preconditioner_ai.o

solvers/iter_ctrl.o: solvers/iter_ctrl.cpp solvers/iter_ctrl.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/iter_ctrl.cpp -o solvers/iter_ctrl.o

solvers/mixed_precision.o: solvers/mixed_precision.cpp solvers/mixed_precision.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/mixed_precision.cpp -o solvers/mixed_precision.o

solvers/direct/inversion.o: solvers/direct/inversion.cpp solvers/direct/inversion.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/direct/inversion.cpp -o solvers/direct/inversion.o

solvers/direct/lu.o: solvers/direct/lu.cpp solvers/direct/lu.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/direct/lu.cpp -o solvers/direct/lu.o

solvers/direct/qr.o: solvers/direct/qr.cpp solvers/direct/qr.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c solvers/direct/qr.cpp -o solvers/direct/qr.o

#
# utils
#
utils/allocate_free.o: utils/allocate_free.cpp utils/allocate_free.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c utils/allocate_free.cpp -o utils/allocate_free.o

utils/math_functions.o: utils/math_functions.cpp utils/math_functions.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c utils/math_functions.cpp -o utils/math_functions.o

utils/time_functions.o: utils/time_functions.cpp utils/time_functions.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c utils/time_functions.cpp -o utils/time_functions.o

utils/log.o: utils/log.cpp utils/log.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c utils/log.cpp -o utils/log.o

#
# cuda
#
base/gpu/gpu_allocate_free.o: base/gpu/gpu_allocate_free.cu base/gpu/gpu_allocate_free.hpp
	$(NCC) $(NCFLAGS) $(SUPPORT_LIB) -c base/gpu/gpu_allocate_free.cu -o base/gpu/gpu_allocate_free.o

base/gpu/backend_gpu.o: base/gpu/backend_gpu.cu base/gpu/backend_gpu.hpp
	$(NCC) $(NCFLAGS) $(SUPPORT_LIB) -c base/gpu/backend_gpu.cu -o base/gpu/backend_gpu.o

base/gpu/gpu_vector.o: base/gpu/gpu_vector.cu base/gpu/gpu_vector.hpp
	$(NCC) $(NCFLAGS) $(SUPPORT_LIB) -c base/gpu/gpu_vector.cu -o base/gpu/gpu_vector.o

base/gpu/gpu_matrix_csr.o: base/gpu/gpu_matrix_csr.cu base/gpu/gpu_matrix_csr.hpp
	$(NCC) $(NCFLAGS) $(SUPPORT_LIB) -c base/gpu/gpu_matrix_csr.cu -o base/gpu/gpu_matrix_csr.o 

base/gpu/gpu_matrix_coo.o: base/gpu/gpu_matrix_coo.cu base/gpu/gpu_matrix_coo.hpp
	$(NCC) $(NCFLAGS) $(SUPPORT_LIB) -c base/gpu/gpu_matrix_coo.cu -o base/gpu/gpu_matrix_coo.o

base/gpu/gpu_matrix_dense.o: base/gpu/gpu_matrix_dense.cu base/gpu/gpu_matrix_dense.hpp
	$(NCC) $(NCFLAGS) $(SUPPORT_LIB) -c base/gpu/gpu_matrix_dense.cu -o base/gpu/gpu_matrix_dense.o

base/gpu/gpu_matrix_mcsr.o: base/gpu/gpu_matrix_mcsr.cu base/gpu/gpu_matrix_mcsr.hpp
	$(NCC) $(NCFLAGS) $(SUPPORT_LIB) -c base/gpu/gpu_matrix_mcsr.cu -o base/gpu/gpu_matrix_mcsr.o

base/gpu/gpu_matrix_bcsr.o: base/gpu/gpu_matrix_bcsr.cu base/gpu/gpu_matrix_bcsr.hpp
	$(NCC) $(NCFLAGS) $(SUPPORT_LIB) -c base/gpu/gpu_matrix_bcsr.cu -o base/gpu/gpu_matrix_bcsr.o

base/gpu/gpu_matrix_ell.o: base/gpu/gpu_matrix_ell.cu base/gpu/gpu_matrix_ell.hpp
	$(NCC) $(NCFLAGS) $(SUPPORT_LIB) -c base/gpu/gpu_matrix_ell.cu -o base/gpu/gpu_matrix_ell.o

base/gpu/gpu_matrix_dia.o: base/gpu/gpu_matrix_dia.cu base/gpu/gpu_matrix_dia.hpp
	$(NCC) $(NCFLAGS) $(SUPPORT_LIB) -c base/gpu/gpu_matrix_dia.cu -o base/gpu/gpu_matrix_dia.o

base/gpu/gpu_matrix_hyb.o: base/gpu/gpu_matrix_hyb.cu base/gpu/gpu_matrix_hyb.hpp
	$(NCC) $(NCFLAGS) $(SUPPORT_LIB) -c base/gpu/gpu_matrix_hyb.cu -o base/gpu/gpu_matrix_hyb.o


#
# OpenCL
#
base/ocl/ocl_allocate_free.o: base/ocl/ocl_allocate_free.cpp base/ocl/ocl_allocate_free.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/ocl/ocl_allocate_free.cpp -o base/ocl/ocl_allocate_free.o

base/ocl/backend_ocl.o: base/ocl/backend_ocl.cpp base/ocl/backend_ocl.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/ocl/backend_ocl.cpp -o base/ocl/backend_ocl.o

base/ocl/ocl_vector.o: base/ocl/ocl_vector.cpp base/ocl/ocl_vector.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/ocl/ocl_vector.cpp -o base/ocl/ocl_vector.o

base/ocl/ocl_matrix_csr.o: base/ocl/ocl_matrix_csr.cpp base/ocl/ocl_matrix_csr.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/ocl/ocl_matrix_csr.cpp -o base/ocl/ocl_matrix_csr.o

base/ocl/ocl_matrix_dense.o: base/ocl/ocl_matrix_dense.cpp base/ocl/ocl_matrix_dense.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/ocl/ocl_matrix_dense.cpp -o base/ocl/ocl_matrix_dense.o

base/ocl/ocl_matrix_mcsr.o: base/ocl/ocl_matrix_mcsr.cpp base/ocl/ocl_matrix_mcsr.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/ocl/ocl_matrix_mcsr.cpp -o base/ocl/ocl_matrix_mcsr.o

base/ocl/ocl_matrix_ell.o: base/ocl/ocl_matrix_ell.cpp base/ocl/ocl_matrix_ell.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/ocl/ocl_matrix_ell.cpp -o base/ocl/ocl_matrix_ell.o

base/ocl/ocl_matrix_dia.o: base/ocl/ocl_matrix_dia.cpp base/ocl/ocl_matrix_dia.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/ocl/ocl_matrix_dia.cpp -o base/ocl/ocl_matrix_dia.o

base/ocl/ocl_matrix_bcsr.o: base/ocl/ocl_matrix_bcsr.cpp base/ocl/ocl_matrix_bcsr.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/ocl/ocl_matrix_bcsr.cpp -o base/ocl/ocl_matrix_bcsr.o

base/ocl/ocl_matrix_coo.o: base/ocl/ocl_matrix_coo.cpp base/ocl/ocl_matrix_coo.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/ocl/ocl_matrix_coo.cpp -o base/ocl/ocl_matrix_coo.o

base/ocl/ocl_matrix_hyb.o: base/ocl/ocl_matrix_hyb.cpp base/ocl/ocl_matrix_hyb.hpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) -c base/ocl/ocl_matrix_hyb.cpp -o base/ocl/ocl_matrix_hyb.o

utils/ocl_check_hw: utils/ocl_check_hw.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) utils/ocl_check_hw.cpp -o utils/ocl_check_hw $(LFLAGS)
	utils/ocl_check_hw

#
# mic
#
base/mic/mic_allocate_free.o: base/mic/mic_allocate_free.cpp base/mic/mic_allocate_free.hpp
	$(MIC_CC) $(MIC_FLAGS) $(SUPPORT_LIB) -c base/mic/mic_allocate_free.cpp -o base/mic/mic_allocate_free.o

base/mic/backend_mic.o: base/mic/backend_mic.cpp base/mic/backend_mic.hpp
	$(MIC_CC) $(MIC_FLAGS) $(SUPPORT_LIB) -c base/mic/backend_mic.cpp -o base/mic/backend_mic.o

base/mic/mic_vector.o: base/mic/mic_vector.cpp base/mic/mic_vector.hpp
	$(MIC_CC) $(MIC_FLAGS) $(SUPPORT_LIB) -c base/mic/mic_vector.cpp -o base/mic/mic_vector.o

base/mic/mic_matrix_csr.o: base/mic/mic_matrix_csr.cpp base/mic/mic_matrix_csr.hpp
	$(MIC_CC) $(MIC_FLAGS) $(SUPPORT_LIB) -c base/mic/mic_matrix_csr.cpp -o base/mic/mic_matrix_csr.o 

base/mic/mic_matrix_coo.o: base/mic/mic_matrix_coo.cpp base/mic/mic_matrix_coo.hpp
	$(MIC_CC) $(MIC_FLAGS) $(SUPPORT_LIB) -c base/mic/mic_matrix_coo.cpp -o base/mic/mic_matrix_coo.o

base/mic/mic_matrix_dense.o: base/mic/mic_matrix_dense.cpp base/mic/mic_matrix_dense.hpp
	$(MIC_CC) $(MIC_FLAGS) $(SUPPORT_LIB) -c base/mic/mic_matrix_dense.cpp -o base/mic/mic_matrix_dense.o

base/mic/mic_matrix_mcsr.o: base/mic/mic_matrix_mcsr.cpp base/mic/mic_matrix_mcsr.hpp
	$(MIC_CC) $(MIC_FLAGS) $(SUPPORT_LIB) -c base/mic/mic_matrix_mcsr.cpp -o base/mic/mic_matrix_mcsr.o

base/mic/mic_matrix_bcsr.o: base/mic/mic_matrix_bcsr.cpp base/mic/mic_matrix_bcsr.hpp
	$(MIC_CC) $(MIC_FLAGS) $(SUPPORT_LIB) -c base/mic/mic_matrix_bcsr.cpp -o base/mic/mic_matrix_bcsr.o

base/mic/mic_matrix_ell.o: base/mic/mic_matrix_ell.cpp base/mic/mic_matrix_ell.hpp
	$(MIC_CC) $(MIC_FLAGS) $(SUPPORT_LIB) -c base/mic/mic_matrix_ell.cpp -o base/mic/mic_matrix_ell.o

base/mic/mic_matrix_dia.o: base/mic/mic_matrix_dia.cpp base/mic/mic_matrix_dia.hpp
	$(MIC_CC) $(MIC_FLAGS) $(SUPPORT_LIB) -c base/mic/mic_matrix_dia.cpp -o base/mic/mic_matrix_dia.o

base/mic/mic_matrix_hyb.o: base/mic/mic_matrix_hyb.cpp base/mic/mic_matrix_hyb.hpp
	$(MIC_CC) $(MIC_FLAGS) $(SUPPORT_LIB) -c base/mic/mic_matrix_hyb.cpp -o base/mic/mic_matrix_hyb.o

#
# examples
#
examples/amg.o: examples/amg.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/amg.cpp -o examples/amg.o

amg: examples/amg.o
	$(CCP) -o examples/amg examples/amg.o $(PARALUTION_LIB) $(LFLAGS)

examples/cg-amg.o: examples/cg-amg.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/cg-amg.cpp -o examples/cg-amg.o

cg-amg: examples/cg-amg.o
	$(CCP) -o examples/cg-amg examples/cg-amg.o $(PARALUTION_LIB) $(LFLAGS)

examples/cg.o: examples/cg.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/cg.cpp -o examples/cg.o

cg: examples/cg.o
	$(CCP) -o examples/cg examples/cg.o $(PARALUTION_LIB) $(LFLAGS)

examples/gmres.o: examples/gmres.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/gmres.cpp -o examples/gmres.o

gmres: examples/gmres.o
	$(CCP) -o examples/gmres examples/gmres.o $(PARALUTION_LIB) $(LFLAGS)

examples/sp-precond.o: examples/sp-precond.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/sp-precond.cpp -o examples/sp-precond.o

sp-precond: examples/sp-precond.o
	$(CCP) -o examples/sp-precond examples/sp-precond.o $(PARALUTION_LIB) $(LFLAGS)

examples/assembly.o: examples/assembly.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/assembly.cpp -o examples/assembly.o

assembly: examples/assembly.o
	$(CCP) -o examples/assembly examples/assembly.o $(PARALUTION_LIB) $(LFLAGS)

examples/assembly2.o: examples/assembly2.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/assembly2.cpp -o examples/assembly2.o

assembly2: examples/assembly2.o
	$(CCP) -o examples/assembly2 examples/assembly2.o $(PARALUTION_LIB) $(LFLAGS)

examples/as-precond.o: examples/as-precond.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/as-precond.cpp -o examples/as-precond.o

as-precond: examples/as-precond.o
	$(CCP) -o examples/as-precond examples/as-precond.o $(PARALUTION_LIB) $(LFLAGS)

examples/direct.o: examples/direct.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/direct.cpp -o examples/direct.o

direct: examples/direct.o
	$(CCP) -o examples/direct examples/direct.o $(PARALUTION_LIB) $(LFLAGS)

examples/block-precond.o: examples/block-precond.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/block-precond.cpp -o examples/block-precond.o

block-precond: examples/block-precond.o
	$(CCP) -o examples/block-precond examples/block-precond.o $(PARALUTION_LIB) $(LFLAGS)

examples/idr.o: examples/idr.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/idr.cpp -o examples/idr.o

idr: examples/idr.o
	$(CCP) -o examples/idr examples/idr.o $(PARALUTION_LIB) $(LFLAGS)

examples/dpcg.o: examples/dpcg.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/dpcg.cpp -o examples/dpcg.o

dpcg: examples/dpcg.o
	$(CCP) -o examples/dpcg examples/dpcg.o $(PARALUTION_LIB) $(LFLAGS)

examples/bicgstab.o: examples/bicgstab.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/bicgstab.cpp -o examples/bicgstab.o

bicgstab: examples/bicgstab.o
	$(CCP) -o examples/bicgstab examples/bicgstab.o $(PARALUTION_LIB) $(LFLAGS) 

examples/fixed-point.o: examples/fixed-point.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/fixed-point.cpp -o examples/fixed-point.o

fixed-point: examples/fixed-point.o
	$(CCP) -o examples/fixed-point examples/fixed-point.o $(PARALUTION_LIB) $(LFLAGS) 

examples/benchmark.o: examples/benchmark.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/benchmark.cpp -o examples/benchmark.o

benchmark: examples/benchmark.o
	$(CCP) -o examples/benchmark examples/benchmark.o $(PARALUTION_LIB) $(LFLAGS) 

examples/power-method.o: examples/power-method.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/power-method.cpp -o examples/power-method.o

power-method: examples/power-method.o
	$(CCP) -o examples/power-method examples/power-method.o $(PARALUTION_LIB) $(LFLAGS) 

examples/simple-spmv.o: examples/simple-spmv.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/simple-spmv.cpp -o examples/simple-spmv.o

simple-spmv: examples/simple-spmv.o
	$(CCP) -o examples/simple-spmv examples/simple-spmv.o $(PARALUTION_LIB) $(LFLAGS)

examples/me-preconditioner.o: examples/me-preconditioner.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/me-preconditioner.cpp -o examples/me-preconditioner.o

me-preconditioner: examples/me-preconditioner.o
	$(CCP) -o examples/me-preconditioner examples/me-preconditioner.o $(PARALUTION_LIB) $(LFLAGS)

examples/mixed-precision.o: examples/mixed-precision.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/mixed-precision.cpp -o examples/mixed-precision.o

mixed-precision: examples/mixed-precision.o
	$(CCP) -o examples/mixed-precision examples/mixed-precision.o $(PARALUTION_LIB) $(LFLAGS)

examples/async.o: examples/async.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/async.cpp -o examples/async.o

async: examples/async.o
	$(CCP) -o examples/async examples/async.o $(PARALUTION_LIB) $(LFLAGS)


examples/cmk.o: examples/cmk.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/cmk.cpp -o examples/cmk.o

cmk: examples/cmk.o
	$(CCP) -o examples/cmk examples/cmk.o $(PARALUTION_LIB) $(LFLAGS)

examples/complex.o: examples/complex.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(SUPPORT_LIB) $(PARALUTION_INC) -c examples/complex.cpp -o examples/complex.o

complex: examples/complex.o
	$(CCP) -o examples/complex examples/complex.o $(PARALUTION_LIB) $(LFLAGS)

examples/sira.o: examples/sira.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/sira.cpp -o examples/sira.o

sira: examples/sira.o
	$(CCP) -o examples/sira examples/sira.o $(PARALUTION_LIB) $(LFLAGS)

examples/stencil.o: examples/stencil.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/stencil.cpp -o examples/stencil.o

stencil: examples/stencil.o
	$(CCP) -o examples/stencil examples/stencil.o $(PARALUTION_LIB) $(LFLAGS)

examples/tns.o: examples/tns.cpp
	$(CCP) $(CPFLAGS) $(OMPFLAGS) $(PARALUTION_INC) -c examples/tns.cpp -o examples/tns.o

tns: examples/tns.o
	$(CCP) -o examples/tns examples/tns.o $(PARALUTION_LIB) $(LFLAGS)
