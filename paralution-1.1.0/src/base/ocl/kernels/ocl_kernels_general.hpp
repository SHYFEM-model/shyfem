// **************************************************************************
//
//    PARALUTION   www.paralution.com
//
//    Copyright (C) 2015  PARALUTION Labs UG (haftungsbeschränkt) & Co. KG
//                        Am Hasensprung 6, 76571 Gaggenau
//                        Handelsregister: Amtsgericht Mannheim, HRA 706051
//                        Vertreten durch:
//                        PARALUTION Labs Verwaltungs UG (haftungsbeschränkt)
//                        Am Hasensprung 6, 76571 Gaggenau
//                        Handelsregister: Amtsgericht Mannheim, HRB 721277
//                        Geschäftsführer: Dimitar Lukarski, Nico Trost
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// **************************************************************************



// PARALUTION version 1.1.0 


#ifndef PARALUTION_OCL_KERNELS_GENERAL_HPP_
#define PARALUTION_OCL_KERNELS_GENERAL_HPP_

namespace paralution {

const char *ocl_kernels_general =
	"__kernel void kernel_red_recurse(__global IndexType *dst, __global const IndexType *src, const IndexType numElems) {\n"
	"\n"
	"  IndexType index = BLOCK_SIZE * get_global_id(0);\n"
	"\n"
	"  if (index >= numElems)\n"
	"    return;\n"
	"\n"
	"  IndexType i = index;\n"
	"\n"
	"  if (i < BLOCK_SIZE)\n"
	"    return;\n"
	"\n"
	"  IndexType a = 0;\n"
	"\n"
	"  while (i >= BLOCK_SIZE) {\n"
	"    a += src[i];\n"
	"    i -= BLOCK_SIZE;\n"
	"  }\n"
	"\n"
	"  dst[index] = a;\n"
	"\n"
	"}\n"
	"\n"
	"__kernel void kernel_red_partial_sum(__global       IndexType *dst,\n"
	"                                     __global const IndexType *src,\n"
	"                                              const IndexType numElems,\n"
	"                                              const IndexType shift) {\n"
	"\n"
	"  IndexType index = get_global_id(0);\n"
	"  IndexType tid   = get_local_id(0);\n"
	"  IndexType gid   = get_group_id(0);\n"
	"\n"
	"  if (index < numElems) {\n"
	"\n"
	"    __local IndexType data[BLOCK_SIZE];\n"
	"\n"
	"    data[tid] = src[index];\n"
	"\n"
	"    barrier(CLK_LOCAL_MEM_FENCE);\n"
	"\n"
	"    for (IndexType i = BLOCK_SIZE/2; i > 0; i/=2) {\n"
	"\n"
	"      if (tid < i)\n"
	"        data[tid] = data[tid] + data[tid+i];\n"
	"\n"
	"      barrier(CLK_LOCAL_MEM_FENCE);\n"
	"\n"
	"    }\n"
	"\n"
	"    if (tid == 0 && BLOCK_SIZE*(1+gid)-1<numElems)\n"
	"      dst[BLOCK_SIZE*(1+gid)-1+shift] = data[0];\n"
	"\n"
	"  }\n"
	"\n"
	"}\n"
	"\n"
	"__kernel void kernel_red_extrapolate(__global       IndexType *dst,\n"
	"                                     __global const IndexType *srcBorder,\n"
	"                                     __global const IndexType *srcData,\n"
	"                                                    IndexType  numElems,\n"
	"                                              const IndexType  shift) {\n"
	"\n"
	"  IndexType index = get_local_size(0) * get_local_id(0);\n"
	"\n"
	"  if (index < numElems-1) {\n"
	"\n"
	"    IndexType sum = srcBorder[index];\n"
	"\n"
	"    for(IndexType i = 0; i < get_local_size(0) && index+i<numElems; ++i) {\n"
	"      sum += srcData[index+i];\n"
	"      dst[index+i+shift] = sum;\n"
	"    }\n"
	"\n"
	"  }\n"
	"\n"
	"}\n"
	"\n"
;
}

#endif // PARALUTION_OCL_KERNELS_GENERAL_HPP_
