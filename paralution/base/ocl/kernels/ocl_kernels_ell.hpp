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


#ifndef PARALUTION_OCL_KERNELS_ELL_HPP_
#define PARALUTION_OCL_KERNELS_ELL_HPP_

namespace paralution {

const char *ocl_kernels_ell =
// Nathan Bell and Michael Garland
// Efficient Sparse Matrix-Vector Multiplication on {CUDA}
// NVR-2008-004 / NVIDIA Technical Report
	"__kernel void kernel_ell_spmv(         const IndexType num_rows,\n"
	"                                       const IndexType num_cols,\n"
	"                                       const IndexType num_cols_per_row,\n"
	"                              __global const IndexType *Acol,\n"
	"                              __global const ValueType *Aval,\n"
	"                              __global const ValueType *x,\n"
	"                              __global       ValueType *y) {\n"
	"\n"
	"  IndexType row = get_global_id(0);\n"
	"\n"
	"  if (row < num_rows) {\n"
	"\n"
	"    ValueType sum = (ValueType) 0;\n"
	"\n"
	"    for (IndexType n=0; n<num_cols_per_row; ++n) {\n"
	"\n"
	"      const IndexType ind = n * num_rows + row;\n"
	"      const IndexType col = Acol[ind];\n"
	"\n"
	"      if ((col >= 0) && (col < num_cols))\n"
	"        sum += Aval[ind] * x[col];\n"
	"\n"
	"    }\n"
	"\n"
	"    y[row] = sum;\n"
	"\n"
	"  }\n"
	"\n"
	"}\n"
	"\n"
// Nathan Bell and Michael Garland
// Efficient Sparse Matrix-Vector Multiplication on {CUDA}
// NVR-2008-004 / NVIDIA Technical Report
	"__kernel void kernel_ell_add_spmv(         const IndexType num_rows,\n"
	"                                           const IndexType num_cols,\n"
	"                                           const IndexType num_cols_per_row,\n"
	"                                  __global const IndexType *Acol,\n"
	"                                  __global const ValueType *Aval,\n"
	"                                           const ValueType scalar,\n"
	"                                  __global const ValueType *x,\n"
	"                                  __global       ValueType *y) {\n"
	"\n"
	"  IndexType row = get_global_id(0);\n"
	"\n"
	"  if (row < num_rows) {\n"
	"\n"
	"    ValueType sum = (ValueType) 0;\n"
	"\n"
	"    for (IndexType n=0; n<num_cols_per_row; ++n) {\n"
	"\n"
	"      const IndexType ind = n * num_rows + row;\n"
	"      const IndexType col = Acol[ind];\n"
	"      \n"
	"      if ((col >= 0) && (col < num_cols))\n"
	"        sum += Aval[ind] * x[col];\n"
	"\n"
	"    }\n"
	"        \n"
	"    y[row] += scalar * sum;\n"
	"\n"
	"  }\n"
	"\n"
	"}\n"
	"\n"
	"__kernel void kernel_ell_max_row(         const IndexType nrow,\n"
	"                                 __global const IndexType *data,\n"
	"                                 __global       IndexType *out,\n"
	"                                          const IndexType  GROUP_SIZE,\n"
	"                                          const IndexType  LOCAL_SIZE) {\n"
	"\n"
	"    IndexType tid = get_local_id(0);\n"
	"\n"
	"    __local IndexType sdata[BLOCK_SIZE];\n"
	"\n"
	"    sdata[tid] = 0;\n"
	"\n"
	"    IndexType max;\n"
	"\n"
	"    IndexType gid = GROUP_SIZE * get_group_id(0) + tid;\n"
	"\n"
	"    for (IndexType i = 0; i < LOCAL_SIZE; ++i, gid += BLOCK_SIZE) {\n"
	"\n"
	"      if (gid < nrow) {\n"
	"        max = data[gid+1] - data[gid];\n"
	"        if (max > sdata[tid])\n"
	"          sdata[tid] = max;\n"
	"      }\n"
	"\n"
	"    }\n"
	"\n"
	"    barrier(CLK_LOCAL_MEM_FENCE);\n"
	"\n"
	"    for (IndexType i = BLOCK_SIZE/2; i > 0; i /= 2) {\n"
	"\n"
	"      if (tid < i)\n"
	"        if (sdata[tid+i] > sdata[tid]) sdata[tid] = sdata[tid+i];\n"
	"\n"
	"      barrier(CLK_LOCAL_MEM_FENCE);\n"
	"\n"
	"    }\n"
	"\n"
	"    if (tid == 0)\n"
	"      out[get_group_id(0)] = sdata[tid];\n"
	"\n"
	"}\n"
	"\n"
	"__kernel void kernel_ell_csr_to_ell(         const IndexType nrow,\n"
	"                                             const IndexType max_row,\n"
	"                                    __global const IndexType *src_row_offset,\n"
	"                                    __global const IndexType *src_col,\n"
	"                                    __global const ValueType *src_val,\n"
	"                                    __global       IndexType *ell_col,\n"
	"                                    __global       ValueType *ell_val) {\n"
	"\n"
	"  IndexType ai = get_global_id(0);\n"
	"  IndexType aj;\n"
	"  IndexType n = 0;\n"
	"  IndexType ell_ind;\n"
	"\n"
	"  if (ai < nrow) {\n"
	"\n"
	"    for (aj=src_row_offset[ai]; aj<src_row_offset[ai+1]; ++aj) {\n"
	"\n"
	"      ell_ind = n * nrow + ai;\n"
	"\n"
	"      ell_col[ell_ind] = src_col[aj];\n"
	"      ell_val[ell_ind] = src_val[aj];\n"
	"\n"
	"      ++n;\n"
	"\n"
	"    }\n"
	"\n"
	"    for (aj=src_row_offset[ai+1]-src_row_offset[ai]; aj<max_row; ++aj) {\n"
	"\n"
	"      ell_ind = n * nrow + ai;\n"
	"\n"
	"      ell_col[ell_ind] = (int) -1;\n"
	"      ell_val[ell_ind] = (ValueType) 0;\n"
	"\n"
	"      ++n;\n"
	"\n"
	"    }\n"
	"\n"
	"  }\n"
	"\n"
	"}\n"
	"\n"
	"\n"
;
}

#endif // PARALUTION_OCL_KERNELS_ELL_HPP_
