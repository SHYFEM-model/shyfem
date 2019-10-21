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


#ifndef PARALUTION_OCL_KERNELS_HYB_HPP_
#define PARALUTION_OCL_KERNELS_HYB_HPP_

namespace paralution {

const char *ocl_kernels_hyb =
	"__kernel void kernel_ell_nnz_coo(         const IndexType nrow,\n"
	"                                          const IndexType max_row,\n"
	"                                 __global const IndexType *row_offset,\n"
	"                                 __global       IndexType *nnz_coo) {\n"
	"\n"
	"  IndexType gid = get_global_id(0);\n"
	"\n"
	"  if (gid < nrow) {\n"
	"\n"
	"    nnz_coo[gid] = 0;\n"
	"    IndexType nnz_per_row = row_offset[gid+1] - row_offset[gid];\n"
	"\n"
	"    if (nnz_per_row > max_row)\n"
	"      nnz_coo[gid] = nnz_per_row - max_row;\n"
	"\n"
	"  }\n"
	"\n"
	"}\n"
	"\n"
	"__kernel void kernel_ell_fill_ell(         const IndexType nrow,\n"
	"                                           const IndexType max_row,\n"
	"                                  __global const IndexType *row_offset,\n"
	"                                  __global const IndexType *col,\n"
	"                                  __global const ValueType *val,\n"
	"                                  __global       IndexType *ELL_col,\n"
	"                                  __global       ValueType *ELL_val,\n"
	"                                  __global       IndexType *nnz_ell) {\n"
	"\n"
	"  IndexType gid = get_global_id(0);\n"
	"\n"
	"  if (gid < nrow) {\n"
	"\n"
	"    IndexType n = 0;\n"
	"\n"
	"    for (IndexType i=row_offset[gid]; i<row_offset[gid+1]; ++i) {\n"
	"\n"
	"      if (n >= max_row) break;\n"
	"\n"
	"      IndexType idx = n * nrow + gid;\n"
	"\n"
	"      ELL_col[idx] = col[i];\n"
	"      ELL_val[idx] = val[i];\n"
	"\n"
	"      ++n;\n"
	"\n"
	"    }\n"
	"\n"
	"    nnz_ell[gid] = n;\n"
	"\n"
	"  }\n"
	"\n"
	"}\n"
	"\n"
	"__kernel void kernel_ell_fill_coo(         const IndexType nrow,\n"
	"                                  __global const IndexType *row_offset,\n"
	"                                  __global const IndexType *col,\n"
	"                                  __global const ValueType *val,\n"
	"                                  __global const IndexType *nnz_coo,\n"
	"                                  __global const IndexType *nnz_ell,\n"
	"                                  __global       IndexType *COO_row,\n"
	"                                  __global       IndexType *COO_col,\n"
	"                                  __global       ValueType *COO_val) {\n"
	"\n"
	"  IndexType gid = get_global_id(0);\n"
	"\n"
	"  if (gid < nrow) {\n"
	"\n"
	"    IndexType row_ptr = row_offset[gid+1];\n"
	"\n"
	"    for (IndexType i=row_ptr - nnz_coo[gid]; i<row_ptr; ++i) {\n"
	"\n"
	"      IndexType idx = i - nnz_ell[gid];\n"
	"\n"
	"      COO_row[idx] = gid;\n"
	"      COO_col[idx] = col[i];\n"
	"      COO_val[idx] = val[i];\n"
	"\n"
	"    }\n"
	"\n"
	"  }\n"
	"\n"
	"}\n"
	"\n"
;
}

#endif // PARALUTION_OCL_KERNELS_HYB_HPP_
