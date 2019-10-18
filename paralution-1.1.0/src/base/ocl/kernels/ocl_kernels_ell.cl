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


// Nathan Bell and Michael Garland
// Efficient Sparse Matrix-Vector Multiplication on {CUDA}
// NVR-2008-004 / NVIDIA Technical Report
__kernel void kernel_ell_spmv(         const IndexType num_rows,
                                       const IndexType num_cols,
                                       const IndexType num_cols_per_row,
                              __global const IndexType *Acol,
                              __global const ValueType *Aval,
                              __global const ValueType *x,
                              __global       ValueType *y) {

  IndexType row = get_global_id(0);

  if (row < num_rows) {

    ValueType sum = (ValueType) 0;

    for (IndexType n=0; n<num_cols_per_row; ++n) {

      const IndexType ind = n * num_rows + row;
      const IndexType col = Acol[ind];

      if ((col >= 0) && (col < num_cols))
        sum += Aval[ind] * x[col];

    }

    y[row] = sum;

  }

}

// Nathan Bell and Michael Garland
// Efficient Sparse Matrix-Vector Multiplication on {CUDA}
// NVR-2008-004 / NVIDIA Technical Report
__kernel void kernel_ell_add_spmv(         const IndexType num_rows,
                                           const IndexType num_cols,
                                           const IndexType num_cols_per_row,
                                  __global const IndexType *Acol,
                                  __global const ValueType *Aval,
                                           const ValueType scalar,
                                  __global const ValueType *x,
                                  __global       ValueType *y) {

  IndexType row = get_global_id(0);

  if (row < num_rows) {

    ValueType sum = (ValueType) 0;

    for (IndexType n=0; n<num_cols_per_row; ++n) {

      const IndexType ind = n * num_rows + row;
      const IndexType col = Acol[ind];
      
      if ((col >= 0) && (col < num_cols))
        sum += Aval[ind] * x[col];

    }
        
    y[row] += scalar * sum;

  }

}

__kernel void kernel_ell_max_row(         const IndexType nrow,
                                 __global const IndexType *data,
                                 __global       IndexType *out,
                                          const IndexType  GROUP_SIZE,
                                          const IndexType  LOCAL_SIZE) {

    IndexType tid = get_local_id(0);

    __local IndexType sdata[BLOCK_SIZE];

    sdata[tid] = 0;

    IndexType max;

    IndexType gid = GROUP_SIZE * get_group_id(0) + tid;

    for (IndexType i = 0; i < LOCAL_SIZE; ++i, gid += BLOCK_SIZE) {

      if (gid < nrow) {
        max = data[gid+1] - data[gid];
        if (max > sdata[tid])
          sdata[tid] = max;
      }

    }

    barrier(CLK_LOCAL_MEM_FENCE);

    for (IndexType i = BLOCK_SIZE/2; i > 0; i /= 2) {

      if (tid < i)
        if (sdata[tid+i] > sdata[tid]) sdata[tid] = sdata[tid+i];

      barrier(CLK_LOCAL_MEM_FENCE);

    }

    if (tid == 0)
      out[get_group_id(0)] = sdata[tid];

}

__kernel void kernel_ell_csr_to_ell(         const IndexType nrow,
                                             const IndexType max_row,
                                    __global const IndexType *src_row_offset,
                                    __global const IndexType *src_col,
                                    __global const ValueType *src_val,
                                    __global       IndexType *ell_col,
                                    __global       ValueType *ell_val) {

  IndexType ai = get_global_id(0);
  IndexType aj;
  IndexType n = 0;
  IndexType ell_ind;

  if (ai < nrow) {

    for (aj=src_row_offset[ai]; aj<src_row_offset[ai+1]; ++aj) {

      ell_ind = n * nrow + ai;

      ell_col[ell_ind] = src_col[aj];
      ell_val[ell_ind] = src_val[aj];

      ++n;

    }

    for (aj=src_row_offset[ai+1]-src_row_offset[ai]; aj<max_row; ++aj) {

      ell_ind = n * nrow + ai;

      ell_col[ell_ind] = (int) -1;
      ell_val[ell_ind] = (ValueType) 0;

      ++n;

    }

  }

}

