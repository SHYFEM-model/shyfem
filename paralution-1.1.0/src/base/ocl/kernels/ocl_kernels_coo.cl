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


__kernel void kernel_coo_permute(         const IndexType  nnz,
                                 __global const IndexType *in_row,
                                 __global const IndexType *in_col,
                                 __global const IndexType *perm,
                                 __global       IndexType *out_row,
                                 __global       IndexType *out_col) {

  IndexType ind = get_global_id(0);

  for (IndexType i=ind; i<nnz; i+=get_local_size(0)) {

    out_row[i] = perm[ in_row[i] ];
    out_col[i] = perm[ in_col[i] ];

  }

}

// ----------------------------------------------------------
// function segreduce_warp(...)
// ----------------------------------------------------------
// Modified and adapted from CUSP 0.3.1, 
// http://code.google.com/p/cusp-library/
// NVIDIA, APACHE LICENSE 2.0
// ----------------------------------------------------------
// CHANGELOG
// - adapted interface
// ----------------------------------------------------------
inline ValueType segreduce_warp(const IndexType thread_lane, IndexType row, ValueType val,
                                __local IndexType *rows, __local ValueType *vals) {

  IndexType tid = get_local_id(0);

  rows[tid] = row;
  vals[tid] = val;
  
  if( thread_lane >=  1 && row == rows[tid -  1] ) { vals[tid] = val = val + vals[tid -  1]; } 
  if( thread_lane >=  2 && row == rows[tid -  2] ) { vals[tid] = val = val + vals[tid -  2]; }
  if( thread_lane >=  4 && row == rows[tid -  4] ) { vals[tid] = val = val + vals[tid -  4]; }
  if( thread_lane >=  8 && row == rows[tid -  8] ) { vals[tid] = val = val + vals[tid -  8]; }
  if( thread_lane >= 16 && row == rows[tid - 16] ) { vals[tid] = val = val + vals[tid - 16]; }
  
  return val;

}

// ----------------------------------------------------------
// function segreduce_block(...)
// ----------------------------------------------------------
// Modified and adapted from CUSP 0.3.1, 
// http://code.google.com/p/cusp-library/
// NVIDIA, APACHE LICENSE 2.0
// ----------------------------------------------------------
// CHANGELOG
// - adapted interface
// ----------------------------------------------------------
inline void segreduce_block(__local const IndexType *idx, __local ValueType *val) {

  ValueType left = (ValueType) 0;
  IndexType tid = get_local_id(0);

  if( tid >=   1 && idx[tid] == idx[tid -   1] ) { left = val[tid -   1]; } barrier(CLK_LOCAL_MEM_FENCE); val[tid] += left; left = 0; barrier(CLK_LOCAL_MEM_FENCE);  
  if( tid >=   2 && idx[tid] == idx[tid -   2] ) { left = val[tid -   2]; } barrier(CLK_LOCAL_MEM_FENCE); val[tid] += left; left = 0; barrier(CLK_LOCAL_MEM_FENCE);
  if( tid >=   4 && idx[tid] == idx[tid -   4] ) { left = val[tid -   4]; } barrier(CLK_LOCAL_MEM_FENCE); val[tid] += left; left = 0; barrier(CLK_LOCAL_MEM_FENCE);
  if( tid >=   8 && idx[tid] == idx[tid -   8] ) { left = val[tid -   8]; } barrier(CLK_LOCAL_MEM_FENCE); val[tid] += left; left = 0; barrier(CLK_LOCAL_MEM_FENCE);
  if( tid >=  16 && idx[tid] == idx[tid -  16] ) { left = val[tid -  16]; } barrier(CLK_LOCAL_MEM_FENCE); val[tid] += left; left = 0; barrier(CLK_LOCAL_MEM_FENCE);
  if( tid >=  32 && idx[tid] == idx[tid -  32] ) { left = val[tid -  32]; } barrier(CLK_LOCAL_MEM_FENCE); val[tid] += left; left = 0; barrier(CLK_LOCAL_MEM_FENCE);  
  if( tid >=  64 && idx[tid] == idx[tid -  64] ) { left = val[tid -  64]; } barrier(CLK_LOCAL_MEM_FENCE); val[tid] += left; left = 0; barrier(CLK_LOCAL_MEM_FENCE);
  if( tid >= 128 && idx[tid] == idx[tid - 128] ) { left = val[tid - 128]; } barrier(CLK_LOCAL_MEM_FENCE); val[tid] += left; left = 0; barrier(CLK_LOCAL_MEM_FENCE);
  if( tid >= 256 && idx[tid] == idx[tid - 256] ) { left = val[tid - 256]; } barrier(CLK_LOCAL_MEM_FENCE); val[tid] += left; left = 0; barrier(CLK_LOCAL_MEM_FENCE);

}

// ----------------------------------------------------------
// function kernel_spmv_coo_flat(...)
// ----------------------------------------------------------
// Modified and adapted from CUSP 0.3.1, 
// http://code.google.com/p/cusp-library/
// NVIDIA, APACHE LICENSE 2.0
// ----------------------------------------------------------
// CHANGELOG
// - adapted interface
// ----------------------------------------------------------
__kernel void kernel_coo_spmv_flat(         const IndexType  num_nonzeros,
                                            const IndexType  interval_size,
                                   __global const IndexType *I, 
                                   __global const IndexType *J, 
                                   __global const ValueType *V, 
                                            const ValueType  scalar,
                                   __global const ValueType *x, 
                                   __global       ValueType *y,
                                   __global       IndexType *temp_rows,
                                   __global       ValueType *temp_vals) {

  __local volatile IndexType rows[48 * (BLOCK_SIZE/32)];
  __local volatile ValueType vals[BLOCK_SIZE];

        IndexType tid         = get_local_id(0);
  const IndexType thread_id   = BLOCK_SIZE * get_group_id(0) + tid;
  const IndexType thread_lane = tid & (WARP_SIZE-1);
  const IndexType warp_id     = thread_id   / WARP_SIZE;

  const IndexType interval_begin = warp_id * interval_size;
  IndexType interval_end2 = interval_begin + interval_size;
  if (interval_end2 > num_nonzeros)
    interval_end2 = num_nonzeros;

  const IndexType interval_end = interval_end2;

  const IndexType idx = 16 * (tid/32 + 1) + tid;

  rows[idx - 16] = -1;

  if(interval_begin >= interval_end)
    return;

  if (thread_lane == 31) {
    rows[idx] = I[interval_begin]; 
    vals[tid] = (ValueType) 0;
  }

  for(IndexType n = interval_begin + thread_lane; n < interval_end; n += WARP_SIZE) {
    IndexType row = I[n];
    ValueType val = scalar * V[n] * x[J[n]];
    
    if (thread_lane == 0) {
      if(row == rows[idx + 31])
        val += vals[tid + 31];
      else
        y[rows[idx + 31]] += vals[tid + 31];
    }

    rows[idx] = row;
    vals[tid] = val;

    if (row == rows[idx -  1]) { vals[tid] = val = val + vals[tid -  1]; } 
    if (row == rows[idx -  2]) { vals[tid] = val = val + vals[tid -  2]; }
    if (row == rows[idx -  4]) { vals[tid] = val = val + vals[tid -  4]; }
    if (row == rows[idx -  8]) { vals[tid] = val = val + vals[tid -  8]; }
    if (row == rows[idx - 16]) { vals[tid] = val = val + vals[tid - 16]; }

    if(thread_lane < 31 && row != rows[idx + 1])
      y[row] += vals[tid];
  }

  if (thread_lane == 31) {
    temp_rows[warp_id] = rows[idx];
    temp_vals[warp_id] = vals[tid];
  }

}

// ----------------------------------------------------------
// function kernel_spmv_coo_reduce_update(...)
// ----------------------------------------------------------
// Modified and adapted from CUSP 0.3.1, 
// http://code.google.com/p/cusp-library/
// NVIDIA, APACHE LICENSE 2.0
// ----------------------------------------------------------
// CHANGELOG
// - adapted interface
// ----------------------------------------------------------
__kernel void kernel_coo_spmv_reduce_update(         const IndexType  num_warps,
                                            __global const IndexType *temp_rows,
                                            __global const ValueType *temp_vals,
                                            __global       ValueType *y) {

  __local IndexType rows[BLOCK_SIZE + 1];
  __local ValueType vals[BLOCK_SIZE + 1];

  IndexType tid = get_local_id(0);

  const IndexType end = num_warps - (num_warps & (BLOCK_SIZE - 1));

  if (tid == 0) {
    rows[BLOCK_SIZE] = (int) -1;
    vals[BLOCK_SIZE] = (ValueType) 0;
  }
  
  barrier(CLK_LOCAL_MEM_FENCE);

  IndexType i = tid;

  while (i < end) {
    rows[tid] = temp_rows[i];
    vals[tid] = temp_vals[i];

    barrier(CLK_LOCAL_MEM_FENCE);

    segreduce_block(rows, vals);

    if (rows[tid] != rows[tid + 1])
      y[rows[tid]] += vals[tid];

    barrier(CLK_LOCAL_MEM_FENCE);

    i += BLOCK_SIZE; 
  }

  if (end < num_warps) {
    if (i < num_warps) {
      rows[tid] = temp_rows[i];
      vals[tid] = temp_vals[i];
    } else {
      rows[tid] = (int) -1;
      vals[tid] = (ValueType) 0;
    }

    barrier(CLK_LOCAL_MEM_FENCE);
 
    segreduce_block(rows, vals);

    if (i < num_warps)
      if (rows[tid] != rows[tid + 1])
        y[rows[tid]] += vals[tid];
  }

}

// ----------------------------------------------------------
// function spmv_coo_serial_kernel(...)
// ----------------------------------------------------------
// Modified and adapted from CUSP 0.3.1, 
// http://code.google.com/p/cusp-library/
// NVIDIA, APACHE LICENSE 2.0
// ----------------------------------------------------------
// CHANGELOG
// - adapted interface
// ----------------------------------------------------------
__kernel void kernel_coo_spmv_serial(         const IndexType  num_entries,
                                     __global const IndexType *I, 
                                     __global const IndexType *J, 
                                     __global const ValueType *V, 
                                              const ValueType  scalar,
                                     __global const ValueType *x, 
                                     __global       ValueType *y,
                                              const IndexType  shift) {

  for(IndexType n = 0; n < num_entries-shift; n++)
    y[I[n+shift]] += scalar * V[n+shift] * x[J[n+shift]];

}

__kernel void kernel_coo_csr_to_coo(const IndexType nrow, __global const IndexType *row_offset, __global IndexType *row) {

  IndexType ai = get_global_id(0);
  IndexType aj;

  if (ai < nrow) {

    for (aj=row_offset[ai]; aj<row_offset[ai+1]; ++aj)
      row[aj] = ai;

  }

}
