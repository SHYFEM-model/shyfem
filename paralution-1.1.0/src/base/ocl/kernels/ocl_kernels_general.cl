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


__kernel void kernel_red_recurse(__global       IndexType *dst,
                                 __global const IndexType *src,
                                          const IndexType numElems) {

  IndexType index = BLOCK_SIZE * get_global_id(0);

  if (index >= numElems)
    return;

  IndexType i = index;

  if (i < BLOCK_SIZE)
    return;

  IndexType a = 0;

  while (i >= BLOCK_SIZE) {
    a += src[i];
    i -= BLOCK_SIZE;
  }

  dst[index] = a;

}

__kernel void kernel_red_partial_sum(__global       IndexType *dst,
                                     __global const IndexType *src,
                                              const IndexType numElems,
                                              const IndexType shift) {

  IndexType index = get_global_id(0);
  IndexType tid   = get_local_id(0);
  IndexType gid   = get_group_id(0);

  if (index < numElems) {

    __local IndexType data[BLOCK_SIZE];

    data[tid] = src[index];

    barrier(CLK_LOCAL_MEM_FENCE);

    for (IndexType i = BLOCK_SIZE/2; i > 0; i/=2) {

      if (tid < i)
        data[tid] = data[tid] + data[tid+i];

      barrier(CLK_LOCAL_MEM_FENCE);

    }

    if (tid == 0 && BLOCK_SIZE*(1+gid)-1<numElems)
      dst[BLOCK_SIZE*(1+gid)-1+shift] = data[0];

  }

}

__kernel void kernel_red_extrapolate(__global       IndexType *dst,
                                     __global const IndexType *srcBorder,
                                     __global const IndexType *srcData,
                                                    IndexType  numElems,
                                              const IndexType  shift) {

  IndexType index = get_local_size(0) * get_local_id(0);

  if (index < numElems-1) {

    IndexType sum = srcBorder[index];

    for(IndexType i = 0; i < get_local_size(0) && index+i<numElems; ++i) {
      sum += srcData[index+i];
      dst[index+i+shift] = sum;
    }

  }

}
