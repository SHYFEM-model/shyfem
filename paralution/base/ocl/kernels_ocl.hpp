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


#ifndef PARALUTION_KERNELS_OCL_HPP_
#define PARALUTION_KERNELS_OCL_HPP_

#define KERNELCOUNT 74

#if defined(__APPLE__) && defined(__MACH__)
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

static const std::string kernels_ocl[KERNELCOUNT] = {
                "kernel_reverse_index",
                "kernel_buffer_addscalar",
                "kernel_scaleadd",
                "kernel_scaleaddscale",
                "kernel_scaleadd2",
                "kernel_pointwisemult",
                "kernel_pointwisemult2",
                "kernel_copy_offset_from",
                "kernel_permute",
                "kernel_permute_backward",
                "kernel_dot",
                "kernel_axpy",
                "kernel_csr_spmv_scalar",
                "kernel_csr_add_spmv_scalar",
                "kernel_csr_scale_diagonal",
                "kernel_csr_scale_offdiagonal",
                "kernel_csr_add_diagonal",
                "kernel_csr_add_offdiagonal",
                "kernel_csr_extract_diag",
                "kernel_csr_extract_inv_diag",
                "kernel_csr_extract_submatrix_row_nnz",
                "kernel_csr_extract_submatrix_copy",
                "kernel_csr_diagmatmult_r",
                "kernel_csr_add_csr_same_struct",
                "kernel_scale",
                "kernel_mcsr_spmv_scalar",
                "kernel_mcsr_add_spmv_scalar",
                "kernel_ell_spmv",
                "kernel_ell_add_spmv",
                "kernel_dia_spmv",
                "kernel_dia_add_spmv",
                "kernel_coo_permute",
                "kernel_coo_spmv_flat",
                "kernel_coo_spmv_reduce_update",
                "kernel_coo_spmv_serial",
                "kernel_red_recurse",
                "kernel_red_partial_sum",
                "kernel_red_extrapolate",
                "kernel_csr_permute_rows",
                "kernel_csr_permute_cols",
                "kernel_csr_calc_row_nnz",
                "kernel_csr_permute_row_nnz",
                "kernel_reduce",
                "kernel_ell_max_row",
                "kernel_ell_csr_to_ell",
                "kernel_asum",
                "kernel_amax",
                "kernel_dense_spmv",
                "kernel_csr_extract_l_triangular",
                "kernel_csr_slower_nnz_per_row",
                "kernel_csr_supper_nnz_per_row",
                "kernel_csr_lower_nnz_per_row",
                "kernel_csr_upper_nnz_per_row",
                "kernel_csr_compress_count_nrow",
                "kernel_csr_compress_copy",
                "kernel_scaleaddscale_offset",
                "kernel_csr_extract_u_triangular",
                "kernel_norm",
                "kernel_dotc",
                "kernel_csr_diagmatmult_l",
                "kernel_power",
                "kernel_ell_nnz_coo",
                "kernel_ell_fill_ell",
                "kernel_ell_fill_coo",
                "kernel_copy_from_float",
                "kernel_copy_from_double",
                "kernel_dense_replace_column_vector",
                "kernel_dense_replace_row_vector",
                "kernel_dense_extract_column_vector",
                "kernel_dense_extract_row_vector",
                "kernel_csr_extract_column_vector",
                "kernel_csr_replace_column_vector_offset",
                "kernel_csr_replace_column_vector",
                "kernel_coo_csr_to_coo"
};

#define CL_KERNEL_REVERSE_INDEX                     paralution_get_kernel_ocl<ValueType>(  0)
#define CL_KERNEL_BUFFER_ADDSCALAR                  paralution_get_kernel_ocl<ValueType>(  4)
#define CL_KERNEL_SCALEADD                          paralution_get_kernel_ocl<ValueType>(  8)
#define CL_KERNEL_SCALEADDSCALE                     paralution_get_kernel_ocl<ValueType>( 12)
#define CL_KERNEL_SCALEADD2                         paralution_get_kernel_ocl<ValueType>( 16)
#define CL_KERNEL_POINTWISEMULT                     paralution_get_kernel_ocl<ValueType>( 20)
#define CL_KERNEL_POINTWISEMULT2                    paralution_get_kernel_ocl<ValueType>( 24)
#define CL_KERNEL_COPY_OFFSET_FROM                  paralution_get_kernel_ocl<ValueType>( 28)
#define CL_KERNEL_PERMUTE                           paralution_get_kernel_ocl<ValueType>( 32)
#define CL_KERNEL_PERMUTE_BACKWARD                  paralution_get_kernel_ocl<ValueType>( 36)
#define CL_KERNEL_DOT                               paralution_get_kernel_ocl<ValueType>( 40)
#define CL_KERNEL_AXPY                              paralution_get_kernel_ocl<ValueType>( 44)
#define CL_KERNEL_CSR_SPMV_SCALAR                   paralution_get_kernel_ocl<ValueType>( 48)
#define CL_KERNEL_CSR_ADD_SPMV_SCALAR               paralution_get_kernel_ocl<ValueType>( 52)
#define CL_KERNEL_CSR_SCALE_DIAGONAL                paralution_get_kernel_ocl<ValueType>( 56)
#define CL_KERNEL_CSR_SCALE_OFFDIAGONAL             paralution_get_kernel_ocl<ValueType>( 60)
#define CL_KERNEL_CSR_ADD_DIAGONAL                  paralution_get_kernel_ocl<ValueType>( 64)
#define CL_KERNEL_CSR_ADD_OFFDIAGONAL               paralution_get_kernel_ocl<ValueType>( 68)
#define CL_KERNEL_CSR_EXTRACT_DIAG                  paralution_get_kernel_ocl<ValueType>( 72)
#define CL_KERNEL_CSR_EXTRACT_INV_DIAG              paralution_get_kernel_ocl<ValueType>( 76)
#define CL_KERNEL_CSR_EXTRACT_SUBMATRIX_ROW_NNZ     paralution_get_kernel_ocl<ValueType>( 80)
#define CL_KERNEL_CSR_EXTRACT_SUBMATRIX_COPY        paralution_get_kernel_ocl<ValueType>( 84)
#define CL_KERNEL_CSR_DIAGMATMULT_R                 paralution_get_kernel_ocl<ValueType>( 88)
#define CL_KERNEL_CSR_ADD_CSR_SAME_STRUCT           paralution_get_kernel_ocl<ValueType>( 92)
#define CL_KERNEL_SCALE                             paralution_get_kernel_ocl<ValueType>( 96)
#define CL_KERNEL_MCSR_SPMV_SCALAR                  paralution_get_kernel_ocl<ValueType>(100)
#define CL_KERNEL_MCSR_ADD_SPMV_SCALAR              paralution_get_kernel_ocl<ValueType>(104)
#define CL_KERNEL_ELL_SPMV                          paralution_get_kernel_ocl<ValueType>(108)
#define CL_KERNEL_ELL_ADD_SPMV                      paralution_get_kernel_ocl<ValueType>(112)
#define CL_KERNEL_DIA_SPMV                          paralution_get_kernel_ocl<ValueType>(116)
#define CL_KERNEL_DIA_ADD_SPMV                      paralution_get_kernel_ocl<ValueType>(120)
#define CL_KERNEL_COO_PERMUTE                       paralution_get_kernel_ocl<ValueType>(124)
#define CL_KERNEL_COO_SPMV_FLAT                     paralution_get_kernel_ocl<ValueType>(128)
#define CL_KERNEL_COO_SPMV_REDUCE_UPDATE            paralution_get_kernel_ocl<ValueType>(132)
#define CL_KERNEL_COO_SPMV_SERIAL                   paralution_get_kernel_ocl<ValueType>(136)
#define CL_KERNEL_RED_RECURSE                       paralution_get_kernel_ocl<ValueType>(140)
#define CL_KERNEL_RED_PARTIAL_SUM                   paralution_get_kernel_ocl<ValueType>(144)
#define CL_KERNEL_RED_EXTRAPOLATE                   paralution_get_kernel_ocl<ValueType>(148)
#define CL_KERNEL_CSR_PERMUTE_ROWS                  paralution_get_kernel_ocl<ValueType>(152)
#define CL_KERNEL_CSR_PERMUTE_COLS                  paralution_get_kernel_ocl<ValueType>(156)
#define CL_KERNEL_CSR_CALC_ROW_NNZ                  paralution_get_kernel_ocl<ValueType>(160)
#define CL_KERNEL_CSR_PERMUTE_ROW_NNZ               paralution_get_kernel_ocl<ValueType>(164)
#define CL_KERNEL_REDUCE                            paralution_get_kernel_ocl<ValueType>(168)
#define CL_KERNEL_ELL_MAX_ROW                       paralution_get_kernel_ocl<ValueType>(172)
#define CL_KERNEL_ELL_CSR_TO_ELL                    paralution_get_kernel_ocl<ValueType>(176)
#define CL_KERNEL_ASUM                              paralution_get_kernel_ocl<ValueType>(180)
#define CL_KERNEL_AMAX                              paralution_get_kernel_ocl<ValueType>(184)
#define CL_KERNEL_DENSE_SPMV                        paralution_get_kernel_ocl<ValueType>(188)
#define CL_KERNEL_CSR_EXTRACT_L_TRIANGULAR          paralution_get_kernel_ocl<ValueType>(192)
#define CL_KERNEL_CSR_SLOWER_NNZ_PER_ROW            paralution_get_kernel_ocl<ValueType>(196)
#define CL_KERNEL_CSR_SUPPER_NNZ_PER_ROW            paralution_get_kernel_ocl<ValueType>(200)
#define CL_KERNEL_CSR_LOWER_NNZ_PER_ROW             paralution_get_kernel_ocl<ValueType>(204)
#define CL_KERNEL_CSR_UPPER_NNZ_PER_ROW             paralution_get_kernel_ocl<ValueType>(208)
#define CL_KERNEL_CSR_COMPRESS_COUNT_NROW           paralution_get_kernel_ocl<ValueType>(212)
#define CL_KERNEL_CSR_COMPRESS_COPY                 paralution_get_kernel_ocl<ValueType>(216)
#define CL_KERNEL_SCALEADDSCALE_OFFSET              paralution_get_kernel_ocl<ValueType>(220)
#define CL_KERNEL_CSR_EXTRACT_U_TRIANGULAR          paralution_get_kernel_ocl<ValueType>(224)
#define CL_KERNEL_NORM                              paralution_get_kernel_ocl<ValueType>(228)
#define CL_KERNEL_DOTC                              paralution_get_kernel_ocl<ValueType>(232)
#define CL_KERNEL_CSR_DIAGMATMULT_L                 paralution_get_kernel_ocl<ValueType>(236)
#define CL_KERNEL_POWER                             paralution_get_kernel_ocl<ValueType>(240)
#define CL_KERNEL_ELL_NNZ_COO                       paralution_get_kernel_ocl<ValueType>(244)
#define CL_KERNEL_ELL_FILL_ELL                      paralution_get_kernel_ocl<ValueType>(248)
#define CL_KERNEL_ELL_FILL_COO                      paralution_get_kernel_ocl<ValueType>(252)
#define CL_KERNEL_COPY_FROM_FLOAT                   paralution_get_kernel_ocl<double>   (256)
#define CL_KERNEL_COPY_FROM_DOUBLE                  paralution_get_kernel_ocl<float>    (260)
#define CL_KERNEL_DENSE_REPLACE_COLUMN_VECTOR       paralution_get_kernel_ocl<ValueType>(264)
#define CL_KERNEL_DENSE_REPLACE_ROW_VECTOR          paralution_get_kernel_ocl<ValueType>(268)
#define CL_KERNEL_DENSE_EXTRACT_COLUMN_VECTOR       paralution_get_kernel_ocl<ValueType>(272)
#define CL_KERNEL_DENSE_EXTRACT_ROW_VECTOR          paralution_get_kernel_ocl<ValueType>(276)
#define CL_KERNEL_CSR_EXTRACT_COLUMN_VECTOR         paralution_get_kernel_ocl<ValueType>(280)
#define CL_KERNEL_CSR_REPLACE_COLUMN_VECTOR_OFFSET  paralution_get_kernel_ocl<ValueType>(284)
#define CL_KERNEL_CSR_REPLACE_COLUMN_VECTOR         paralution_get_kernel_ocl<ValueType>(288)
#define CL_KERNEL_COO_CSR_TO_COO                    paralution_get_kernel_ocl<ValueType>(292)

template <typename ValueType, typename A0, typename A1, typename A2>
cl_int ocl_kernel(cl_kernel kernel, cl_command_queue cmdQueue, const size_t LocalSize, const size_t GlobalSize, A0 arg0, A1 arg1, A2 arg2, bool sync = true) {

  cl_int err;
  cl_event event;

  err  = clSetKernelArg(kernel, 0, sizeof(A0), (void*) &arg0);
  err |= clSetKernelArg(kernel, 1, sizeof(A1), (void*) &arg1);
  err |= clSetKernelArg(kernel, 2, sizeof(A2), (void*) &arg2);

  if (err != CL_SUCCESS)
    return err;

  size_t szLocalWorkSize[1];
  size_t szGlobalWorkSize[1];

  szLocalWorkSize[0]  = LocalSize;
  szGlobalWorkSize[0] = GlobalSize;

  err = clEnqueueNDRangeKernel(cmdQueue, kernel, 1, NULL, szGlobalWorkSize, szLocalWorkSize, 0, NULL, &event);

  if (err != CL_SUCCESS)
    return err;

  if (sync == true) {

    // Wait for kernel run to finish
    err = clWaitForEvents(1, &event);

    if (err != CL_SUCCESS)
      return err;

    // Release event when kernel run finished
    err = clReleaseEvent(event);

    if (err != CL_SUCCESS)
      return err;

  }

  return CL_SUCCESS;

}

template <typename ValueType, typename A0, typename A1, typename A2, typename A3>
cl_int ocl_kernel(cl_kernel kernel, cl_command_queue cmdQueue, const size_t LocalSize, const size_t GlobalSize, A0 arg0, A1 arg1, A2 arg2, A3 arg3, bool sync = true) {

  cl_int err;
  cl_event event;

  err  = clSetKernelArg(kernel, 0, sizeof(A0), (void*) &arg0);
  err |= clSetKernelArg(kernel, 1, sizeof(A1), (void*) &arg1);
  err |= clSetKernelArg(kernel, 2, sizeof(A2), (void*) &arg2);
  err |= clSetKernelArg(kernel, 3, sizeof(A3), (void*) &arg3);

  if (err != CL_SUCCESS)
    return err;

  size_t szLocalWorkSize[1];
  size_t szGlobalWorkSize[1];

  szLocalWorkSize[0]  = LocalSize;
  szGlobalWorkSize[0] = GlobalSize;

  err = clEnqueueNDRangeKernel(cmdQueue, kernel, 1, NULL, szGlobalWorkSize, szLocalWorkSize, 0, NULL, &event);

  if (err != CL_SUCCESS)
    return err;

  if (sync == true) {

    // Wait for kernel run to finish
    err = clWaitForEvents(1, &event);

    if (err != CL_SUCCESS)
      return err;

    // Release event when kernel run finished
    err = clReleaseEvent(event);

    if (err != CL_SUCCESS)
      return err;

  }

  return CL_SUCCESS;

}

template <typename ValueType, typename A0, typename A1, typename A2, typename A3, typename A4>
cl_int ocl_kernel(cl_kernel kernel, cl_command_queue cmdQueue, const size_t LocalSize, const size_t GlobalSize, A0 arg0, A1 arg1, A2 arg2, A3 arg3, A4 arg4, bool sync = true) {

  cl_int err;
  cl_event event;

  err  = clSetKernelArg(kernel, 0, sizeof(A0), (void*) &arg0);
  err |= clSetKernelArg(kernel, 1, sizeof(A1), (void*) &arg1);
  err |= clSetKernelArg(kernel, 2, sizeof(A2), (void*) &arg2);
  err |= clSetKernelArg(kernel, 3, sizeof(A3), (void*) &arg3);
  err |= clSetKernelArg(kernel, 4, sizeof(A4), (void*) &arg4);

  if (err != CL_SUCCESS)
    return err;

  size_t szLocalWorkSize[1];
  size_t szGlobalWorkSize[1];

  szLocalWorkSize[0]  = LocalSize;
  szGlobalWorkSize[0] = GlobalSize;

  err = clEnqueueNDRangeKernel(cmdQueue, kernel, 1, NULL, szGlobalWorkSize, szLocalWorkSize, 0, NULL, &event);

  if (err != CL_SUCCESS)
    return err;

  if (sync == true) {

    // Wait for kernel run to finish
    err = clWaitForEvents(1, &event);

    if (err != CL_SUCCESS)
      return err;

    // Release event when kernel run finished
    err = clReleaseEvent(event);

    if (err != CL_SUCCESS)
      return err;

  }

  return CL_SUCCESS;

}

template <typename ValueType, typename A0, typename A1, typename A2, typename A3, typename A4, typename A5>
cl_int ocl_kernel(cl_kernel kernel, cl_command_queue cmdQueue, const size_t LocalSize, const size_t GlobalSize, A0 arg0, A1 arg1, A2 arg2, A3 arg3, A4 arg4, A5 arg5, bool sync = true) {

  cl_int err;
  cl_event event;

  err  = clSetKernelArg(kernel, 0, sizeof(A0), (void*) &arg0);
  err |= clSetKernelArg(kernel, 1, sizeof(A1), (void*) &arg1);
  err |= clSetKernelArg(kernel, 2, sizeof(A2), (void*) &arg2);
  err |= clSetKernelArg(kernel, 3, sizeof(A3), (void*) &arg3);
  err |= clSetKernelArg(kernel, 4, sizeof(A4), (void*) &arg4);
  err |= clSetKernelArg(kernel, 5, sizeof(A5), (void*) &arg5);

  if (err != CL_SUCCESS)
    return err;

  size_t szLocalWorkSize[1];
  size_t szGlobalWorkSize[1];

  szLocalWorkSize[0]  = LocalSize;
  szGlobalWorkSize[0] = GlobalSize;

  err = clEnqueueNDRangeKernel(cmdQueue, kernel, 1, NULL, szGlobalWorkSize, szLocalWorkSize, 0, NULL, &event);

  if (err != CL_SUCCESS)
    return err;

  if (sync == true) {

    // Wait for kernel run to finish
    err = clWaitForEvents(1, &event);

    if (err != CL_SUCCESS)
      return err;

    // Release event when kernel run finished
    err = clReleaseEvent(event);

    if (err != CL_SUCCESS)
      return err;

  }

  return CL_SUCCESS;

}

template <typename ValueType, typename A0, typename A1, typename A2, typename A3, typename A4, typename A5, typename A6>
cl_int ocl_kernel(cl_kernel kernel, cl_command_queue cmdQueue, const size_t LocalSize, const size_t GlobalSize, A0 arg0, A1 arg1, A2 arg2, A3 arg3, A4 arg4, A5 arg5, A6 arg6, bool sync = true) {

  cl_int err;
  cl_event event;

  err  = clSetKernelArg(kernel, 0, sizeof(A0), (void*) &arg0);
  err |= clSetKernelArg(kernel, 1, sizeof(A1), (void*) &arg1);
  err |= clSetKernelArg(kernel, 2, sizeof(A2), (void*) &arg2);
  err |= clSetKernelArg(kernel, 3, sizeof(A3), (void*) &arg3);
  err |= clSetKernelArg(kernel, 4, sizeof(A4), (void*) &arg4);
  err |= clSetKernelArg(kernel, 5, sizeof(A5), (void*) &arg5);
  err |= clSetKernelArg(kernel, 6, sizeof(A6), (void*) &arg6);

  if (err != CL_SUCCESS)
    return err;

  size_t szLocalWorkSize[1];
  size_t szGlobalWorkSize[1];

  szLocalWorkSize[0]  = LocalSize;
  szGlobalWorkSize[0] = GlobalSize;

  err = clEnqueueNDRangeKernel(cmdQueue, kernel, 1, NULL, szGlobalWorkSize, szLocalWorkSize, 0, NULL, &event);

  if (err != CL_SUCCESS)
    return err;

  if (sync == true) {

    // Wait for kernel run to finish
    err = clWaitForEvents(1, &event);

    if (err != CL_SUCCESS)
      return err;

    // Release event when kernel run finished
    err = clReleaseEvent(event);

    if (err != CL_SUCCESS)
      return err;

  }

  return CL_SUCCESS;

}

template <typename ValueType, typename A0, typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7>
cl_int ocl_kernel(cl_kernel kernel, cl_command_queue cmdQueue, const size_t LocalSize, const size_t GlobalSize, A0 arg0, A1 arg1, A2 arg2, A3 arg3, A4 arg4, A5 arg5, A6 arg6, A7 arg7, bool sync = true) {

  cl_int err;
  cl_event event;

  err  = clSetKernelArg(kernel, 0, sizeof(A0), (void*) &arg0);
  err |= clSetKernelArg(kernel, 1, sizeof(A1), (void*) &arg1);
  err |= clSetKernelArg(kernel, 2, sizeof(A2), (void*) &arg2);
  err |= clSetKernelArg(kernel, 3, sizeof(A3), (void*) &arg3);
  err |= clSetKernelArg(kernel, 4, sizeof(A4), (void*) &arg4);
  err |= clSetKernelArg(kernel, 5, sizeof(A5), (void*) &arg5);
  err |= clSetKernelArg(kernel, 6, sizeof(A6), (void*) &arg6);
  err |= clSetKernelArg(kernel, 7, sizeof(A7), (void*) &arg7);

  if (err != CL_SUCCESS)
    return err;

  size_t szLocalWorkSize[1];
  size_t szGlobalWorkSize[1];

  szLocalWorkSize[0]  = LocalSize;
  szGlobalWorkSize[0] = GlobalSize;

  err = clEnqueueNDRangeKernel(cmdQueue, kernel, 1, NULL, szGlobalWorkSize, szLocalWorkSize, 0, NULL, &event);

  if (err != CL_SUCCESS)
    return err;

  if (sync == true) {

    // Wait for kernel run to finish
    err = clWaitForEvents(1, &event);

    if (err != CL_SUCCESS)
      return err;

    // Release event when kernel run finished
    err = clReleaseEvent(event);

    if (err != CL_SUCCESS)
      return err;

  }

  return CL_SUCCESS;

}

template <typename ValueType, typename A0, typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7, typename A8>
cl_int ocl_kernel(cl_kernel kernel, cl_command_queue cmdQueue, const size_t LocalSize, const size_t GlobalSize, A0 arg0, A1 arg1, A2 arg2, A3 arg3, A4 arg4, A5 arg5, A6 arg6, A7 arg7, A8 arg8, bool sync = true) {

  cl_int err;
  cl_event event;

  err  = clSetKernelArg(kernel, 0, sizeof(A0), (void*) &arg0);
  err |= clSetKernelArg(kernel, 1, sizeof(A1), (void*) &arg1);
  err |= clSetKernelArg(kernel, 2, sizeof(A2), (void*) &arg2);
  err |= clSetKernelArg(kernel, 3, sizeof(A3), (void*) &arg3);
  err |= clSetKernelArg(kernel, 4, sizeof(A4), (void*) &arg4);
  err |= clSetKernelArg(kernel, 5, sizeof(A5), (void*) &arg5);
  err |= clSetKernelArg(kernel, 6, sizeof(A6), (void*) &arg6);
  err |= clSetKernelArg(kernel, 7, sizeof(A7), (void*) &arg7);
  err |= clSetKernelArg(kernel, 8, sizeof(A8), (void*) &arg8);

  if (err != CL_SUCCESS)
    return err;

  size_t szLocalWorkSize[1];
  size_t szGlobalWorkSize[1];

  szLocalWorkSize[0]  = LocalSize;
  szGlobalWorkSize[0] = GlobalSize;

  err = clEnqueueNDRangeKernel(cmdQueue, kernel, 1, NULL, szGlobalWorkSize, szLocalWorkSize, 0, NULL, &event);

  if (err != CL_SUCCESS)
    return err;

  if (sync == true) {

    // Wait for kernel run to finish
    err = clWaitForEvents(1, &event);

    if (err != CL_SUCCESS)
      return err;

    // Release event when kernel run finished
    err = clReleaseEvent(event);

    if (err != CL_SUCCESS)
      return err;

  }

  return CL_SUCCESS;

}

template <typename ValueType, typename A0, typename A1, typename A2, typename A3, typename A4, typename A5, typename A6, typename A7, typename A8, typename A9>
cl_int ocl_kernel(cl_kernel kernel, cl_command_queue cmdQueue, const size_t LocalSize, const size_t GlobalSize, A0 arg0, A1 arg1, A2 arg2, A3 arg3, A4 arg4, A5 arg5, A6 arg6, A7 arg7, A8 arg8, A9 arg9, bool sync = true) {

  cl_int err;
  cl_event event;

  err  = clSetKernelArg(kernel, 0, sizeof(A0), (void*) &arg0);
  err |= clSetKernelArg(kernel, 1, sizeof(A1), (void*) &arg1);
  err |= clSetKernelArg(kernel, 2, sizeof(A2), (void*) &arg2);
  err |= clSetKernelArg(kernel, 3, sizeof(A3), (void*) &arg3);
  err |= clSetKernelArg(kernel, 4, sizeof(A4), (void*) &arg4);
  err |= clSetKernelArg(kernel, 5, sizeof(A5), (void*) &arg5);
  err |= clSetKernelArg(kernel, 6, sizeof(A6), (void*) &arg6);
  err |= clSetKernelArg(kernel, 7, sizeof(A7), (void*) &arg7);
  err |= clSetKernelArg(kernel, 8, sizeof(A8), (void*) &arg8);
  err |= clSetKernelArg(kernel, 9, sizeof(A9), (void*) &arg9);

  if (err != CL_SUCCESS)
    return err;

  size_t szLocalWorkSize[1];
  size_t szGlobalWorkSize[1];

  szLocalWorkSize[0]  = LocalSize;
  szGlobalWorkSize[0] = GlobalSize;

  err = clEnqueueNDRangeKernel(cmdQueue, kernel, 1, NULL, szGlobalWorkSize, szLocalWorkSize, 0, NULL, &event);

  if (err != CL_SUCCESS)
    return err;

  if (sync == true) {

    // Wait for kernel run to finish
    err = clWaitForEvents(1, &event);

    if (err != CL_SUCCESS)
      return err;

    // Release event when kernel run finished
    err = clReleaseEvent(event);

    if (err != CL_SUCCESS)
      return err;

  }

  return CL_SUCCESS;

}

#endif // PARALUTION_KERNELS_OCL_HPP_
