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


#include "../../utils/def.hpp"
#include "../../utils/log.hpp"
#include "kernels_ocl.hpp"
#include "ocl_allocate_free.hpp"
#include "ocl_utils.hpp"
#include "ocl_matrix_csr.hpp"
#include "ocl_matrix_coo.hpp"
#include "ocl_vector.hpp"
#include "../host/host_matrix_coo.hpp"
#include "../backend_manager.hpp"

#include <algorithm>

namespace paralution {

template <typename ValueType>
OCLAcceleratorMatrixCOO<ValueType>::OCLAcceleratorMatrixCOO() {

  // no default constructors
  LOG_INFO("no default constructor");
  FATAL_ERROR(__FILE__, __LINE__);

}

template <typename ValueType>
OCLAcceleratorMatrixCOO<ValueType>::OCLAcceleratorMatrixCOO(const Paralution_Backend_Descriptor local_backend) {

  LOG_DEBUG(this, "OCLAcceleratorMatrixCOO::OCLAcceleratorMatrixCOO()",
            "constructor with local_backend");

  this->mat_.row = NULL;
  this->mat_.col = NULL;
  this->mat_.val = NULL;

  this->set_backend(local_backend);

}

template <typename ValueType>
OCLAcceleratorMatrixCOO<ValueType>::~OCLAcceleratorMatrixCOO() {

  LOG_DEBUG(this, "OCLAcceleratorMatrixCOO::~OCLAcceleratorMatrixCOO()",
            "destructor");

  this->Clear();

}

template <typename ValueType>
void OCLAcceleratorMatrixCOO<ValueType>::info(void) const {

  LOG_INFO("OCLAcceleratorMatrixCOO<ValueType>");

}

template <typename ValueType>
void OCLAcceleratorMatrixCOO<ValueType>::AllocateCOO(const int nnz, const int nrow, const int ncol) {

  assert (nnz  >= 0);
  assert (ncol >= 0);
  assert (nrow >= 0);

  if (this->nnz_ > 0)
    this->Clear();

  if (nnz > 0) {

    allocate_ocl(nnz, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &this->mat_.row);
    allocate_ocl(nnz, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &this->mat_.col);
    allocate_ocl(nnz, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &this->mat_.val);

    // Set entries of device object to zero
    ocl_set_to(nnz, (int) 0, this->mat_.row, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    ocl_set_to(nnz, (int) 0, this->mat_.col, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    ocl_set_to(nnz, (ValueType) 0, this->mat_.val, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    this->nrow_ = nrow;
    this->ncol_ = ncol;
    this->nnz_  = nnz;

  }

}

template <typename ValueType>
void OCLAcceleratorMatrixCOO<ValueType>::SetDataPtrCOO(int **row, int **col, ValueType **val,
                                                       const int nnz, const int nrow, const int ncol) {

  assert (*row != NULL);
  assert (*col != NULL);
  assert (*val != NULL);
  assert (nnz  > 0);
  assert (nrow > 0);
  assert (ncol > 0);

  this->Clear();

  this->nrow_ = nrow;
  this->ncol_ = ncol;
  this->nnz_  = nnz;

  cl_int err = clFinish(OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
  CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  this->mat_.row = *row;
  this->mat_.col = *col;
  this->mat_.val = *val;

}

template <typename ValueType>
void OCLAcceleratorMatrixCOO<ValueType>::LeaveDataPtrCOO(int **row, int **col, ValueType **val) {

  assert (this->nrow_ > 0);
  assert (this->ncol_ > 0);
  assert (this->nnz_  > 0);

  cl_int err = clFinish(OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
  CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  // see free_host function for details
  *row = this->mat_.row;
  *col = this->mat_.col;
  *val = this->mat_.val;

  this->mat_.row = NULL;
  this->mat_.col = NULL;
  this->mat_.val = NULL;

  this->nrow_ = 0;
  this->ncol_ = 0;
  this->nnz_  = 0;

}

template <typename ValueType>
void OCLAcceleratorMatrixCOO<ValueType>::Clear(void) {

  if (this->nnz_ > 0) {

    free_ocl(&this->mat_.row);
    free_ocl(&this->mat_.col);
    free_ocl(&this->mat_.val);

    this->nrow_ = 0;
    this->ncol_ = 0;
    this->nnz_  = 0;

  }

}

template <typename ValueType>
void OCLAcceleratorMatrixCOO<ValueType>::CopyFromHost(const HostMatrix<ValueType> &src) {

  assert (&src != NULL);

  const HostMatrixCOO<ValueType> *cast_mat;

  // copy only in the same format
  assert (this->get_mat_format() == src.get_mat_format());

  // CPU to OCL copy
  if ((cast_mat = dynamic_cast<const HostMatrixCOO<ValueType>*> (&src)) != NULL) {

    if (this->nnz_ == 0)
      this->AllocateCOO(cast_mat->nnz_, cast_mat->nrow_, cast_mat->ncol_);

    if (this->nnz_ > 0) {

      assert (this->nnz_  == cast_mat->nnz_);
      assert (this->nrow_ == cast_mat->nrow_);
      assert (this->ncol_ == cast_mat->ncol_);

      // Copy object from host to device memory
      ocl_host2dev(this->nnz_,         // size
                   cast_mat->mat_.row, // src
                   this->mat_.row,     // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      // Copy object from host to device memory
      ocl_host2dev(this->nnz_,         // size
                   cast_mat->mat_.col, // src
                   this->mat_.col,     // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      // Copy object from host to device memory
      ocl_host2dev(this->nnz_,         // size
                   cast_mat->mat_.val, // src
                   this->mat_.val,     // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    }

  } else {

    LOG_INFO("Error unsupported OCL matrix type");
    this->info();
    src.info();
    FATAL_ERROR(__FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorMatrixCOO<ValueType>::CopyToHost(HostMatrix<ValueType> *dst) const {

  assert (dst != NULL);

  HostMatrixCOO<ValueType> *cast_mat;

  // copy only in the same format
  assert (this->get_mat_format() == dst->get_mat_format());

  // OCL to CPU copy
  if ((cast_mat = dynamic_cast<HostMatrixCOO<ValueType>*> (dst)) != NULL) {

    cast_mat->set_backend(this->local_backend_);

    if (cast_mat->nnz_ == 0)
      cast_mat->AllocateCOO(this->nnz_, this->nrow_, this->ncol_);

    if (this->nnz_ > 0) {

      assert (this->nnz_  == cast_mat->nnz_);
      assert (this->nrow_ == cast_mat->nrow_);
      assert (this->ncol_ == cast_mat->ncol_);

      // Copy object from device to host memory
      ocl_dev2host(this->nnz_,         // size
                   this->mat_.row,     // src
                   cast_mat->mat_.row, // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      // Copy object from device to host memory
      ocl_dev2host(this->nnz_,         // size
                   this->mat_.col,     // src
                   cast_mat->mat_.col, // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      // Copy object from device to host memory
      ocl_dev2host(this->nnz_,         // size
                   this->mat_.val,     // src
                   cast_mat->mat_.val, // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    }

  } else {

    LOG_INFO("Error unsupported OCL matrix type");
    this->info();
    dst->info();
    FATAL_ERROR(__FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorMatrixCOO<ValueType>::CopyFrom(const BaseMatrix<ValueType> &src) {

  assert (&src != NULL);

  const OCLAcceleratorMatrixCOO<ValueType> *ocl_cast_mat;
  const HostMatrix<ValueType> *host_cast_mat;

  // copy only in the same format
  assert (this->get_mat_format() == src.get_mat_format());

  // OCL to OCL copy
  if ((ocl_cast_mat = dynamic_cast<const OCLAcceleratorMatrixCOO<ValueType>*> (&src)) != NULL) {

    if (this->nnz_ == 0)
      this->AllocateCOO(ocl_cast_mat->nnz_, ocl_cast_mat->nrow_, ocl_cast_mat->ncol_);

    assert (this->nnz_  == ocl_cast_mat->nnz_);
    assert (this->nrow_ == ocl_cast_mat->nrow_);
    assert (this->ncol_ == ocl_cast_mat->ncol_);

    if (this->nnz_ > 0) {

      // Copy object from device to device memory (internal copy)
      ocl_dev2dev(this->nnz_,             // size
                  ocl_cast_mat->mat_.row, // src
                  this->mat_.row,         // dst
                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      // Copy object from device to device memory (internal copy)
      ocl_dev2dev(this->nnz_,             // size
                  ocl_cast_mat->mat_.col, // src
                  this->mat_.col,         // dst
                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      // Copy object from device to device memory (internal copy)
      ocl_dev2dev(this->nnz_,             // size
                  ocl_cast_mat->mat_.val, // src
                  this->mat_.val,         // dst
                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    }

  } else {

    //CPU to OCL
    if ((host_cast_mat = dynamic_cast<const HostMatrix<ValueType>*> (&src)) != NULL) {

      this->CopyFromHost(*host_cast_mat);

    } else {

      LOG_INFO("Error unsupported OCL matrix type");
      this->info();
      src.info();
      FATAL_ERROR(__FILE__, __LINE__);

    }

  }

}

template <typename ValueType>
void OCLAcceleratorMatrixCOO<ValueType>::CopyTo(BaseMatrix<ValueType> *dst) const {

  assert (dst != NULL);

  OCLAcceleratorMatrixCOO<ValueType> *ocl_cast_mat;
  HostMatrix<ValueType> *host_cast_mat;

  // copy only in the same format
  assert (this->get_mat_format() == dst->get_mat_format());

  // OCL to OCL copy
  if ((ocl_cast_mat = dynamic_cast<OCLAcceleratorMatrixCOO<ValueType>*> (dst)) != NULL) {

    ocl_cast_mat->set_backend(this->local_backend_);

    if (ocl_cast_mat->nnz_ == 0)
      ocl_cast_mat->AllocateCOO(this->nnz_, this->nrow_, this->ncol_);

    assert (this->nnz_  == ocl_cast_mat->nnz_);
    assert (this->nrow_ == ocl_cast_mat->nrow_);
    assert (this->ncol_ == ocl_cast_mat->ncol_);

    if (this->nnz_ > 0) {

      // Copy object from device to device memory (internal copy)
      ocl_dev2dev(this->nnz_,             // size
                  this->mat_.row,         // src
                  ocl_cast_mat->mat_.row, // dst
                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      // Copy object from device to device memory (internal copy)
      ocl_dev2dev(this->nnz_,             // size
                  this->mat_.col,         // src
                  ocl_cast_mat->mat_.col, // dst
                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      // Copy object from device to device memory (internal copy)
      ocl_dev2dev(this->nnz_,             // size
                  this->mat_.val,         // src
                  ocl_cast_mat->mat_.val, // dst
                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    }

  } else {

    //OCL to CPU
    if ((host_cast_mat = dynamic_cast<HostMatrix<ValueType>*> (dst)) != NULL) {

      this->CopyToHost(host_cast_mat);

    } else {

      LOG_INFO("Error unsupported OCL matrix type");
      this->info();
      dst->info();
      FATAL_ERROR(__FILE__, __LINE__);

    }

  }

}

template <typename ValueType>
bool OCLAcceleratorMatrixCOO<ValueType>::ConvertFrom(const BaseMatrix<ValueType> &mat) {

  assert (&mat != NULL);

  this->Clear();

  // empty matrix is empty matrix
  if (mat.get_nnz() == 0)
    return true;

  const OCLAcceleratorMatrixCOO<ValueType> *cast_mat_coo;

  if ((cast_mat_coo = dynamic_cast<const OCLAcceleratorMatrixCOO<ValueType>*> (&mat)) != NULL) {

      this->CopyFrom(*cast_mat_coo);
      return true;

  }

  const OCLAcceleratorMatrixCSR<ValueType> *cast_mat_csr;
  if ((cast_mat_csr = dynamic_cast<const OCLAcceleratorMatrixCSR<ValueType>*> (&mat)) != NULL) {

    this->Clear();

    assert (cast_mat_csr->nrow_ > 0);
    assert (cast_mat_csr->ncol_ > 0);
    assert (cast_mat_csr->nnz_  > 0);

    this->AllocateCOO(cast_mat_csr->nnz_, cast_mat_csr->nrow_, cast_mat_csr->ncol_);

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = (cast_mat_csr->nrow_ / LocalSize + 1) * LocalSize;

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_COO_CSR_TO_COO,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       cast_mat_csr->nrow_, cast_mat_csr->mat_.row_offset, this->mat_.row);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    ocl_dev2dev(this->nnz_,             // size
                cast_mat_csr->mat_.col, // src
                this->mat_.col,         // dst
                OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    ocl_dev2dev(this->nnz_,             // size
                cast_mat_csr->mat_.val, // src
                this->mat_.val,         // dst
                OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    return true;

  }

  return false;

}

template <typename ValueType>
void OCLAcceleratorMatrixCOO<ValueType>::Apply(const BaseVector<ValueType> &in, BaseVector<ValueType> *out) const {

  // TODO some devices hang up while waiting for COO_SPMV_FLAT event to finish
  // this is a bug we are fixing in some future release
  if (this->nnz_ > 0) {

    assert (&in != NULL);
    assert (out != NULL);
    assert (in.  get_size() >= 0);
    assert (out->get_size() >= 0);
    assert (in.  get_size() == this->ncol_);
    assert (out->get_size() == this->nrow_);

    const OCLAcceleratorVector<ValueType> *cast_in = dynamic_cast<const OCLAcceleratorVector<ValueType>*> (&in);
    OCLAcceleratorVector<ValueType> *cast_out      = dynamic_cast<      OCLAcceleratorVector<ValueType>*> (out);

    assert (cast_in  != NULL);
    assert (cast_out != NULL);

    cast_out->Zeros();

    // do not support super small matrices
    assert (this->nnz_ > this->local_backend_.OCL_warp_size);

    // ----------------------------------------------------------
    // Modified and adapted from CUSP 0.3.1, 
    // http://code.google.com/p/cusp-library/
    // NVIDIA, APACHE LICENSE 2.0
    // ----------------------------------------------------------
    // see __spmv_coo_flat(...)
    // ----------------------------------------------------------
    // CHANGELOG
    // - adapted interface
    // ----------------------------------------------------------  

    //TODO
    //move in extra file -  max_active_blocks, warp_size, block_size

    const int BLOCK_SIZE = (int) this->local_backend_.OCL_max_work_group_size;
    //    const unsigned int MAX_BLOCKS = this->local_backend_.GPU_max_blocks;

    const int MAX_BLOCKS = 32; //  cusp::detail::device::arch::max_active_blocks(spmv_coo_flat_kernel<IndexType, ValueType, BLOCK_SIZE, UseCache>, BLOCK_SIZE, (size_t) 0);

    const int WARPS_PER_BLOCK = BLOCK_SIZE / this->local_backend_.OCL_warp_size;


    const int num_units = this->nnz_ / this->local_backend_.OCL_warp_size; 
    const unsigned int num_warps  = std::min(num_units, WARPS_PER_BLOCK * MAX_BLOCKS);
    const unsigned int num_blocks = (num_warps + (WARPS_PER_BLOCK-1)) / WARPS_PER_BLOCK; // (N + (granularity - 1)) / granularity
    const unsigned int num_iters  = (num_units +  (num_warps-1)) / num_warps;
    
    const unsigned int interval_size = this->local_backend_.OCL_warp_size * num_iters;

    const int tail = num_units * this->local_backend_.OCL_warp_size; // do the last few nonzeros separately (fewer than this->local_backend_.GPU_warp elements)

    const unsigned int active_warps = (interval_size == 0) ? 0 : ((tail + (interval_size-1))/interval_size);

    int *temp_rows = NULL;
    ValueType *temp_vals = NULL;

    allocate_ocl(active_warps, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &temp_rows);
    allocate_ocl(active_warps, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &temp_vals);

    //  LOG_INFO("active_warps = " << active_warps);
    //  LOG_INFO("tail =" << tail);
    //  LOG_INFO("interval_size =" << interval_size);
    //  LOG_INFO("num_iters =" << num_iters);
    //  LOG_INFO("num_blocks =" << num_blocks);
    //  LOG_INFO("num_warps =" << num_warps);
    //  LOG_INFO("num_units =" << num_units);
    //  LOG_INFO("WARPS_PER_BLOCK =" << WARPS_PER_BLOCK);
    //  LOG_INFO("MAX_BLOCKS =" << MAX_BLOCKS);
    //  LOG_INFO("BLOCK_SIZE =" << BLOCK_SIZE);
    //  LOG_INFO("WARP_SIZE =" << WARP_SIZE);
    //  LOG_INFO("WARP_SIZE =" << this->local_backend_.GPU_warp);

    cl_int err;
    ValueType scalar = (ValueType) 1;

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = num_blocks * LocalSize;

    err = ocl_kernel<ValueType>(CL_KERNEL_COO_SPMV_FLAT,
                                OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                LocalSize, GlobalSize,
                                tail, interval_size, this->mat_.row, this->mat_.col, this->mat_.val,
                                scalar, cast_in->vec_, cast_out->vec_, temp_rows, temp_vals);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    err = ocl_kernel<ValueType>(CL_KERNEL_COO_SPMV_REDUCE_UPDATE,
                                OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                LocalSize, LocalSize,
                                active_warps, temp_rows, temp_vals, cast_out->vec_);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    free_ocl(&temp_rows);
    free_ocl(&temp_vals);

    err = ocl_kernel<ValueType>(CL_KERNEL_COO_SPMV_SERIAL,
                                OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                1, 1,
                                this->nnz_, this->mat_.row, this->mat_.col, this->mat_.val,
                                scalar, cast_in->vec_, cast_out->vec_, tail);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorMatrixCOO<ValueType>::ApplyAdd(const BaseVector<ValueType> &in, const ValueType scalar,
                                                        BaseVector<ValueType> *out) const {

  // TODO some devices hang up while waiting for COO_SPMV_FLAT event to finish
  // this is a bug we are fixing in some future release
  if (this->nnz_ > 0) {

    assert (&in != NULL);
    assert (out != NULL);
    assert (in.  get_size() >= 0);
    assert (out->get_size() >= 0);
    assert (in.  get_size() == this->ncol_);
    assert (out->get_size() == this->nrow_);
    
    
    const OCLAcceleratorVector<ValueType> *cast_in = dynamic_cast<const OCLAcceleratorVector<ValueType>*> (&in) ; 
    OCLAcceleratorVector<ValueType> *cast_out      = dynamic_cast<      OCLAcceleratorVector<ValueType>*> (out) ; 
    
    assert (cast_in  != NULL);
    assert (cast_out != NULL);

    // do not support super small matrices
    assert (this->nnz_ > this->local_backend_.OCL_warp_size); 

    // ----------------------------------------------------------
    // Modified and adapted from CUSP 0.3.1, 
    // http://code.google.com/p/cusp-library/
    // NVIDIA, APACHE LICENSE 2.0
    // ----------------------------------------------------------
    // see __spmv_coo_flat(...)
    // ----------------------------------------------------------
    // CHANGELOG
    // - adapted interface
    // ----------------------------------------------------------  

    //TODO
    //move in extra file -  max_active_blocks, warp_size, block_size

    const int BLOCK_SIZE = (int) this->local_backend_.OCL_max_work_group_size;
    //    const unsigned int MAX_BLOCKS = this->local_backend_.GPU_max_blocks;

    const unsigned int MAX_BLOCKS = 32; //  cusp::detail::device::arch::max_active_blocks(spmv_coo_flat_kernel<IndexType, ValueType, BLOCK_SIZE, UseCache>, BLOCK_SIZE, (size_t) 0);

    const unsigned int WARPS_PER_BLOCK = BLOCK_SIZE / this->local_backend_.OCL_warp_size;


    const unsigned int num_units  = this->nnz_ / this->local_backend_.OCL_warp_size; 
    const unsigned int num_warps  = std::min(num_units, WARPS_PER_BLOCK * MAX_BLOCKS);
    const unsigned int num_blocks = (num_warps + (WARPS_PER_BLOCK-1)) / WARPS_PER_BLOCK; // (N + (granularity - 1)) / granularity
    const unsigned int num_iters  = (num_units +  (num_warps-1)) / num_warps;
    
    const unsigned int interval_size = this->local_backend_.OCL_warp_size * num_iters;

    const int tail = num_units * this->local_backend_.OCL_warp_size; // do the last few nonzeros separately (fewer than this->local_backend_.GPU_warp elements)

    const unsigned int active_warps = (interval_size == 0) ? 0 : ((tail + (interval_size-1))/interval_size);

    int *temp_rows = NULL;
    ValueType *temp_vals = NULL;

    allocate_ocl(active_warps, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &temp_rows);
    allocate_ocl(active_warps, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &temp_vals);

    //  LOG_INFO("active_warps = " << active_warps);
    //  LOG_INFO("tail =" << tail);
    //  LOG_INFO("interval_size =" << interval_size);
    //  LOG_INFO("num_iters =" << num_iters);
    //  LOG_INFO("num_blocks =" << num_blocks);
    //  LOG_INFO("num_warps =" << num_warps);
    //  LOG_INFO("num_units =" << num_units);
    //  LOG_INFO("WARPS_PER_BLOCK =" << WARPS_PER_BLOCK);
    //  LOG_INFO("MAX_BLOCKS =" << MAX_BLOCKS);
    //  LOG_INFO("BLOCK_SIZE =" << BLOCK_SIZE);
    //  LOG_INFO("WARP_SIZE =" << WARP_SIZE);
    //  LOG_INFO("WARP_SIZE =" << this->local_backend_.GPU_warp);

    cl_int err;

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = num_blocks * LocalSize;

    err = ocl_kernel<ValueType>(CL_KERNEL_COO_SPMV_FLAT,
                                OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                LocalSize, GlobalSize,
                                tail, interval_size, this->mat_.row, this->mat_.col, this->mat_.val,
                                scalar, cast_in->vec_, cast_out->vec_, temp_rows, temp_vals);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    err = ocl_kernel<ValueType>(CL_KERNEL_COO_SPMV_REDUCE_UPDATE,
                                OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                LocalSize, LocalSize,
                                active_warps, temp_rows, temp_vals, cast_out->vec_);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    free_ocl(&temp_rows);
    free_ocl(&temp_vals);

    err = ocl_kernel<ValueType>(CL_KERNEL_COO_SPMV_SERIAL,
                                OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                1, 1,
                                this->nnz_, this->mat_.row, this->mat_.col, this->mat_.val,
                                scalar, cast_in->vec_, cast_out->vec_, tail);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

template <typename ValueType>
bool OCLAcceleratorMatrixCOO<ValueType>::Permute(const BaseVector<int> &permutation) {
/*
  assert (&permutation != NULL);

  // symmetric permutation only
  assert (permutation.get_size() == this->nrow_);
  assert (permutation.get_size() == this->ncol_);

  if (this->nnz_ > 0) {

    const OCLAcceleratorVector<int> *cast_perm = dynamic_cast<const OCLAcceleratorVector<int>*> (&permutation);
    assert (cast_perm != NULL);

    OCLAcceleratorMatrixCOO<ValueType> src(this->local_backend_);
    src.AllocateCOO(this->nnz_, this->nrow_, this->ncol_);
    src.CopyFrom(*this);

    LOG_INFO("OCLAcceleratorMatrixCOO::Permute NYI")
    FATAL_ERROR(__FILE__, __LINE__);

  }
*/
  return false;

}

template <typename ValueType>
bool OCLAcceleratorMatrixCOO<ValueType>::PermuteBackward(const BaseVector<int> &permutation) {
/*
  assert (&permutation != NULL);

  // symmetric permutation only
  assert (permutation.get_size() == this->nrow_);
  assert (permutation.get_size() == this->ncol_);

  if (this->nnz_ > 0) {

    const OCLAcceleratorVector<int> *cast_perm = dynamic_cast<const OCLAcceleratorVector<int>*> (&permutation);
    assert (cast_perm != NULL);

    LOG_INFO("OCLAcceleratorMatrixCOO::PermuteBackward NYI");
    FATAL_ERROR(__FILE__, __LINE__);

  }
*/
  return false;

}


template class OCLAcceleratorMatrixCOO<double>;
template class OCLAcceleratorMatrixCOO<float>;

}
