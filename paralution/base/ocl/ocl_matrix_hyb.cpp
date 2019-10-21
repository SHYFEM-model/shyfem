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
#include "../../utils/allocate_free.hpp"
#include "kernels_ocl.hpp"
#include "ocl_allocate_free.hpp"
#include "ocl_utils.hpp"
#include "ocl_matrix_csr.hpp"
#include "ocl_matrix_coo.hpp"
#include "ocl_matrix_ell.hpp"
#include "ocl_matrix_hyb.hpp"
#include "ocl_vector.hpp"
#include "../host/host_matrix_hyb.hpp"
#include "../backend_manager.hpp"

#include <algorithm>

namespace paralution {

template <typename ValueType>
OCLAcceleratorMatrixHYB<ValueType>::OCLAcceleratorMatrixHYB() {

  // no default constructors
  LOG_INFO("no default constructor");
  FATAL_ERROR(__FILE__, __LINE__);

}

template <typename ValueType>
OCLAcceleratorMatrixHYB<ValueType>::OCLAcceleratorMatrixHYB(const Paralution_Backend_Descriptor local_backend) {

  LOG_DEBUG(this, "OCLAcceleratorMatrixHYB::OCLAcceleratorMatrixHYB()",
            "constructor with local_backend");

  this->mat_.ELL.val = NULL;
  this->mat_.ELL.col = NULL;
  this->mat_.ELL.max_row = 0;

  this->mat_.COO.row = NULL;
  this->mat_.COO.col = NULL;
  this->mat_.COO.val = NULL;

  this->ell_nnz_ = 0;
  this->coo_nnz_ = 0;

  this->set_backend(local_backend);

}

template <typename ValueType>
OCLAcceleratorMatrixHYB<ValueType>::~OCLAcceleratorMatrixHYB() {

  LOG_DEBUG(this, "OCLAcceleratorMatrixHYB::~OCLAcceleratorMatrixHYB()",
            "destructor");

  this->Clear();

}

template <typename ValueType>
void OCLAcceleratorMatrixHYB<ValueType>::info(void) const {

  LOG_INFO("OCLAcceleratorMatrixHYB<ValueType>");

}

template <typename ValueType>
void OCLAcceleratorMatrixHYB<ValueType>::AllocateHYB(const int ell_nnz, const int coo_nnz, const int ell_max_row,
                                                     const int nrow, const int ncol) {

  assert (ell_nnz >= 0);
  assert (coo_nnz >= 0);
  assert (ell_max_row >= 0);

  assert (ncol >= 0);
  assert (nrow >= 0);

  if (this->nnz_ > 0)
    this->Clear();

  if (ell_nnz + coo_nnz > 0) {

    // ELL
    assert (ell_nnz == ell_max_row * nrow);

    allocate_ocl(ell_nnz, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &this->mat_.ELL.col);
    allocate_ocl(ell_nnz, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &this->mat_.ELL.val);

    ocl_set_to(ell_nnz, (int) 0, this->mat_.ELL.col, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    ocl_set_to(ell_nnz, (ValueType) 0, this->mat_.ELL.val, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    this->mat_.ELL.max_row = ell_max_row;
    this->ell_nnz_ = ell_nnz;

    // COO
    allocate_ocl(coo_nnz, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &this->mat_.COO.row);
    allocate_ocl(coo_nnz, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &this->mat_.COO.col);
    allocate_ocl(coo_nnz, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &this->mat_.COO.val);

    ocl_set_to(coo_nnz, (int) 0, this->mat_.COO.row, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    ocl_set_to(coo_nnz, (int) 0, this->mat_.COO.col, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    ocl_set_to(coo_nnz, (ValueType) 0, this->mat_.COO.val, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    this->nrow_ = nrow;
    this->ncol_ = ncol;

    this->coo_nnz_ = coo_nnz;
    this->nnz_  = ell_nnz + coo_nnz;

  }

}

template <typename ValueType>
void OCLAcceleratorMatrixHYB<ValueType>::Clear(void) {

  if (this->nnz_ > 0) {

    free_ocl(&this->mat_.COO.row);
    free_ocl(&this->mat_.COO.col);
    free_ocl(&this->mat_.COO.val);

    free_ocl(&this->mat_.ELL.val);
    free_ocl(&this->mat_.ELL.col);

    this->ell_nnz_ = 0;
    this->coo_nnz_ = 0;
    this->mat_.ELL.max_row = 0;

    this->nrow_ = 0;
    this->ncol_ = 0;
    this->nnz_  = 0;

  }

}

template <typename ValueType>
void OCLAcceleratorMatrixHYB<ValueType>::CopyFromHost(const HostMatrix<ValueType> &src) {

  assert (&src != NULL);

  const HostMatrixHYB<ValueType> *cast_mat;

  // copy only in the same format
  assert (this->get_mat_format() == src.get_mat_format());

  // CPU to OCL copy
  if ((cast_mat = dynamic_cast<const HostMatrixHYB<ValueType>*> (&src)) != NULL) {

    if (this->nnz_ == 0)
      this->AllocateHYB(cast_mat->ell_nnz_, cast_mat->coo_nnz_, cast_mat->mat_.ELL.max_row,
                        cast_mat->nrow_, cast_mat->ncol_);

    assert (this->nnz_  == cast_mat->nnz_);
    assert (this->nrow_ == cast_mat->nrow_);
    assert (this->ncol_ == cast_mat->ncol_);

    if (this->ell_nnz_ > 0) {

      // ELL
      ocl_host2dev(this->ell_nnz_,         // size
                   cast_mat->mat_.ELL.col, // src
                   this->mat_.ELL.col,     // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      ocl_host2dev(this->ell_nnz_,         // size
                   cast_mat->mat_.ELL.val, // src
                   this->mat_.ELL.val,     // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    }

    if (this->coo_nnz_ > 0) {

      // COO
      ocl_host2dev(this->coo_nnz_,         // size
                   cast_mat->mat_.COO.row, // src
                   this->mat_.COO.row,     // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      ocl_host2dev(this->coo_nnz_,         // size
                   cast_mat->mat_.COO.col, // src
                   this->mat_.COO.col,     // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      ocl_host2dev(this->coo_nnz_,         // size
                   cast_mat->mat_.COO.val, // src
                   this->mat_.COO.val,     // dst
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
void OCLAcceleratorMatrixHYB<ValueType>::CopyToHost(HostMatrix<ValueType> *dst) const {

  assert (dst != NULL);

  HostMatrixHYB<ValueType> *cast_mat;

  // copy only in the same format
  assert (this->get_mat_format() == dst->get_mat_format());

  // OCL to CPU copy
  if ((cast_mat = dynamic_cast<HostMatrixHYB<ValueType>*> (dst)) != NULL) {

    cast_mat->set_backend(this->local_backend_);

    if (cast_mat->nnz_ == 0)
      cast_mat->AllocateHYB(this->ell_nnz_, this->coo_nnz_, this->mat_.ELL.max_row,
                            this->nrow_, this->ncol_);

    assert (this->nnz_  == cast_mat->nnz_);
    assert (this->nrow_ == cast_mat->nrow_);
    assert (this->ncol_ == cast_mat->ncol_);

    if (this->ell_nnz_ > 0) {

      // ELL
      ocl_dev2host(this->ell_nnz_,         // size
                   this->mat_.ELL.col,     // src
                   cast_mat->mat_.ELL.col, // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      ocl_dev2host(this->ell_nnz_,         // size
                   this->mat_.ELL.val,     // src
                   cast_mat->mat_.ELL.val, // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    }

    if (this->coo_nnz_ > 0) {

      // COO
      ocl_dev2host(this->coo_nnz_,         // size
                   this->mat_.COO.row,     // src
                   cast_mat->mat_.COO.row, // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      ocl_dev2host(this->coo_nnz_,         // size
                   this->mat_.COO.col,     // src
                   cast_mat->mat_.COO.col, // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      ocl_dev2host(this->coo_nnz_,         // size
                   this->mat_.COO.val,     // src
                   cast_mat->mat_.COO.val, // dst
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
void OCLAcceleratorMatrixHYB<ValueType>::CopyFrom(const BaseMatrix<ValueType> &src) {

  assert (&src != NULL);

  const OCLAcceleratorMatrixHYB<ValueType> *ocl_cast_mat;
  const HostMatrix<ValueType> *host_cast_mat;

  // copy only in the same format
  assert (this->get_mat_format() == src.get_mat_format());

  // OCL to OCL copy
  if ((ocl_cast_mat = dynamic_cast<const OCLAcceleratorMatrixHYB<ValueType>*> (&src)) != NULL) {

    if (this->nnz_ == 0)
      this->AllocateHYB(ocl_cast_mat->ell_nnz_, ocl_cast_mat->coo_nnz_, ocl_cast_mat->mat_.ELL.max_row,
                        ocl_cast_mat->nrow_, ocl_cast_mat->ncol_);

    assert (this->nnz_  == ocl_cast_mat->nnz_);
    assert (this->nrow_ == ocl_cast_mat->nrow_);
    assert (this->ncol_ == ocl_cast_mat->ncol_);

    if (this->ell_nnz_ > 0) {

      // ELL
      // must be within same opencl context
      ocl_dev2dev(this->ell_nnz_,             // size
                  ocl_cast_mat->mat_.ELL.col, // src
                  this->mat_.ELL.col,         // dst
                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      ocl_dev2dev(this->ell_nnz_,             // size
                  ocl_cast_mat->mat_.ELL.val, // src
                  this->mat_.ELL.val,         // dst
                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    }

    if (this->coo_nnz_ > 0) {

      // COO
      // must be within same opencl context
      ocl_dev2dev(this->coo_nnz_,             // size
                  ocl_cast_mat->mat_.COO.row, // src
                  this->mat_.COO.row,         // dst
                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      ocl_dev2dev(this->coo_nnz_,             // size
                  ocl_cast_mat->mat_.COO.col, // src
                  this->mat_.COO.col,         // dst
                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      ocl_dev2dev(this->coo_nnz_,             // size
                  ocl_cast_mat->mat_.COO.val, // src
                  this->mat_.COO.val,         // dst
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
void OCLAcceleratorMatrixHYB<ValueType>::CopyTo(BaseMatrix<ValueType> *dst) const {

  assert (dst != NULL);

  OCLAcceleratorMatrixHYB<ValueType> *ocl_cast_mat;
  HostMatrix<ValueType> *host_cast_mat;

  // copy only in the same format
  assert (this->get_mat_format() == dst->get_mat_format());

  // OCL to OCL copy
  if ((ocl_cast_mat = dynamic_cast<OCLAcceleratorMatrixHYB<ValueType>*> (dst)) != NULL) {

    ocl_cast_mat->set_backend(this->local_backend_);

    if (ocl_cast_mat->nnz_ == 0)
      ocl_cast_mat->AllocateHYB(this->ell_nnz_, this->coo_nnz_, this->mat_.ELL.max_row,
                                this->nrow_, this->ncol_);

    assert (this->nnz_  == ocl_cast_mat->nnz_);
    assert (this->nrow_ == ocl_cast_mat->nrow_);
    assert (this->ncol_ == ocl_cast_mat->ncol_);

    if (this->ell_nnz_ > 0) {

      // ELL
      // must be within same opencl context
      ocl_dev2dev(this->ell_nnz_,             // size
                  this->mat_.ELL.col,         // src
                  ocl_cast_mat->mat_.ELL.col, // dst
                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      ocl_dev2dev(this->ell_nnz_,             // size
                  this->mat_.ELL.val,         // src
                  ocl_cast_mat->mat_.ELL.val, // dst
                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    }

    if (this->coo_nnz_ > 0) {

      // COO
      // must be within same opencl context
      ocl_dev2dev(this->coo_nnz_,             // size
                  this->mat_.COO.row,         // src
                  ocl_cast_mat->mat_.COO.row, // dst
                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      ocl_dev2dev(this->coo_nnz_,             // size
                  this->mat_.COO.col,         // src
                  ocl_cast_mat->mat_.COO.col, // dst
                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      ocl_dev2dev(this->coo_nnz_,             // size
                  this->mat_.COO.val,         // src
                  ocl_cast_mat->mat_.COO.val, // dst
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
bool OCLAcceleratorMatrixHYB<ValueType>::ConvertFrom(const BaseMatrix<ValueType> &mat) {

  assert (&mat != NULL);

  this->Clear();

  // empty matrix is empty matrix
  if (mat.get_nnz() == 0)
    return true;

  const OCLAcceleratorMatrixHYB<ValueType> *cast_mat_hyb;

  if ((cast_mat_hyb = dynamic_cast<const OCLAcceleratorMatrixHYB<ValueType>*> (&mat)) != NULL) {

      this->CopyFrom(*cast_mat_hyb);
      return true;

  }

  const OCLAcceleratorMatrixCSR<ValueType> *cast_mat_csr;

  if ((cast_mat_csr = dynamic_cast<const OCLAcceleratorMatrixCSR<ValueType>*> (&mat)) != NULL) {

    this->Clear();

    assert (cast_mat_csr->nrow_ > 0);
    assert (cast_mat_csr->ncol_ > 0);
    assert (cast_mat_csr->nnz_ > 0);

    int nrow = cast_mat_csr->nrow_;
    int ncol = cast_mat_csr->ncol_;
    int max_row = cast_mat_csr->nnz_ / nrow;

    // get nnz per row for COO part
    int *nnz_coo = NULL;
    allocate_ocl(nrow, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &nnz_coo);

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = (nrow / LocalSize + 1) * LocalSize;

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_ELL_NNZ_COO,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       nrow, max_row, cast_mat_csr->mat_.row_offset, nnz_coo);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    // TODO fix int kernel for reduce
    int *hostBuffer = NULL;
    allocate_host(nrow, &hostBuffer);
    ocl_dev2host(nrow, nnz_coo, hostBuffer, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    int num_nnz_coo = 0;

    for (int i=0; i<nrow; ++i)
      num_nnz_coo += hostBuffer[i];
    free_host(&hostBuffer);
/*
    // get nnz for COO part by summing up nnz per row array
    int reducesize = (int) this->local_backend_.OCL_computeUnits * 4;
    int *deviceBuffer = NULL;
    int *hostBuffer = NULL;

    allocate_ocl(reducesize, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &deviceBuffer);

    GlobalSize = reducesize * LocalSize;

    int GROUP_SIZE = ((nrow / reducesize + 1) / (int) LocalSize + 1) * (int) LocalSize;
    int LOCAL_SIZE = GROUP_SIZE / (int) LocalSize;

    err = ocl_kernel<int>(CL_KERNEL_REDUCE,
                          OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                          LocalSize, GlobalSize,
                          nrow, nnz_coo, deviceBuffer, GROUP_SIZE, LOCAL_SIZE);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    allocate_host(reducesize, &hostBuffer);

    ocl_dev2host(reducesize, deviceBuffer, hostBuffer, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    free_ocl(&deviceBuffer);

    int num_nnz_coo = 0;
    for (int i=0; i<reducesize; ++i)
      num_nnz_coo += hostBuffer[i];

    free_host(&hostBuffer);
*/
    // END TODO

    // allocate ELL and COO matrices
    int num_nnz_ell = max_row * nrow;

    if (num_nnz_ell <= 0 || num_nnz_coo <= 0) {
      free_ocl(&nnz_coo);
      return false;
    }

    this->AllocateHYB(num_nnz_ell, num_nnz_coo, max_row, nrow, ncol);

    ocl_set_to(num_nnz_ell, (int) -1, this->mat_.ELL.col, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    // copy up to num_cols_per_row values of row i into the ELL
    int *nnz_ell = NULL;

    allocate_ocl(nrow, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &nnz_ell);

    GlobalSize = (nrow / LocalSize + 1) * LocalSize;

    err = ocl_kernel<ValueType>(CL_KERNEL_ELL_FILL_ELL,
                                OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                LocalSize, GlobalSize,
                                nrow, max_row,
                                cast_mat_csr->mat_.row_offset, cast_mat_csr->mat_.col, cast_mat_csr->mat_.val,
                                this->mat_.ELL.col, this->mat_.ELL.val, nnz_ell);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    // TODO currently performing partial sum on host
    allocate_host(nrow, &hostBuffer);
    ocl_dev2host(nrow, nnz_ell, hostBuffer, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    for (int i=1; i<nrow; ++i)
      hostBuffer[i] += hostBuffer[i-1];

    ocl_host2dev(nrow, hostBuffer, nnz_ell, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    free_host(&hostBuffer);
    // end TODO

    // copy any remaining values in row i into the COO
    err = ocl_kernel<ValueType>(CL_KERNEL_ELL_FILL_COO,
                                OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                LocalSize, GlobalSize,
                                nrow, cast_mat_csr->mat_.row_offset, cast_mat_csr->mat_.col, cast_mat_csr->mat_.val,
                                nnz_coo, nnz_ell, this->mat_.COO.row, this->mat_.COO.col, this->mat_.COO.val);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    free_ocl(&nnz_ell);
    free_ocl(&nnz_coo);

    this->nrow_ = cast_mat_csr->nrow_;
    this->ncol_ = cast_mat_csr->ncol_;
    this->nnz_  = num_nnz_ell + num_nnz_coo;
    this->mat_.ELL.max_row = max_row;
    this->ell_nnz_ = num_nnz_ell;
    this->coo_nnz_ = num_nnz_coo;

    return true;

  }

  return false;

}

template <typename ValueType>
void OCLAcceleratorMatrixHYB<ValueType>::Apply(const BaseVector<ValueType> &in, BaseVector<ValueType> *out) const {

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

    // ELL
    if (this->ell_nnz_ > 0) {

      size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
      size_t GlobalSize = (this->nrow_ / LocalSize + 1) * LocalSize;

      cl_int err = ocl_kernel<ValueType>(CL_KERNEL_ELL_SPMV,
                                         OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                         LocalSize, GlobalSize,
                                         this->nrow_, this->ncol_, this->mat_.ELL.max_row, this->mat_.ELL.col, this->mat_.ELL.val,
                                         cast_in->vec_, cast_out->vec_);
      CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    }

    // COO
    if (this->coo_nnz_ > 0) {

      // do not support super small matrices
      assert (this->coo_nnz_ > this->local_backend_.OCL_warp_size);

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

      const int BLOCK_SIZE = int(this->local_backend_.OCL_max_work_group_size);
      //    const unsigned int MAX_BLOCKS = this->local_backend_.GPU_max_blocks;

      const unsigned int MAX_BLOCKS = 32; //  cusp::detail::device::arch::max_active_blocks(spmv_coo_flat_kernel<IndexType, ValueType, BLOCK_SIZE, UseCache>, BLOCK_SIZE, (size_t) 0);

      const unsigned int WARPS_PER_BLOCK = BLOCK_SIZE / this->local_backend_.OCL_warp_size;

      const unsigned int num_units  = this->coo_nnz_ / this->local_backend_.OCL_warp_size;
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

      cl_int err;

      size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
      size_t GlobalSize = num_blocks * LocalSize;

      ValueType scalar = (ValueType) 1;

      err = ocl_kernel<ValueType>(CL_KERNEL_COO_SPMV_FLAT,
                                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                  LocalSize, GlobalSize,
                                  tail, interval_size, this->mat_.COO.row, this->mat_.COO.col, this->mat_.COO.val,
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
                                  this->coo_nnz_, this->mat_.COO.row, this->mat_.COO.col, this->mat_.COO.val,
                                  scalar, cast_in->vec_, cast_out->vec_, tail);
      CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    }

  }

}

template <typename ValueType>
void OCLAcceleratorMatrixHYB<ValueType>::ApplyAdd(const BaseVector<ValueType> &in, const ValueType scalar,
                                                  BaseVector<ValueType> *out) const {

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

    // ELL
    if (this->ell_nnz_ > 0) {

      size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
      size_t GlobalSize = (this->nrow_ / LocalSize + 1) * LocalSize;

      cl_int err = ocl_kernel<ValueType>(CL_KERNEL_ELL_ADD_SPMV,
                                         OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                         LocalSize, GlobalSize,
                                         this->nrow_, this->ncol_, this->mat_.ELL.max_row, this->mat_.ELL.col, this->mat_.ELL.val,
                                         scalar, cast_in->vec_, cast_out->vec_);
      CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    }

    if (this->coo_nnz_ > 0) {

      // do not support super small matrices
      assert (this->coo_nnz_ > this->local_backend_.OCL_warp_size);

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

      const int BLOCK_SIZE = int(this->local_backend_.OCL_max_work_group_size);
      //    const unsigned int MAX_BLOCKS = this->local_backend_.GPU_max_blocks;

      const unsigned int MAX_BLOCKS = 32; //  cusp::detail::device::arch::max_active_blocks(spmv_coo_flat_kernel<IndexType, ValueType, BLOCK_SIZE, UseCache>, BLOCK_SIZE, (size_t) 0);

      const unsigned int WARPS_PER_BLOCK = BLOCK_SIZE / this->local_backend_.OCL_warp_size;

      const unsigned int num_units  = this->coo_nnz_ / this->local_backend_.OCL_warp_size;
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

      cl_int err;

      size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
      size_t GlobalSize = num_blocks * LocalSize;

      err = ocl_kernel<ValueType>(CL_KERNEL_COO_SPMV_FLAT,
                                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                  LocalSize, GlobalSize,
                                  tail, interval_size, this->mat_.COO.row, this->mat_.COO.col, this->mat_.COO.val,
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
                                  this->coo_nnz_, this->mat_.COO.row, this->mat_.COO.col, this->mat_.COO.val,
                                  scalar, cast_in->vec_, cast_out->vec_, tail);
      CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    }

  }

}


template class OCLAcceleratorMatrixHYB<double>;
template class OCLAcceleratorMatrixHYB<float>;

}
