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
#include "ocl_matrix_ell.hpp"
#include "ocl_vector.hpp"
#include "../host/host_matrix_ell.hpp"
#include "../backend_manager.hpp"

namespace paralution {

template <typename ValueType>
OCLAcceleratorMatrixELL<ValueType>::OCLAcceleratorMatrixELL() {

  // no default constructors
  LOG_INFO("no default constructor");
  FATAL_ERROR(__FILE__, __LINE__);

}

template <typename ValueType>
OCLAcceleratorMatrixELL<ValueType>::OCLAcceleratorMatrixELL(const Paralution_Backend_Descriptor local_backend) {

  LOG_DEBUG(this, "OCLAcceleratorMatrixELL::OCLAcceleratorMatrixELL()",
            "constructor with local_backend");

  this->mat_.val = NULL;
  this->mat_.col = NULL;

  this->mat_.max_row = 0;

  this->set_backend(local_backend);

}

template <typename ValueType>
OCLAcceleratorMatrixELL<ValueType>::~OCLAcceleratorMatrixELL() {

  LOG_DEBUG(this, "OCLAcceleratorMatrixELL::~OCLAcceleratorMatrixELL()",
            "destructor");

  this->Clear();

}

template <typename ValueType>
void OCLAcceleratorMatrixELL<ValueType>::info(void) const {

  LOG_INFO("OCLAcceleratorMatrixELL<ValueType>");

}

template <typename ValueType>
void OCLAcceleratorMatrixELL<ValueType>::AllocateELL(const int nnz, const int nrow, const int ncol, const int max_row) {

  assert (nnz   >= 0);
  assert (ncol  >= 0);
  assert (nrow  >= 0);
  assert (max_row >= 0);

  if (this->nnz_ > 0)
    this->Clear();

  if (nnz > 0) {

    assert (nnz == max_row * nrow);

    allocate_ocl(nnz, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &this->mat_.col);
    allocate_ocl(nnz, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &this->mat_.val);

    ocl_set_to(nnz, (int) 0, this->mat_.col, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    ocl_set_to(nnz, (ValueType) 0, this->mat_.val, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    this->mat_.max_row = max_row;
    this->nrow_ = nrow;
    this->ncol_ = ncol;
    this->nnz_  = nnz;

  }

}

template <typename ValueType>
void OCLAcceleratorMatrixELL<ValueType>::SetDataPtrELL(int **col, ValueType **val,
                                                       const int nnz, const int nrow, const int ncol,
                                                       const int max_row) {

  assert (*col != NULL);
  assert (*val != NULL);
  assert (nnz  > 0);
  assert (nrow > 0);
  assert (ncol > 0);
  assert (max_row > 0);
  assert (max_row * nrow == nnz);

  this->Clear();

  cl_int err = clFinish(OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
  CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  this->mat_.max_row = max_row;
  this->nrow_ = nrow;
  this->ncol_ = ncol;
  this->nnz_  = nnz;

  this->mat_.col = *col;
  this->mat_.val = *val;

}

template <typename ValueType>
void OCLAcceleratorMatrixELL<ValueType>::LeaveDataPtrELL(int **col, ValueType **val, int &max_row) {

  assert (this->nrow_ > 0);
  assert (this->ncol_ > 0);
  assert (this->nnz_  > 0);
  assert (this->mat_.max_row > 0);
  assert (this->mat_.max_row * this->nrow_ == this->nnz_);

  cl_int err = clFinish(OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
  CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  // see free_host function for details
  *col = this->mat_.col;
  *val = this->mat_.val;

  this->mat_.col = NULL;
  this->mat_.val = NULL;

  max_row = this->mat_.max_row;

  this->mat_.max_row = 0;
  this->nrow_ = 0;
  this->ncol_ = 0;
  this->nnz_  = 0;

}

template <typename ValueType>
void OCLAcceleratorMatrixELL<ValueType>::Clear(void) {

  if (this->nnz_ > 0) {

    free_ocl(&this->mat_.val);
    free_ocl(&this->mat_.col);

    this->nrow_ = 0;
    this->ncol_ = 0;
    this->nnz_  = 0;

  }

}

template <typename ValueType>
void OCLAcceleratorMatrixELL<ValueType>::CopyFromHost(const HostMatrix<ValueType> &src) {

  assert (&src != NULL);

  const HostMatrixELL<ValueType> *cast_mat;

  // copy only in the same format
//  assert (this->get_mat_format() == src.get_mat_format());

  // CPU to OCL copy
  if ((cast_mat = dynamic_cast<const HostMatrixELL<ValueType>*> (&src)) != NULL) {

    if (this->nnz_ == 0)
      this->AllocateELL(cast_mat->nnz_, cast_mat->nrow_, cast_mat->ncol_, cast_mat->mat_.max_row);

    assert (this->nnz_  == cast_mat->nnz_);
    assert (this->nrow_ == cast_mat->nrow_);
    assert (this->ncol_ == cast_mat->ncol_);

    if (this->nnz_ > 0) {

      ocl_host2dev(this->nnz_,         // size
                   cast_mat->mat_.col, // src
                   this->mat_.col,     // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

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
void OCLAcceleratorMatrixELL<ValueType>::CopyToHost(HostMatrix<ValueType> *dst) const {

  assert (dst != NULL);

  HostMatrixELL<ValueType> *cast_mat;

  // copy only in the same format
  assert (this->get_mat_format() == dst->get_mat_format());

  // OCL to CPU copy
  if ((cast_mat = dynamic_cast<HostMatrixELL<ValueType>*> (dst)) != NULL) {

    cast_mat->set_backend(this->local_backend_);

    if (cast_mat->nnz_ == 0)
      cast_mat->AllocateELL(this->nnz_, this->nrow_, this->ncol_, this->mat_.max_row);

    assert (this->nnz_  == cast_mat->nnz_);
    assert (this->nrow_ == cast_mat->nrow_);
    assert (this->ncol_ == cast_mat->ncol_);

    if (this->nnz_ > 0) {

      ocl_dev2host(this->nnz_,         // size
                   this->mat_.col,     // src
                   cast_mat->mat_.col, // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

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
void OCLAcceleratorMatrixELL<ValueType>::CopyFrom(const BaseMatrix<ValueType> &src) {

  assert (&src != NULL);

  const OCLAcceleratorMatrixELL<ValueType> *ocl_cast_mat;
  const HostMatrix<ValueType> *host_cast_mat;

  // copy only in the same format
  // assert (this->get_mat_format() == src.get_mat_format());

  // OCL to OCL copy
  if ((ocl_cast_mat = dynamic_cast<const OCLAcceleratorMatrixELL<ValueType>*> (&src)) != NULL) {

    if (this->nnz_ == 0)
      this->AllocateELL(ocl_cast_mat->nnz_, ocl_cast_mat->nrow_, ocl_cast_mat->ncol_, ocl_cast_mat->mat_.max_row);

    assert (this->nnz_  == ocl_cast_mat->nnz_);
    assert (this->nrow_ == ocl_cast_mat->nrow_);
    assert (this->ncol_ == ocl_cast_mat->ncol_);

    if (this->nnz_ > 0) {

      ocl_dev2dev(this->nnz_,             // size
                  ocl_cast_mat->mat_.col, // src
                  this->mat_.col,         // dst
                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

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
void OCLAcceleratorMatrixELL<ValueType>::CopyTo(BaseMatrix<ValueType> *dst) const {

  assert (dst != NULL);

  OCLAcceleratorMatrixELL<ValueType> *ocl_cast_mat;
  HostMatrix<ValueType> *host_cast_mat;

  // copy only in the same format
  assert (this->get_mat_format() == dst->get_mat_format());

  // OCL to OCL copy
  if ((ocl_cast_mat = dynamic_cast<OCLAcceleratorMatrixELL<ValueType>*> (dst)) != NULL) {

    ocl_cast_mat->set_backend(this->local_backend_);

    if (this->nnz_ == 0)
      ocl_cast_mat->AllocateELL(ocl_cast_mat->nnz_, ocl_cast_mat->nrow_, ocl_cast_mat->ncol_, ocl_cast_mat->mat_.max_row);

    assert (this->nnz_  == ocl_cast_mat->nnz_);
    assert (this->nrow_ == ocl_cast_mat->nrow_);
    assert (this->ncol_ == ocl_cast_mat->ncol_);

    if (this->nnz_ > 0) {

      ocl_dev2dev(this->nnz_,             // size
                  this->mat_.col,         // src
                  ocl_cast_mat->mat_.col, // dst
                  OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

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
bool OCLAcceleratorMatrixELL<ValueType>::ConvertFrom(const BaseMatrix<ValueType> &mat) {

  assert (&mat != NULL);

  this->Clear();

  // empty matrix is empty matrix
  if (mat.get_nnz() == 0)
    return true;

  const OCLAcceleratorMatrixELL<ValueType> *cast_mat_ell;

  if ((cast_mat_ell = dynamic_cast<const OCLAcceleratorMatrixELL<ValueType>*> (&mat)) != NULL) {

    this->CopyFrom(*cast_mat_ell);
    return true;

  }

  const OCLAcceleratorMatrixCSR<ValueType> *cast_mat_csr;
  if ((cast_mat_csr = dynamic_cast<const OCLAcceleratorMatrixCSR<ValueType>*> (&mat)) != NULL) {

    this->Clear();

    assert (cast_mat_csr->nrow_ > 0);
    assert (cast_mat_csr->ncol_ > 0);
    assert (cast_mat_csr->nnz_  > 0);

    int max_row = 0;
    int nrow = cast_mat_csr->nrow_;
    int reducesize = (int) this->local_backend_.OCL_computeUnits * 4;

    int *deviceBuffer = NULL;
    int *hostBuffer   = NULL;

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = reducesize * LocalSize;

    int GROUP_SIZE = ((nrow / reducesize + 1) / (int) LocalSize + 1) * (int) LocalSize;
    int LOCAL_SIZE = GROUP_SIZE / (int) LocalSize;

    allocate_ocl(reducesize, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &deviceBuffer);

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_ELL_MAX_ROW,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       nrow, cast_mat_csr->mat_.row_offset, deviceBuffer, GROUP_SIZE, LOCAL_SIZE);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    allocate_host(reducesize, &hostBuffer);
    ocl_dev2host(reducesize, deviceBuffer, hostBuffer, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    free_ocl(&deviceBuffer);

    for (int i=0; i<reducesize; ++i)
      if (max_row < hostBuffer[i])
        max_row = hostBuffer[i];

    free_host(&hostBuffer);

    int nnz_ell = max_row * nrow;

    this->AllocateELL(nnz_ell, nrow, cast_mat_csr->ncol_, max_row);

    ocl_set_to(nnz_ell, (int) 0, this->mat_.col, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    ocl_set_to(nnz_ell, (ValueType) 0, this->mat_.val, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    GlobalSize = (nrow / LocalSize + 1) * LocalSize;

    err = ocl_kernel<ValueType>(CL_KERNEL_ELL_CSR_TO_ELL,
                                OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                LocalSize, GlobalSize,
                                nrow, max_row,
                                cast_mat_csr->mat_.row_offset, cast_mat_csr->mat_.col, cast_mat_csr->mat_.val,
                                this->mat_.col, this->mat_.val);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    this->mat_.max_row = max_row;
    this->nrow_ = cast_mat_csr->nrow_;
    this->ncol_ = cast_mat_csr->ncol_;
    this->nnz_  = max_row * nrow;

    return true;

  }

  return false;

}

template <typename ValueType>
void OCLAcceleratorMatrixELL<ValueType>::Apply(const BaseVector<ValueType> &in, BaseVector<ValueType> *out) const {

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

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = (this->nrow_ / LocalSize + 1) * LocalSize;

    // Nathan Bell and Michael Garland
    // Efficient Sparse Matrix-Vector Multiplication on {CUDA}
    // NVR-2008-004 / NVIDIA Technical Report
    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_ELL_SPMV,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->nrow_, this->ncol_, this->mat_.max_row, this->mat_.col, this->mat_.val,
                                       cast_in->vec_, cast_out->vec_);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorMatrixELL<ValueType>::ApplyAdd(const BaseVector<ValueType> &in, const ValueType scalar,
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

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = (this->nrow_ / LocalSize + 1) * LocalSize;

    // Nathan Bell and Michael Garland
    // Efficient Sparse Matrix-Vector Multiplication on {CUDA}
    // NVR-2008-004 / NVIDIA Technical Report
    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_ELL_ADD_SPMV,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->nrow_, this->ncol_, this->mat_.max_row, this->mat_.col, this->mat_.val,
                                       scalar, cast_in->vec_, cast_out->vec_);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}


template class OCLAcceleratorMatrixELL<double>;
template class OCLAcceleratorMatrixELL<float>;

}
