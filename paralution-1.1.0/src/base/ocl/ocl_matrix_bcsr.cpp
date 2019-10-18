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
#include "ocl_matrix_bcsr.hpp"
#include "ocl_vector.hpp"
#include "../host/host_matrix_bcsr.hpp"
#include "../backend_manager.hpp"

namespace paralution {

template <typename ValueType>
OCLAcceleratorMatrixBCSR<ValueType>::OCLAcceleratorMatrixBCSR() {

  // no default constructors
  LOG_INFO("no default constructor");
  FATAL_ERROR(__FILE__, __LINE__);

}

template <typename ValueType>
OCLAcceleratorMatrixBCSR<ValueType>::OCLAcceleratorMatrixBCSR(const Paralution_Backend_Descriptor local_backend) {

  LOG_DEBUG(this, "OCLAcceleratorMatrixBCSR::OCLAcceleratorMatrixBCSR()",
            "constructor with local_backend");

  this->set_backend(local_backend);

  // this is not working anyway
  FATAL_ERROR(__FILE__, __LINE__);

}

template <typename ValueType>
OCLAcceleratorMatrixBCSR<ValueType>::~OCLAcceleratorMatrixBCSR() {

  LOG_DEBUG(this, "OCLAcceleratorMatrixBCSR::~OCLAcceleratorMatrixBCSR()",
            "destructor");

  this->Clear();

}

template <typename ValueType>
void OCLAcceleratorMatrixBCSR<ValueType>::info(void) const {

  LOG_INFO("OCLAcceleratorMatrixBCSR<ValueType>");

}

template <typename ValueType>
void OCLAcceleratorMatrixBCSR<ValueType>::AllocateBCSR(const int nnz, const int nrow, const int ncol) {

  assert (nnz  >= 0);
  assert (ncol >= 0);
  assert (nrow >= 0);

  if (this->nnz_ > 0)
    this->Clear();

  if (nnz > 0) {

    FATAL_ERROR(__FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorMatrixBCSR<ValueType>::Clear(void) {

  if (this->nnz_ > 0) {

    FATAL_ERROR(__FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorMatrixBCSR<ValueType>::CopyFromHost(const HostMatrix<ValueType> &src) {

  assert (&src != NULL);

  const HostMatrixBCSR<ValueType> *cast_mat;

  // copy only in the same format
  assert (this->get_mat_format() == src.get_mat_format());

  // CPU to OCL copy
  if ((cast_mat = dynamic_cast<const HostMatrixBCSR<ValueType>*> (&src)) != NULL) {

    if (this->nnz_ == 0)
      this->AllocateBCSR(cast_mat->nnz_, cast_mat->nrow_, cast_mat->ncol_);

    assert (this->nnz_  == cast_mat->nnz_);
    assert (this->nrow_ == cast_mat->nrow_);
    assert (this->ncol_ == cast_mat->ncol_);

    FATAL_ERROR(__FILE__, __LINE__);

  } else {

    LOG_INFO("Error unsupported OCL matrix type");
    this->info();
    src.info();
    FATAL_ERROR(__FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorMatrixBCSR<ValueType>::CopyToHost(HostMatrix<ValueType> *dst) const {

  assert (dst != NULL);

  HostMatrixBCSR<ValueType> *cast_mat;

  // copy only in the same format
  assert (this->get_mat_format() == dst->get_mat_format());

  // OCL to CPU copy
  if ((cast_mat = dynamic_cast<HostMatrixBCSR<ValueType>*> (dst)) != NULL) {

    cast_mat->set_backend(this->local_backend_);

    if (cast_mat->nnz_ == 0)
      cast_mat->AllocateBCSR(this->nnz_, this->nrow_, this->ncol_);

    assert (this->nnz_  == cast_mat->nnz_);
    assert (this->nrow_ == cast_mat->nrow_);
    assert (this->ncol_ == cast_mat->ncol_);

    FATAL_ERROR(__FILE__, __LINE__);

  } else {

    LOG_INFO("Error unsupported OCL matrix type");
    this->info();
    dst->info();
    FATAL_ERROR(__FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorMatrixBCSR<ValueType>::CopyFrom(const BaseMatrix<ValueType> &src) {

  assert (&src != NULL);

  const OCLAcceleratorMatrixBCSR<ValueType> *ocl_cast_mat;
  const HostMatrix<ValueType> *host_cast_mat;

  // copy only in the same format
  assert (this->get_mat_format() == src.get_mat_format());

  // OCL to OCL copy
  if ((ocl_cast_mat = dynamic_cast<const OCLAcceleratorMatrixBCSR<ValueType>*> (&src)) != NULL) {

    if (this->nnz_ == 0)
      this->AllocateBCSR(ocl_cast_mat->nnz_, ocl_cast_mat->nrow_, ocl_cast_mat->ncol_);

    assert (this->nnz_  == ocl_cast_mat->nnz_);
    assert (this->nrow_ == ocl_cast_mat->nrow_);
    assert (this->ncol_ == ocl_cast_mat->ncol_);

    FATAL_ERROR(__FILE__, __LINE__);

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
void OCLAcceleratorMatrixBCSR<ValueType>::CopyTo(BaseMatrix<ValueType> *dst) const {

  assert (dst != NULL);

  OCLAcceleratorMatrixBCSR<ValueType> *ocl_cast_mat;
  HostMatrix<ValueType> *host_cast_mat;

  // copy only in the same format
  assert (this->get_mat_format() == dst->get_mat_format());

  // OCL to OCL copy
  if ((ocl_cast_mat = dynamic_cast<OCLAcceleratorMatrixBCSR<ValueType>*> (dst)) != NULL) {

    ocl_cast_mat->set_backend(this->local_backend_);

    if (this->nnz_ == 0)
      ocl_cast_mat->AllocateBCSR(this->nnz_, this->nrow_, this->ncol_);

    assert (this->nnz_  == ocl_cast_mat->nnz_);
    assert (this->nrow_ == ocl_cast_mat->nrow_);
    assert (this->ncol_ == ocl_cast_mat->ncol_);

    FATAL_ERROR(__FILE__, __LINE__);

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
bool OCLAcceleratorMatrixBCSR<ValueType>::ConvertFrom(const BaseMatrix<ValueType> &mat) {

  assert (&mat != NULL);

  this->Clear();

  // empty matrix is empty matrix
  if (mat.get_nnz() == 0)
    return true;

  const OCLAcceleratorMatrixBCSR<ValueType> *cast_mat_bcsr;

  if ((cast_mat_bcsr = dynamic_cast<const OCLAcceleratorMatrixBCSR<ValueType>*> (&mat)) != NULL) {

      this->CopyFrom(*cast_mat_bcsr);
      return true;

  }

  /*
    const OCLAcceleratorMatrixCSR<ValueType>  *cast_mat_csr;
    if ((cast_mat_csr = dynamic_cast<const OCLAcceleratorMatrixCSR<ValueType>*> (&mat)) != NULL) {

      this->Clear();

      FATAL_ERROR(__FILE__, __LINE__);

      this->nrow_ = cast_mat_csr->nrow_;
      this->ncol_ = cast_mat_csr->ncol_;
      this->nnz_  = cast_mat_csr->nnz_;

    return 0;

  }
  */

  return false;

}

template <typename ValueType>
void OCLAcceleratorMatrixBCSR<ValueType>::Apply(const BaseVector<ValueType> &in, BaseVector<ValueType> *out) const {
/*
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
*/
  FATAL_ERROR(__FILE__, __LINE__);

  // to avoid compiler warnings
  int err;
  CHECK_OCL_ERROR(err, __FILE__, __LINE__);

}

template <typename ValueType>
void OCLAcceleratorMatrixBCSR<ValueType>::ApplyAdd(const BaseVector<ValueType> &in, const ValueType scalar,
                                                   BaseVector<ValueType> *out) const {
  FATAL_ERROR(__FILE__, __LINE__);
}


template class OCLAcceleratorMatrixBCSR<double>;
template class OCLAcceleratorMatrixBCSR<float>;

}
