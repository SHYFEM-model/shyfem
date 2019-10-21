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
#include "../../utils/math_functions.hpp"
#include "kernels_ocl.hpp"
#include "ocl_allocate_free.hpp"
#include "ocl_utils.hpp"
#include "ocl_vector.hpp"
#include "../host/host_vector.hpp"
#include "../backend_manager.hpp"

#include <math.h>

namespace paralution {

template <typename ValueType>
OCLAcceleratorVector<ValueType>::OCLAcceleratorVector() {

  // no default constructors
  LOG_INFO("no default constructor");
  FATAL_ERROR(__FILE__, __LINE__);

}

template <typename ValueType>
OCLAcceleratorVector<ValueType>::OCLAcceleratorVector(const Paralution_Backend_Descriptor local_backend) {

  LOG_DEBUG(this, "OCLAcceleratorVector::OCLAcceleratorVector()",
            "constructor with local_backend");

  this->vec_ = NULL;

  this->set_backend(local_backend);

}

template <typename ValueType>
OCLAcceleratorVector<ValueType>::~OCLAcceleratorVector() {

  LOG_DEBUG(this, "OCLAcceleratorVector::~OCLAcceleratorVector()",
            "destructor");

  this->Clear();

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::info(void) const {

  LOG_INFO("OCLAcceleratorVector<ValueType>");

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::Allocate(const int n) {

  assert (n >= 0);

  if (this->size_ > 0)
    this->Clear();

  if (n > 0) {

    allocate_ocl(n, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &this->vec_);

    ocl_set_to(n, (ValueType) 0, this->vec_, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    this->size_ = n;

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::SetDataPtr(ValueType **ptr, const int size) {

  assert (*ptr != NULL);
  assert (size > 0);

  cl_int err = clFinish(OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
  CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  this->vec_ = *ptr;
  this->size_ = size;

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::LeaveDataPtr(ValueType **ptr) {

  assert (this->size_ > 0);

  cl_int err = clFinish(OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
  CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  *ptr = this->vec_;
  this->vec_ = NULL;

  this->size_ = 0;

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::Clear(void) {

  if (this->size_ > 0) {

    free_ocl(&this->vec_);

    this->size_ = 0;

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::CopyFromHost(const HostVector<ValueType> &src) {

  assert (&src != NULL);

  // CPU to OCL copy
  const HostVector<ValueType> *cast_vec;
  if ((cast_vec = dynamic_cast<const HostVector<ValueType>*> (&src)) != NULL) {

    if (this->size_ == 0)
      this->Allocate(cast_vec->size_);

    assert (cast_vec->size_ == this->size_);

    if (this->size_ > 0) {

      ocl_host2dev(this->size_,    // size
                   cast_vec->vec_, // src
                   this->vec_,     // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

    }

  } else {

    LOG_INFO("Error unsupported OpenCL vector type");
    this->info();
    src.info();
    FATAL_ERROR(__FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::CopyToHost(HostVector<ValueType> *dst) const {

  assert (dst != NULL);

  // OCL to CPU copy
  HostVector<ValueType> *cast_vec;
  if ((cast_vec = dynamic_cast<HostVector<ValueType>*> (dst)) != NULL) {

    if (cast_vec->size_ == 0)
      cast_vec->Allocate(this->size_);

    assert (cast_vec->size_ == this->size_);

    if (this->size_ > 0) {

      ocl_dev2host(this->size_,    // size
                   this->vec_,     // src
                   cast_vec->vec_, // dst
                   OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    }

  } else {

    LOG_INFO("Error unsupported OCL vector type");
    this->info();
    dst->info();
    FATAL_ERROR(__FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::CopyFrom(const BaseVector<ValueType> &src) {

  assert (&src != NULL);

  const OCLAcceleratorVector<ValueType> *ocl_cast_vec;
  const HostVector<ValueType> *host_cast_vec;

  // OCL to OCL copy
  if ((ocl_cast_vec = dynamic_cast<const OCLAcceleratorVector<ValueType>*> (&src)) != NULL) {

    if (this->size_ == 0)
      this->Allocate(ocl_cast_vec->size_);

    assert (ocl_cast_vec->size_ == this->size_);

    if (this != ocl_cast_vec) {

      if (this->size_ > 0) {

        // must be within same opencl context
        ocl_dev2dev(this->size_,        // size
                    ocl_cast_vec->vec_, // src
                    this->vec_,         // dst
                    OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

      }

    }

  } else {

    //OCL to CPU copy
    if ((host_cast_vec = dynamic_cast<const HostVector<ValueType>*> (&src)) != NULL) {

      this->CopyFromHost(*host_cast_vec);

    } else {

      LOG_INFO("Error unsupported OpenCL vector type");
      this->info();
      src.info();
      FATAL_ERROR(__FILE__, __LINE__);

    }

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::CopyFrom(const BaseVector<ValueType> &src,
                                               const int src_offset,
                                               const int dst_offset,
                                               const int size) {

  assert (&src != this);
  assert (this->size_ > 0);
  assert (src.get_size() > 0);
  assert (size > 0);

  assert (src_offset + size <= src.get_size());
  assert (dst_offset + size <= this->size_);

  const OCLAcceleratorVector<ValueType> *cast_src = dynamic_cast<const OCLAcceleratorVector<ValueType>*> (&src);

  assert (cast_src != NULL);

  size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
  size_t GlobalSize = (size / LocalSize + 1) * LocalSize;

  cl_int err = ocl_kernel<ValueType>(CL_KERNEL_COPY_OFFSET_FROM,
                                     OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                     LocalSize, GlobalSize,
                                     size, src_offset, dst_offset, cast_src->vec_, this->vec_);
  CHECK_OCL_ERROR(err, __FILE__, __LINE__);

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::CopyTo(BaseVector<ValueType> *dst) const{

  assert (dst != NULL);

  OCLAcceleratorVector<ValueType> *ocl_cast_vec;
  HostVector<ValueType> *host_cast_vec;

  // OCL to OCL copy
  if ((ocl_cast_vec = dynamic_cast<OCLAcceleratorVector<ValueType>*> (dst)) != NULL) {

    if (this != ocl_cast_vec) {

      if (ocl_cast_vec->size_ == 0)
        ocl_cast_vec->Allocate(this->size_);

      assert (ocl_cast_vec->size_ == this->size_);

      if (this->size_ > 0) {

        // must be within same opencl context
        ocl_dev2dev(this->size_,        // size
                    this->vec_,         // src
                    ocl_cast_vec->vec_, // dst
                    OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
      }

    }

  } else {

    //OCL to CPU copy
    if ((host_cast_vec = dynamic_cast<HostVector<ValueType>*> (dst)) != NULL) {

      this->CopyToHost(host_cast_vec);

    } else {

      LOG_INFO("Error unsupported OpenCL vector type");
      this->info();
      dst->info();
      FATAL_ERROR(__FILE__, __LINE__);

    }

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::CopyFromFloat(const BaseVector<float> &src) {

  assert (&src != NULL);

  const OCLAcceleratorVector<float> *ocl_cast_vec;

  // OCL to OCL copy
  if ((ocl_cast_vec = dynamic_cast<const OCLAcceleratorVector<float>*> (&src)) != NULL) {

    if (this->size_ == 0)
      this->Allocate(ocl_cast_vec->size_);

    assert (ocl_cast_vec->size_ == this->size_);

    if (this->size_ > 0) {

      size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
      size_t GlobalSize = (this->size_ / LocalSize + 1) * LocalSize;

      cl_int err = ocl_kernel<ValueType>(CL_KERNEL_COPY_FROM_FLOAT,
                                         OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                         LocalSize, GlobalSize,
                                         this->size_, ocl_cast_vec->vec_, this->vec_);
      CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    }

  } else {

    LOG_INFO("Error unsupported OCL vector type");
    FATAL_ERROR(__FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::CopyFromDouble(const BaseVector<double> &src) {

  assert (&src != NULL);

  const OCLAcceleratorVector<double> *ocl_cast_vec;

  // GPU to GPU copy
  if ((ocl_cast_vec = dynamic_cast<const OCLAcceleratorVector<double>*> (&src)) != NULL) {

    if (this->size_ == 0)
      this->Allocate(ocl_cast_vec->size_);

    assert (ocl_cast_vec->size_ == this->size_);

    if (this->size_ > 0) {

      size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
      size_t GlobalSize = (this->size_ / LocalSize + 1) * LocalSize;

      cl_int err = ocl_kernel<ValueType>(CL_KERNEL_COPY_FROM_DOUBLE,
                                         OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                         LocalSize, GlobalSize,
                                         this->size_, ocl_cast_vec->vec_, this->vec_);
      CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    }

  } else {
    LOG_INFO("Error unsupported OCL vector type");
    FATAL_ERROR(__FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::Zeros(void) {

  if (this->size_ > 0) {

    ocl_set_to(this->size_, (ValueType) 0, this->vec_, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::Ones(void) {

  if (this->size_ > 0) {

    ocl_set_to(this->size_, (ValueType) 1, this->vec_, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::SetValues(const ValueType val) {

  if (this->size_ > 0) {

    ocl_set_to(this->size_, val, this->vec_, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::AddScale(const BaseVector<ValueType> &x, const ValueType alpha) {

  if (this->size_ > 0) {

    assert (&x != NULL);
    assert (this->size_ == x.get_size());

    const OCLAcceleratorVector<ValueType> *cast_x = dynamic_cast<const OCLAcceleratorVector<ValueType>*> (&x);

    assert (cast_x != NULL);

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = (this->size_ / LocalSize + 1) * LocalSize;

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_AXPY,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, alpha, cast_x->vec_, this->vec_);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::ScaleAdd(const ValueType alpha, const BaseVector<ValueType> &x) {

  if (this->size_ > 0) {

    assert (&x != NULL);
    assert (this->size_ == x.get_size());

    const OCLAcceleratorVector<ValueType> *cast_x = dynamic_cast<const OCLAcceleratorVector<ValueType>*> (&x);

    assert (cast_x != NULL);

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = (this->size_ / LocalSize + 1) * LocalSize;

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_SCALEADD,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, alpha, cast_x->vec_, this->vec_);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::ScaleAddScale(const ValueType alpha,
                                                    const BaseVector<ValueType> &x,
                                                    const ValueType beta) {

  if (this->size_ > 0) {

    assert (&x != NULL);
    assert (this->size_ == x.get_size());

    const OCLAcceleratorVector<ValueType> *cast_x = dynamic_cast<const OCLAcceleratorVector<ValueType>*> (&x);

    assert (cast_x != NULL);

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = (this->size_ / LocalSize + 1) * LocalSize;

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_SCALEADDSCALE,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, alpha, beta, cast_x->vec_, this->vec_);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::ScaleAddScale(const ValueType alpha, const BaseVector<ValueType> &x, const ValueType beta,
                                                    const int src_offset, const int dst_offset,const int size) {

  if (this->size_ > 0) {

    assert (&x != NULL);
    assert (this->size_ > 0);
    assert (x.get_size() > 0);
    assert (size > 0);
    assert (src_offset + size <= x.get_size());
    assert (dst_offset + size <= this->size_);

    const OCLAcceleratorVector<ValueType> *cast_x = dynamic_cast<const OCLAcceleratorVector<ValueType>*> (&x);

    assert (cast_x != NULL);

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = (size / LocalSize + 1) * LocalSize;

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_SCALEADDSCALE_OFFSET,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       size, src_offset, dst_offset, alpha, beta, cast_x->vec_, this->vec_);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::ScaleAdd2(const ValueType alpha, const BaseVector<ValueType> &x,
                                                const ValueType beta,  const BaseVector<ValueType> &y,
                                                const ValueType gamma) {

  if (this->size_ > 0) {

    assert (&x != NULL);
    assert (&y != NULL);
    assert (this->size_ == x.get_size());
    assert (this->size_ == y.get_size());

    const OCLAcceleratorVector<ValueType> *cast_x = dynamic_cast<const OCLAcceleratorVector<ValueType>*> (&x);
    const OCLAcceleratorVector<ValueType> *cast_y = dynamic_cast<const OCLAcceleratorVector<ValueType>*> (&y);

    assert (cast_x != NULL);
    assert (cast_y != NULL);

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = (this->size_ / LocalSize + 1) * LocalSize;

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_SCALEADD2,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, alpha, beta, gamma, cast_x->vec_, cast_y->vec_, this->vec_);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::Scale(const ValueType alpha) {

  if (this->size_ > 0) {

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = (this->size_ / LocalSize + 1) * LocalSize;

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_SCALE,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, alpha, this->vec_);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::ExclusiveScan(const BaseVector<ValueType> &x) {

  LOG_INFO("OCLAcceleratorVector::ExclusiveScan() NYI");
  FATAL_ERROR(__FILE__, __LINE__); 

}

template <typename ValueType>
ValueType OCLAcceleratorVector<ValueType>::Dot(const BaseVector<ValueType> &x) const {

  assert (&x != NULL);
  assert (this->size_ == x.get_size());

  const OCLAcceleratorVector<ValueType> *cast_x = dynamic_cast<const OCLAcceleratorVector<ValueType>*> (&x);

  assert (cast_x != NULL);

  ValueType res = (ValueType) 0;

  if (this->size_ > 0) {

    int finalsize = (int) this->local_backend_.OCL_computeUnits * 4;

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = this->local_backend_.OCL_computeUnits * 4 * LocalSize;

    int GROUP_SIZE = ((this->size_ / finalsize + 1) / (int) LocalSize + 1) * (int) LocalSize;
    int LOCAL_SIZE = GROUP_SIZE / (int) LocalSize;

    ValueType *deviceBuffer = NULL;
    ValueType *hostBuffer = NULL;

    allocate_ocl(finalsize, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &deviceBuffer);

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_DOTC,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, this->vec_, cast_x->vec_, deviceBuffer, GROUP_SIZE, LOCAL_SIZE);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    allocate_host(finalsize, &hostBuffer);
    ocl_dev2host(finalsize, deviceBuffer, hostBuffer, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    free_ocl(&deviceBuffer);

    for (int i=0; i<finalsize; ++i)
      res += hostBuffer[i];

    free_host(&hostBuffer);

  }

  return res;

}

template <typename ValueType>
ValueType OCLAcceleratorVector<ValueType>::DotNonConj(const BaseVector<ValueType> &x) const {

  assert (&x != NULL);
  assert (this->size_ == x.get_size());

  const OCLAcceleratorVector<ValueType> *cast_x = dynamic_cast<const OCLAcceleratorVector<ValueType>*> (&x);

  assert (cast_x != NULL);

  ValueType res = (ValueType) 0;

  if (this->size_ > 0) {

    int finalsize = (int) this->local_backend_.OCL_computeUnits * 4;

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = this->local_backend_.OCL_computeUnits * 4 * LocalSize;

    int GROUP_SIZE = ((this->size_ / finalsize + 1) / (int) LocalSize + 1) * (int) LocalSize;
    int LOCAL_SIZE = GROUP_SIZE / (int) LocalSize;

    ValueType *deviceBuffer = NULL;
    ValueType *hostBuffer = NULL;

    allocate_ocl(finalsize, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &deviceBuffer);

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_DOT,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, this->vec_, cast_x->vec_, deviceBuffer, GROUP_SIZE, LOCAL_SIZE);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    allocate_host(finalsize, &hostBuffer);
    ocl_dev2host(finalsize, deviceBuffer, hostBuffer, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    free_ocl(&deviceBuffer);

    for (int i=0; i<finalsize; ++i)
      res += hostBuffer[i];

    free_host(&hostBuffer);

  }

  return res;

}

template <typename ValueType>
ValueType OCLAcceleratorVector<ValueType>::Norm(void) const {

  ValueType res = (ValueType) 0;

  if (this->size_ > 0) {

    int finalsize = (int) this->local_backend_.OCL_computeUnits * 4;

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = finalsize * LocalSize;

    int GROUP_SIZE = ((this->size_ / finalsize + 1) / (int) LocalSize + 1) * (int) LocalSize;
    int LOCAL_SIZE = GROUP_SIZE / (int) LocalSize;

    ValueType *deviceBuffer = NULL;
    ValueType *hostBuffer = NULL;

    allocate_ocl(finalsize, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &deviceBuffer);

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_NORM,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, this->vec_, deviceBuffer, GROUP_SIZE, LOCAL_SIZE);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    allocate_host(finalsize, &hostBuffer);
    ocl_dev2host(finalsize, deviceBuffer, hostBuffer, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    free_ocl(&deviceBuffer);

    for (int i=0; i<finalsize; ++i)
      res += hostBuffer[i];

    free_host(&hostBuffer);

  }

  return sqrt(res);

}

template <>
int OCLAcceleratorVector<int>::Norm(void) const {

  LOG_INFO("What is int OCLAcceleratorVector<ValueType>::Norm(void) const?");
  FATAL_ERROR(__FILE__, __LINE__);

}

template <typename ValueType>
ValueType OCLAcceleratorVector<ValueType>::Reduce(void) const {

  ValueType res = (ValueType) 0;

  if (this->size_ > 0) {

    int finalsize = (int) this->local_backend_.OCL_computeUnits * 4;

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = finalsize * LocalSize;

    int GROUP_SIZE = ((this->size_ / finalsize + 1) / (int) LocalSize + 1) * (int) LocalSize;
    int LOCAL_SIZE = GROUP_SIZE / (int) LocalSize;

    ValueType *deviceBuffer = NULL;
    ValueType *hostBuffer = NULL;

    allocate_ocl(finalsize, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &deviceBuffer);

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_REDUCE,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, this->vec_, deviceBuffer, GROUP_SIZE, LOCAL_SIZE);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    allocate_host(finalsize, &hostBuffer);
    ocl_dev2host(finalsize, deviceBuffer, hostBuffer, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    free_ocl(&deviceBuffer);

    for (int i=0; i<finalsize; ++i)
      res += hostBuffer[i];

    free_host(&hostBuffer);

  }

  return res;

}

template <typename ValueType>
ValueType OCLAcceleratorVector<ValueType>::Asum(void) const {

  ValueType res = (ValueType) 0;

  if (this->size_ > 0) {

    int finalsize = (int) this->local_backend_.OCL_computeUnits * 4;

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = finalsize * LocalSize;

    int GROUP_SIZE = ((this->size_ / finalsize + 1) / (int) LocalSize + 1) * (int) LocalSize;
    int LOCAL_SIZE = GROUP_SIZE / (int) LocalSize;

    ValueType *deviceBuffer = NULL;
    ValueType *hostBuffer = NULL;

    allocate_ocl(finalsize, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &deviceBuffer);

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_ASUM,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, this->vec_, deviceBuffer, GROUP_SIZE, LOCAL_SIZE);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    allocate_host(finalsize, &hostBuffer);
    ocl_dev2host(finalsize, deviceBuffer, hostBuffer, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    free_ocl(&deviceBuffer);

    for (int i=0; i<finalsize; ++i)
      res += paralution_abs(hostBuffer[i]);

    free_host(&hostBuffer);

  }

  return res;

}

template <typename ValueType>
int OCLAcceleratorVector<ValueType>::Amax(ValueType &value) const {

  ValueType res = (ValueType) 0;
  int idx = 0;

  if (this->size_ > 0) {

    int finalsize = (int) this->local_backend_.OCL_computeUnits * 4;

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = finalsize * LocalSize;

    int GROUP_SIZE = ((this->size_ / finalsize + 1) / (int) LocalSize + 1) * (int) LocalSize;
    int LOCAL_SIZE = GROUP_SIZE / (int) LocalSize;

    ValueType *deviceBuffer = NULL;
    int *iDeviceBuffer = NULL;
    ValueType *hostBuffer = NULL;
    int *iHostBuffer = NULL;

    allocate_ocl(finalsize, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &deviceBuffer);
    allocate_ocl(finalsize, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_context, &iDeviceBuffer);

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_AMAX,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, this->vec_, deviceBuffer, iDeviceBuffer, GROUP_SIZE, LOCAL_SIZE);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    allocate_host(finalsize, &hostBuffer);
    allocate_host(finalsize, &iHostBuffer);

    ocl_dev2host(finalsize, deviceBuffer, hostBuffer, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    free_ocl(&deviceBuffer);
    ocl_dev2host(finalsize, iDeviceBuffer, iHostBuffer, OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue);
    free_ocl(&iDeviceBuffer);

    for (int i=0; i<finalsize; ++i) {
      ValueType tmp = paralution_abs(hostBuffer[i]);
      if (res < tmp) {
        res = tmp;
        idx = iHostBuffer[i];
      }
    }

    free_host(&hostBuffer);
    free_host(&iHostBuffer);

  }

  value = res;

  return idx;

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::PointWiseMult(const BaseVector<ValueType> &x) {

  if (this->size_ > 0) {

    assert (&x != NULL);
    assert (this->size_ == x.get_size());

    const OCLAcceleratorVector<ValueType> *cast_x = dynamic_cast<const OCLAcceleratorVector<ValueType>*> (&x);

    assert (cast_x != NULL);

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = (this->size_ / LocalSize + 1) * LocalSize;

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_POINTWISEMULT,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, cast_x->vec_, this->vec_);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::PointWiseMult(const BaseVector<ValueType> &x, const BaseVector<ValueType> &y) {

  if (this->size_ > 0) {

    assert (&x != NULL);
    assert (&y != NULL);
    assert (this->size_ == x.get_size());
    assert (this->size_ == y.get_size());

    const OCLAcceleratorVector<ValueType> *cast_x = dynamic_cast<const OCLAcceleratorVector<ValueType>*> (&x);
    const OCLAcceleratorVector<ValueType> *cast_y = dynamic_cast<const OCLAcceleratorVector<ValueType>*> (&y);

    assert (cast_x != NULL);
    assert (cast_y != NULL);

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = (this->size_ / LocalSize + 1) * LocalSize;

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_POINTWISEMULT2,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, cast_x->vec_, cast_y->vec_, this->vec_);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::Permute(const BaseVector<int> &permutation) {

  if (this->size_ > 0) {

    assert (&permutation != NULL);
    assert (this->size_ == permutation.get_size());

    const OCLAcceleratorVector<int> *cast_perm = dynamic_cast<const OCLAcceleratorVector<int>*> (&permutation);

    assert (cast_perm != NULL);

    OCLAcceleratorVector<ValueType> vec_tmp(this->local_backend_);
    vec_tmp.Allocate(this->size_);
    vec_tmp.CopyFrom(*this);

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = (this->size_ / LocalSize + 1) * LocalSize;

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_PERMUTE,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, cast_perm->vec_, vec_tmp.vec_, this->vec_);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::PermuteBackward(const BaseVector<int> &permutation) {

  if (this->size_ > 0) {

    assert (&permutation != NULL);
    assert (this->size_ == permutation.get_size());

    const OCLAcceleratorVector<int> *cast_perm = dynamic_cast<const OCLAcceleratorVector<int>*> (&permutation);

    assert (cast_perm != NULL);

    OCLAcceleratorVector<ValueType> vec_tmp(this->local_backend_);
    vec_tmp.Allocate(this->size_);
    vec_tmp.CopyFrom(*this);

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = (this->size_ / LocalSize + 1) * LocalSize;

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_PERMUTE_BACKWARD,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, cast_perm->vec_, vec_tmp.vec_, this->vec_);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::CopyFromPermute(const BaseVector<ValueType> &src,
                                                      const BaseVector<int> &permutation) {

  if (this->size_ > 0) {

    assert (&src != NULL);
    assert (&permutation != NULL);
    assert (this != &src);

    const OCLAcceleratorVector<ValueType> *cast_vec = dynamic_cast<const OCLAcceleratorVector<ValueType>*> (&src);
    const OCLAcceleratorVector<int> *cast_perm      = dynamic_cast<const OCLAcceleratorVector<int>*> (&permutation);

    assert (cast_perm != NULL);
    assert (cast_vec  != NULL);

    assert (cast_vec ->size_ == this->size_);
    assert (cast_perm->size_ == this->size_);

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = (this->size_ / LocalSize + 1) * LocalSize;

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_PERMUTE,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, cast_perm->vec_, cast_vec->vec_, this->vec_);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::CopyFromPermuteBackward(const BaseVector<ValueType> &src,
                                                              const BaseVector<int> &permutation) {

  if (this->size_ > 0) {

    assert (&src != NULL);
    assert (&permutation != NULL);
    assert (this != &src);

    const OCLAcceleratorVector<ValueType> *cast_vec = dynamic_cast<const OCLAcceleratorVector<ValueType>*> (&src);
    const OCLAcceleratorVector<int> *cast_perm      = dynamic_cast<const OCLAcceleratorVector<int>*> (&permutation);

    assert (cast_perm != NULL);
    assert (cast_vec  != NULL);

    assert (cast_vec ->size_ == this->size_);
    assert (cast_perm->size_ == this->size_);

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = (this->size_ / LocalSize + 1) * LocalSize;

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_PERMUTE_BACKWARD,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, cast_perm->vec_, cast_vec->vec_, this->vec_);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

template <typename ValueType>
void OCLAcceleratorVector<ValueType>::Power(const double power) {

  if (this->size_ > 0) {

    size_t LocalSize  = this->local_backend_.OCL_max_work_group_size;
    size_t GlobalSize = (this->size_ / LocalSize + 1) * LocalSize;

    cl_int err = ocl_kernel<ValueType>(CL_KERNEL_POWER,
                                       OCL_HANDLE(this->local_backend_.OCL_handle)->OCL_cmdQueue,
                                       LocalSize, GlobalSize,
                                       this->size_, power, this->vec_);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}


template class OCLAcceleratorVector<double>;
template class OCLAcceleratorVector<float>;
template class OCLAcceleratorVector<int>;

}
