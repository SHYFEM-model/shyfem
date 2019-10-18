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
#include "ocl_allocate_free.hpp"
#include "ocl_utils.hpp"

namespace paralution {

// Allocate memory on device
template <typename DataType>
void allocate_ocl(const int size, cl_context ocl_context, DataType **ptr) {

  LOG_DEBUG(0, "allocate_ocl()",
            size);

  if (size > 0) {

    assert (*ptr == NULL);

    cl_int err;

    // Allocate memory on device
    cl_mem data = clCreateBuffer(ocl_context, CL_MEM_READ_WRITE, sizeof(DataType)*size, NULL, &err);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    *ptr = (DataType*) data;

    assert (*ptr != NULL);

  }

}

// Free memory on device
template <typename DataType>
void free_ocl(DataType **ptr) {

  LOG_DEBUG(0, "free_ocl()",
            "");

  // Free memory on device
  cl_int err = clReleaseMemObject((cl_mem) *ptr);
  CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  *ptr = NULL;

}

// Set object on device to specific value (sync)
template<typename DataType>
void ocl_set_to(const int size, const DataType val, DataType *ptr, cl_command_queue ocl_cmdQueue) {

  LOG_DEBUG(0, "ocl_set_to()",
            "size=" << size << " value=" << val);

  assert (ptr != NULL);

  if (size > 0) {

    cl_event event;
    cl_int err = clEnqueueFillBuffer(ocl_cmdQueue, (cl_mem) ptr, &val, sizeof(DataType), 0,
                                     size*sizeof(DataType), 0, NULL, &event);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    err = clWaitForEvents(1, &event);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    err = clReleaseEvent(event);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

// Copy object from host to device memory (sync)
template <typename DataType>
void ocl_host2dev(const int size, const DataType *src, DataType *dst, cl_command_queue ocl_cmdQueue) {

  LOG_DEBUG(0, "ocl_host2dev()",
            size);

  if (size > 0) {

    assert (src != NULL);
    assert (dst != NULL);

    // Copy object from host to device memory
    cl_int err = clEnqueueWriteBuffer(ocl_cmdQueue, (cl_mem) dst, CL_TRUE, 0, size*sizeof(DataType), src, 0, NULL, NULL);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

// Copy object from device to host memory (sync)
template<typename DataType>
void ocl_dev2host(const int size, const DataType *src, DataType *dst, cl_command_queue ocl_cmdQueue) {

  LOG_DEBUG(0, "ocl_dev2host()",
            size);

  if (size > 0) {

    assert (src != NULL);
    assert (dst != NULL);

    // Copy object from device to host memory
    cl_int err = clEnqueueReadBuffer(ocl_cmdQueue, (cl_mem) src, CL_TRUE, 0, size*sizeof(DataType), dst, 0, NULL, NULL);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}

// Copy object from device to device memory (internal copy, sync)
template<typename DataType>
void ocl_dev2dev(const int size, const DataType *src, DataType *dst, cl_command_queue ocl_cmdQueue) {

  LOG_DEBUG(0, "ocl_dev2dev()",
            size);

  if (size > 0) {

    assert (src != NULL);
    assert (dst != NULL);

    // Copy object from device to device memory (internal copy)
    cl_event event;
    cl_int err = clEnqueueCopyBuffer(ocl_cmdQueue, (cl_mem) src, (cl_mem) dst, 0, 0, size*sizeof(DataType), 0, NULL, &event);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    err = clWaitForEvents(1, &event);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

    err = clReleaseEvent(event);
    CHECK_OCL_ERROR(err, __FILE__, __LINE__);

  }

}


template void allocate_ocl<double      >(const int size, cl_context ocl_context, double       **ptr);
template void allocate_ocl<float       >(const int size, cl_context ocl_context, float        **ptr);
template void allocate_ocl<int         >(const int size, cl_context ocl_context, int          **ptr);
template void allocate_ocl<unsigned int>(const int size, cl_context ocl_context, unsigned int **ptr);
template void allocate_ocl<char        >(const int size, cl_context ocl_context, char         **ptr);

template void free_ocl<double      >(double       **ptr);
template void free_ocl<float       >(float        **ptr);
template void free_ocl<int         >(int          **ptr);
template void free_ocl<unsigned int>(unsigned int **ptr);
template void free_ocl<char        >(char         **ptr);

template void ocl_set_to<double      >(const int size, const double       val, double       *ptr, cl_command_queue ocl_cmdQueue);
template void ocl_set_to<float       >(const int size, const float        val, float        *ptr, cl_command_queue ocl_cmdQueue);
template void ocl_set_to<int         >(const int size, const int          val, int          *ptr, cl_command_queue ocl_cmdQueue);
template void ocl_set_to<unsigned int>(const int size, const unsigned int val, unsigned int *ptr, cl_command_queue ocl_cmdQueue);
template void ocl_set_to<char        >(const int size, const char         val, char         *ptr, cl_command_queue ocl_cmdQueue);

template void ocl_host2dev<double      >(const int size, const double       *src, double       *dst, cl_command_queue ocl_cmdQueue);
template void ocl_host2dev<float       >(const int size, const float        *src, float        *dst, cl_command_queue ocl_cmdQueue);
template void ocl_host2dev<int         >(const int size, const int          *src, int          *dst, cl_command_queue ocl_cmdQueue);
template void ocl_host2dev<unsigned int>(const int size, const unsigned int *src, unsigned int *dst, cl_command_queue ocl_cmdQueue);
template void ocl_host2dev<char        >(const int size, const char         *src, char         *dst, cl_command_queue ocl_cmdQueue);

template void ocl_dev2host<double      >(const int size, const double       *src, double       *dst, cl_command_queue ocl_cmdQueue);
template void ocl_dev2host<float       >(const int size, const float        *src, float        *dst, cl_command_queue ocl_cmdQueue);
template void ocl_dev2host<int         >(const int size, const int          *src, int          *dst, cl_command_queue ocl_cmdQueue);
template void ocl_dev2host<unsigned int>(const int size, const unsigned int *src, unsigned int *dst, cl_command_queue ocl_cmdQueue);
template void ocl_dev2host<char        >(const int size, const char         *src, char         *dst, cl_command_queue ocl_cmdQueue);

template void ocl_dev2dev<double      >(const int size, const double       *src, double       *dst, cl_command_queue ocl_cmdQueue);
template void ocl_dev2dev<float       >(const int size, const float        *src, float        *dst, cl_command_queue ocl_cmdQueue);
template void ocl_dev2dev<int         >(const int size, const int          *src, int          *dst, cl_command_queue ocl_cmdQueue);
template void ocl_dev2dev<unsigned int>(const int size, const unsigned int *src, unsigned int *dst, cl_command_queue ocl_cmdQueue);
template void ocl_dev2dev<char        >(const int size, const char         *src, char         *dst, cl_command_queue ocl_cmdQueue);

}
