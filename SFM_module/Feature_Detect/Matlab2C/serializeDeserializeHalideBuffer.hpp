// Copyright 2019-2020 The MathWorks, Inc.

#ifndef MATLAB_HALIDE_CONVERTER
#define MATLAB_HALIDE_CONVERTER

#include <cstring>
#include "rtwtypes.h"

#ifdef HALIDE_CODER
#include "HalideRuntime.h"
#endif

#ifdef _MSC_VER
#pragma warning (disable: 4200)
#pragma warning (disable: 4100)
#endif

template <class T>
void halideBufferToMatlabArray(halide_buffer_t* wrapee, T* out, const int size) {
    out = reinterpret_cast<T*>(wrapee->host);
}

template <class T>
halide_type_t getHalideType(const T* in) {
    halide_type_t htype = halide_type_of<T>();
    return htype;
}

template <class T>
halide_buffer_t matlabArrayToHalideBuffer(const T* wrapee, const int* size, const int dim) {
    halide_buffer_t hbuffer = {};
    hbuffer.dimensions = dim;
    hbuffer.dim = new halide_dimension_t[dim];
    int total_size = (dim > 0) ? 1 : 0;
    for (int i = 0; i < dim; i++) {
        hbuffer.dim[i].min = 0;

        if (i != 0) {
            total_size = total_size * size[i - 1];
        }
       
        hbuffer.dim[i].stride = static_cast<int32_t>(total_size);
        hbuffer.dim[i].extent = static_cast<int32_t>(size[i]);
    }

    hbuffer.host = const_cast<uint8_T*>(reinterpret_cast<const uint8_T*>(wrapee));
    hbuffer.type = getHalideType(wrapee);
    hbuffer.set_host_dirty(true);
    hbuffer.set_device_dirty(false);
    return hbuffer;
}

template <class T>
halide_buffer_t matlabArrayToHalideBuffer(const T* wrapee, const int* min, const int* size, const int* stride, const int dim) {
    halide_buffer_t hbuffer = {};
    hbuffer.dimensions = dim;
    hbuffer.dim = new halide_dimension_t[dim];
    for (int i = 0; i < dim; i++) {
        hbuffer.dim[i].min = static_cast<int32_t>(min[i]);       
        hbuffer.dim[i].extent = static_cast<int32_t>(size[i]);
        hbuffer.dim[i].stride = static_cast<int32_t>(stride[i]);
    }

    hbuffer.host = const_cast<uint8_T*>(reinterpret_cast<const uint8_T*>(wrapee));
    hbuffer.type = getHalideType(wrapee);
    hbuffer.set_host_dirty(true);
    hbuffer.set_device_dirty(false);
    return hbuffer;
}

template <class T>
void deallocateHalideBuffer(T* wrapee) {
    if (wrapee->dim) {
        delete[] (wrapee->dim);
    }
}

#endif
