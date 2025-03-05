#pragma once
#include <cstddef>
#include <array>
#include <vector>
#include <functional>
#include <memory>

#ifdef USE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/embed.h>
namespace py = pybind11;
#endif


#ifdef USE_PYBIND11
template<int Dim>
py::array convert_from_shared_ptr(std::shared_ptr<std::vector<std::array<double, Dim>>> ptr)
{
    double* data_ptr = &((*ptr)[0][0]);
    std::vector<ssize_t> shape = { static_cast<ssize_t>(ptr->size()), static_cast<ssize_t>(Dim) };
    std::vector<ssize_t> strides = { static_cast<ssize_t>(Dim * sizeof(double)), static_cast<ssize_t>(sizeof(double)) };
    return py::array_t<double>(py::buffer_info(
        data_ptr,                            // data pointer
        sizeof(double),                      // size of each element
        py::format_descriptor<double>::format(), // data type
        2,                                   // number of dimensions
        shape,                               // shape, i.e. number of elements in each dimension
        strides                              // strides, i.e. bytes to jump to get to the next element
    ));
}

py::array_t<double> vector_to_array(std::vector<double> &vec) {
    return py::array_t<double>(
        {vec.size()},                // shape
        {sizeof(double)},            // stride
        vec.data(),                  // point to the actual data
        // FIXME: ignore the lifetime of vec
        py::cast(nullptr)            // lifetime management policy
    );
}
#endif