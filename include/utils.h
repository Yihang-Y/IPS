#pragma once
#include <cstddef>
#include <array>
#include <vector>
#include <functional>
#include <memory>
#include <map>
#include <variant>


#ifdef USE_PYBIND11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/embed.h>
namespace py = pybind11;
#endif


#ifdef USE_PYBIND11
PYBIND11_MAKE_OPAQUE(std::vector<std::array<double, 1>>);
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
    ), py::cast(ptr)); // cast to keep the shared_ptr alive, set base
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

using ConfigValue = std::variant<int, double, std::string, bool>;
using ConfigMap = std::map<std::string, ConfigValue>;

struct Config {
    ConfigMap config_map;
    Config() = default;

    ConfigValue& operator[](const std::string& key) {
        return config_map[key];
    }

    const ConfigValue& operator[](const std::string& key) const {
        return config_map.at(key);
    }

    template<typename T>
    T get(const std::string& key, const T& default_value = T{}) const {
        auto it = config_map.find(key);
        if (it == config_map.end()) {
            return default_value;
        }
        return std::get<T>(it->second);
    }
};

using PairForceFactory = std::function<std::function<double(double)>(const Config&)>;
template<size_t dim>
using ConfinementForceFactory = std::function<std::function<vec<dim>(const vec<dim>&)>(const Config&)>;

std::map<std::string, PairForceFactory>& get_PFF_map() {
    static std::map<std::string, PairForceFactory> factoryMap;
    return factoryMap;
}

template<typename T>
std::function<double(double)> create_pair_force(const Config& config) {
    T instance;
    instance.from_config(config);
    return instance;
}

namespace reflect {
    #define REFLECT__PP_FOREACH_1(f, _1) f(_1)
    #define REFLECT__PP_FOREACH_2(f, _1, _2) f(_1) f(_2)
    #define REFLECT__PP_FOREACH_3(f, _1, _2, _3) f(_1) f(_2) f(_3)
    #define REFLECT__PP_FOREACH_4(f, _1, _2, _3, _4) f(_1) f(_2) f(_3) f(_4)
    #define REFLECT__PP_FOREACH_5(f, _1, _2, _3, _4, _5) f(_1) f(_2) f(_3) f(_4) f(_5)
    #define REFLECT__PP_FOREACH_6(f, _1, _2, _3, _4, _5, _6) f(_1) f(_2) f(_3) f(_4) f(_5) f(_6)
    #define REFLECT__PP_FOREACH_7(f, _1, _2, _3, _4, _5, _6, _7) f(_1) f(_2) f(_3) f(_4) f(_5) f(_6) f(_7)
    #define REFLECT__PP_FOREACH_8(f, _1, _2, _3, _4, _5, _6, _7, _8) f(_1) f(_2) f(_3) f(_4) f(_5) f(_6) f(_7) f(_8)

    #define REFLECT__PP_NARGS_IMPL(_1, _2, _3, _4, _5, _6, _7, _8, N, ...) N
    #define REFLECT__PP_NARGS(...) REFLECT__PP_NARGS_IMPL(__VA_ARGS__, 8, 7, 6, 5, 4, 3, 2, 1)

    #define REFLECT__PP_EXPAND_2(...) __VA_ARGS__
    #define REFLECT__PP_EXPAND(...) REFLECT__PP_EXPAND_2(__VA_ARGS__)
    #define REFLECT__PP_CONCAT_2(x, y) x##y
    #define REFLECT__PP_CONCAT(x, y) REFLECT__PP_CONCAT_2(x, y)
    #define REFLECT__PP_FOREACH(f, ...) REFLECT__PP_EXPAND(REFLECT__PP_CONCAT(REFLECT__PP_FOREACH_, REFLECT__PP_NARGS(__VA_ARGS__))(f, __VA_ARGS__))

    #define REFLECT_EACH_MEMBER(member) \
        member = config.get<decltype(this->member)>(#member, decltype(this->member){});


    #define REGISTER_PAIR_FORCE(T) \
    namespace { \
        const bool T##_registered __attribute__((used)) = []() { \
            std::cout << "Registering " << #T << std::endl; \
            get_PFF_map().emplace(#T, create_pair_force<T>); \
            return true; \
        }(); \
    }

    // #define REGISTER_CONFINEMENT(T) \

    #define REFLECT(...) \
        void from_config(const Config& config) { \
            REFLECT__PP_FOREACH(REFLECT_EACH_MEMBER, __VA_ARGS__) \
        }
}