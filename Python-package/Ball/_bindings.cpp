#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

extern "C" {
#include "BD.h"
#include "BI.h"
#include "bcor.h"
#include "kbcov.h"
}

static py::tuple py_bd_test(py::array_t<double, py::array::c_style> xy,
                            py::array_t<int, py::array::c_style> size,
                            int R, int dst, int nthread) {
    auto xy_buf = xy.request();
    auto size_buf = size.request();

    int k = static_cast<int>(size_buf.size);
    int n = 0;
    int *size_ptr = static_cast<int *>(size_buf.ptr);
    for (int i = 0; i < k; i++) {
        n += size_ptr[i];
    }

    int stat_num = (k == 2) ? 3 : 6;
    py::array_t<double> bd_arr(stat_num);
    py::array_t<double> pv_arr(stat_num);

    bd_test(bd_arr.mutable_data(), pv_arr.mutable_data(),
            static_cast<double *>(xy_buf.ptr),
            size_ptr,
            &n, &k, &dst, &R, &nthread);

    return py::make_tuple(bd_arr, pv_arr);
}

static py::tuple py_bcov_test(py::array_t<double, py::array::c_style> x,
                              py::array_t<double, py::array::c_style> y,
                              int n, int R, int dst, int nthread) {
    py::array_t<double> stat_arr(3);
    py::array_t<double> pv_arr(3);

    bcov_test(stat_arr.mutable_data(), pv_arr.mutable_data(),
              x.mutable_data(), y.mutable_data(),
              &n, &R, &dst, &nthread);

    return py::make_tuple(stat_arr, pv_arr);
}

static py::tuple py_kbcov_test(py::array_t<double, py::array::c_style> x,
                               int k, int n, int R, int dst, int nthread) {
    py::array_t<double> stat_arr(3);
    py::array_t<double> pv_arr(3);

    kbcov_test(stat_arr.mutable_data(), pv_arr.mutable_data(),
               x.mutable_data(),
               &k, &n, &R, &dst, &nthread);

    return py::make_tuple(stat_arr, pv_arr);
}

static py::array_t<double> py_bcor_test(py::array_t<double, py::array::c_style> y,
                                        py::array_t<double, py::array::c_style> x,
                                        py::array_t<int, py::array::c_style> x_num,
                                        int f_num, int n, int p, int k,
                                        int dst_y, int dst_x,
                                        int nthread, int complete_flag) {
    py::array_t<double> out(f_num * 3);

    bcor_test(out.mutable_data(),
              y.mutable_data(), x.mutable_data(),
              x_num.mutable_data(),
              &f_num, &n, &p, &k,
              &dst_y, &dst_x, &nthread, &complete_flag);

    return out;
}

PYBIND11_MODULE(_cball, m) {
    m.doc() = "Ball statistics C bindings";
    m.def("py_bd_test", &py_bd_test,
          "Ball Divergence test",
          py::arg("xy"), py::arg("size"),
          py::arg("R"), py::arg("dst"), py::arg("nthread"));
    m.def("py_bcov_test", &py_bcov_test,
          "Ball Covariance test",
          py::arg("x"), py::arg("y"),
          py::arg("n"), py::arg("R"), py::arg("dst"), py::arg("nthread"));
    m.def("py_kbcov_test", &py_kbcov_test,
          "K-sample Ball Covariance test",
          py::arg("x"), py::arg("k"), py::arg("n"),
          py::arg("R"), py::arg("dst"), py::arg("nthread"));
    m.def("py_bcor_test", &py_bcor_test,
          "Ball Correlation test",
          py::arg("y"), py::arg("x"), py::arg("x_num"),
          py::arg("f_num"), py::arg("n"), py::arg("p"), py::arg("k"),
          py::arg("dst_y"), py::arg("dst_x"),
          py::arg("nthread"), py::arg("complete_flag"));
}
