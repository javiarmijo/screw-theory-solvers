#include <MatrixExponential.hpp>
#include <ProductOfExponentials.hpp>
#include <ScrewTheoryIkProblem.hpp>
#include <ConfigurationSelector.hpp>

#include "PySTS.h"

namespace py = pybind11;
using namespace roboticslab;

void init_solvers(pybind11::module &m)
{
    py::class_<MatrixExponential> exp(m, "MatrixExponential");
    py::enum_<MatrixExponential::motion> motion_type(exp, "motion");
    motion_type.value("ROTATION", MatrixExponential::motion::ROTATION);
    motion_type.value("TRANSLATION", MatrixExponential::motion::TRANSLATION);

    exp.def(py::init<MatrixExponential::motion, const KDL::Vector &, const KDL::Vector &>(),
            py::arg("motionType"), py::arg("axis"), py::arg_v("origin", KDL::Vector::Zero(), "Vector.Zero"));
}
