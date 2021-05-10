///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2019-2020, LAAS-CNRS, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

//#include "python/crocoddyl/core/core.hpp"
//#include "crocoddyl/core/costs/control.hpp"
#include <boost/python.hpp>

#include "example-adder/python.hpp"
#include "example-adder/cost_control_p_norm.hpp"

//#include "python/crocoddyl/utils/deprecate.hpp"
//#include "crocoddyl/bindings/python/crocoddyl/utils/deprecate.hpp"

namespace gepetto {
namespace example {

namespace bp = boost::python;

void exposeCostControlPNorm() {
  bp::class_<CostModelControlPNorm, bp::bases<CostModelAbstract> >(
      "CostModelControlPNorm",
      "This cost function defines a residual vector as the p-norm of the control vector ie r = (sum_over_i(u_i^p))^(1/p), "
      "with u the current control vector.",
      bp::init<boost::shared_ptr<StateAbstract>, const std::size_t&, const std::size_t&>(
          bp::args("self", "state", "nu", "p"),
          "Initialize the cost model on the p-norm of the control vector.\n\n"
          ":param state: state description\n"
          ":param nu: dimension of control vector\n"
          ":param p: integer value of p for the p-norm"))
      .def<void (CostModelControlPNorm::*)(const boost::shared_ptr<CostDataAbstract>&,
                                      const Eigen::Ref<const Eigen::VectorXd>&,
                                      const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calc", &CostModelControlPNorm::calc, bp::args("self", "data", "x", "u"),
          "Compute the control cost.\n\n"
          ":param data: cost data\n"
          ":param x: time-discrete state vector\n"
          ":param u: time-discrete control input")
      .def<void (CostModelControlPNorm::*)(const boost::shared_ptr<CostDataAbstract>&,
                                      const Eigen::Ref<const Eigen::VectorXd>&)>("calc", &CostModelAbstract::calc,
                                                                                 bp::args("self", "data", "x"))
      .def<void (CostModelControlPNorm::*)(const boost::shared_ptr<CostDataAbstract>&,
                                      const Eigen::Ref<const Eigen::VectorXd>&,
                                      const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calcDiff", &CostModelControlPNorm::calcDiff, bp::args("self", "data", "x", "u"),
          "Compute the derivatives of the control cost.\n\n"
          "It assumes that calc has been run first.\n"
          ":param data: action data\n"
          ":param x: time-discrete state vector\n"
          ":param u: time-discrete control input\n")
      .def<void (CostModelControlPNorm::*)(const boost::shared_ptr<CostDataAbstract>&,
                                      const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calcDiff", &CostModelAbstract::calcDiff, bp::args("self", "data", "x"))
      .def("createData", &CostModelControlPNorm::createData, bp::with_custodian_and_ward_postcall<0, 2>(),
           bp::args("self", "data"),
           "Create the control cost data.\n\n"
           "Each cost model has its own data that needs to be allocated. This function\n"
           "returns the allocated data for a predefined cost.\n"
           ":param data: shared data\n"
           ":return cost data.");
}

}  // namespace example
}  // namespace gepetto

