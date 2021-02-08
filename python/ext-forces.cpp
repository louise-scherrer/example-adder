///////////////////////////////////////////////////////////////////////////////
// Writting a DAM similar to free-fwddyn but taking into account exterior forces applying on the model
// BSD 3-Clause License
//
// Copyright (C) 2018-2020, LAAS-CNRS, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include "example-adder/python.hpp"
#include "example-adder/ext-forces.hpp"

//#include "python/crocoddyl/multibody/multibody.hpp"
//#include "python/crocoddyl/core/diff-action-base.hpp"

#include <boost/python.hpp>

namespace gepetto {
namespace example {

namespace bp = boost::python;

void exposeDifferentialActionFreeFwdDynamicsExtForces() {
  bp::class_<DifferentialActionModelFreeFwdDynamicsExtForces, bp::bases<DifferentialActionModelAbstract> >(
      "DifferentialActionModelFreeFwdDynamicsExtForces",
      "WORK IN PROGRESS (Feb.2021) \n"
      "Differential action model for free forward dynamics in multibody systems subject to external forces.\n\n"
      "This class implements a the dynamics using Articulate Body Algorithm (ABA),\n"
      "or a custom implementation in case of system with armatures. If you want to\n"
      "include the armature, you need to use setArmature(). On the other hand, the\n"
      "stack of cost functions are implemented in CostModelSum().\n"
      "CAUTION: exterior forces only dealt with in the case of ABA use, not yet for the custom implementation",

      bp::init<boost::shared_ptr<StateMultibody>, boost::shared_ptr<ActuationModelAbstract>,
               boost::shared_ptr<CostModelSum>, 
               boost::shared_ptr<ForceAlignedVector> >(bp::args("self", "state", "actuation", "costs", "extforces"),
                                                 "Initialize the free forward-dynamics action model.\n\n"
                                                 ":param state: multibody state\n"
                                                 ":param actuation: abstract actuation model\n"
                                                 ":param costs: stack of cost functions\n"
                                                 ":param extforces: vector of exterior forces expressed in the local frame of the joints"))
      .def<void (DifferentialActionModelFreeFwdDynamicsExtForces::*)(const boost::shared_ptr<DifferentialActionDataAbstract>&,
                                                            const Eigen::Ref<const Eigen::VectorXd>&,
                                                            const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calc", &DifferentialActionModelFreeFwdDynamicsExtForces::calc, bp::args("self", "data", "x", "u"),
          "Compute the next state and cost value.\n\n"
          "It describes the time-continuous evolution of the multibody system without any contact.\n"
          "Additionally it computes the cost value associated to this state and control pair.\n"
          ":param data: free forward-dynamics action data\n"
          ":param x: time-continuous state vector\n"
          ":param u: time-continuous control input")
      .def<void (DifferentialActionModelFreeFwdDynamicsExtForces::*)(const boost::shared_ptr<DifferentialActionDataAbstract>&,
                                                            const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calc", &DifferentialActionModelAbstract::calc, bp::args("self", "data", "x"))
      .def<void (DifferentialActionModelFreeFwdDynamicsExtForces::*)(const boost::shared_ptr<DifferentialActionDataAbstract>&,
                                                            const Eigen::Ref<const Eigen::VectorXd>&,
                                                            const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calcDiff", &DifferentialActionModelFreeFwdDynamicsExtForces::calcDiff, bp::args("self", "data", "x", "u"),
          "Compute the derivatives of the differential multibody system (free of contact) and\n"
          "its cost functions.\n\n"
          "It computes the partial derivatives of the differential multibody system and the\n"
          "cost function. It assumes that calc has been run first.\n"
          "This function builds a quadratic approximation of the\n"
          "action model (i.e. dynamical system and cost function).\n"
          ":param data: free forward-dynamics action data\n"
          ":param x: time-continuous state vector\n"
          ":param u: time-continuous control input\n")
      .def<void (DifferentialActionModelFreeFwdDynamicsExtForces::*)(const boost::shared_ptr<DifferentialActionDataAbstract>&,
                                                            const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calcDiff", &DifferentialActionModelAbstract::calcDiff, bp::args("self", "data", "x"))
      .def("createData", &DifferentialActionModelFreeFwdDynamicsExtForces::createData, bp::args("self"),
           "Create the free forward dynamics differential action data.")
      .add_property(
          "pinocchio",
          bp::make_function(&DifferentialActionModelFreeFwdDynamicsExtForces::get_pinocchio, bp::return_internal_reference<>()),
          "multibody model (i.e. pinocchio model)")
      .add_property("actuation",
                    bp::make_function(&DifferentialActionModelFreeFwdDynamicsExtForces::get_actuation,
                                      bp::return_value_policy<bp::return_by_value>()),
                    "actuation model")
      .add_property("costs",
                    bp::make_function(&DifferentialActionModelFreeFwdDynamicsExtForces::get_costs,
                                      bp::return_value_policy<bp::return_by_value>()),
                    "total cost model")
      .add_property(
          "armature",
          bp::make_function(&DifferentialActionModelFreeFwdDynamicsExtForces::get_armature, bp::return_internal_reference<>()),
          bp::make_function(&DifferentialActionModelFreeFwdDynamicsExtForces::set_armature),
          "set an armature mechanism in the joints");

  bp::register_ptr_to_python<boost::shared_ptr<DifferentialActionDataFreeFwdDynamicsExtForces> >();

  bp::class_<DifferentialActionDataFreeFwdDynamicsExtForces, bp::bases<DifferentialActionDataAbstract> >(
      "DifferentialActionDataFreeFwdDynamicsExtForces", "Action data for the free forward dynamics system subject to external forces.",
      bp::init<DifferentialActionModelFreeFwdDynamicsExtForces*>(bp::args("self", "model"),
                                                        "Create free forward-dynamics action data.\n\n"
                                                        ":param model: free forward-dynamics action model"))
      .add_property(
          "pinocchio",
          bp::make_getter(&DifferentialActionDataFreeFwdDynamicsExtForces::pinocchio, bp::return_internal_reference<>()),
          "pinocchio data")
      .add_property(
          "multibody",
          bp::make_getter(&DifferentialActionDataFreeFwdDynamicsExtForces::multibody, bp::return_internal_reference<>()),
          "multibody data")
      .add_property("costs",
                    bp::make_getter(&DifferentialActionDataFreeFwdDynamicsExtForces::costs,
                                    bp::return_value_policy<bp::return_by_value>()),
                    "total cost data")
      .add_property("Minv",
                    bp::make_getter(&DifferentialActionDataFreeFwdDynamicsExtForces::Minv, bp::return_internal_reference<>()),
                    "inverse of the joint-space inertia matrix")
      .add_property(
          "u_drift",
          bp::make_getter(&DifferentialActionDataFreeFwdDynamicsExtForces::u_drift, bp::return_internal_reference<>()),
          "force-bias vector that accounts for control, Coriolis and gravitational effects");
}

}  // namespace example
}  // namespace gepetto
