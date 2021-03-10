///////////////////////////////////////////////////////////////////////////////
// Writting a DAM similar to free-fwddyn but taking into account exterior forces applying on the model
//
//
// Copyright (C) 2018-2020, LAAS-CNRS, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include "crocoddyl/core/utils/exception.hpp"

#include "example-adder/ext-forces.hpp"

#include <pinocchio/algorithm/aba.hpp>
#include <pinocchio/algorithm/aba-derivatives.hpp>
#include <pinocchio/algorithm/rnea.hpp>
#include <pinocchio/algorithm/rnea-derivatives.hpp>
#include <pinocchio/algorithm/compute-all-terms.hpp>
#include <pinocchio/algorithm/kinematics.hpp>
#include <pinocchio/algorithm/jacobian.hpp>
#include <pinocchio/algorithm/frames.hpp>
#include <pinocchio/algorithm/cholesky.hpp>

namespace gepetto {
namespace example {
template <typename Scalar>
DifferentialActionModelFreeFwdDynamicsExtForcesTpl<Scalar>::DifferentialActionModelFreeFwdDynamicsExtForcesTpl(
    boost::shared_ptr<StateMultibody> state, boost::shared_ptr<ActuationModelAbstract> actuation,
    boost::shared_ptr<CostModelSum> costs, const ForceAlignedVector& extforces)
    : Base(state, actuation->get_nu(), costs->get_nr()),
      actuation_(actuation),
      costs_(costs),
      pinocchio_(*state->get_pinocchio().get()),
      with_armature_(true),
      armature_(VectorXs::Zero((long)state->get_nv())),
      extforces_(extforces) {
  if (costs_->get_nu() != nu_) {
    throw_pretty("Invalid argument: "
                 << "Costs doesn't have the same control dimension (it should be " + std::to_string(nu_) + ")");
  }
  if (static_cast<int>(extforces_.size()) != pinocchio_.njoints) {
    throw_pretty("Invalid argument: "
                 << "extforces is of wrong dimension: it should be " + std::to_string(pinocchio_.njoints) + ")");
  }
  Base::set_u_lb(Scalar(-1.) * pinocchio_.effortLimit.tail((long)nu_));
  Base::set_u_ub(Scalar(+1.) * pinocchio_.effortLimit.tail((long)nu_));
}

template <typename Scalar>
DifferentialActionModelFreeFwdDynamicsExtForcesTpl<Scalar>::~DifferentialActionModelFreeFwdDynamicsExtForcesTpl() {}

template <typename Scalar>
void DifferentialActionModelFreeFwdDynamicsExtForcesTpl<Scalar>::calc(
    const boost::shared_ptr<DifferentialActionDataAbstract>& data, const Eigen::Ref<const VectorXs>& x,
    const Eigen::Ref<const VectorXs>& u) {
  if (static_cast<std::size_t>(x.size()) != state_->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "x has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
  }
  if (static_cast<std::size_t>(u.size()) != nu_) {
    throw_pretty("Invalid argument: "
                 << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
  }

  Data* d = static_cast<Data*>(data.get());
  const Eigen::VectorBlock<const Eigen::Ref<const VectorXs>, Eigen::Dynamic> q = x.head((long)state_->get_nq());
  const Eigen::VectorBlock<const Eigen::Ref<const VectorXs>, Eigen::Dynamic> v = x.tail((long)state_->get_nv());

  // Expressing extforces_ (expressed in world frame at joint center) in the local frame of the joints
  ForceAlignedVector localFrameExtForces = extforces_;
  VectorXs testForce(6);
  testForce(0) = 0.;
  testForce(1) = 400.;
  testForce(2) = 0.;
  testForce(3) = 0.;
  testForce(4) = 0.;
  testForce(5) = 0.;
  for (pinocchio::JointIndex i = 1; i < (pinocchio::JointIndex)pinocchio_.njoints; ++i) {
  //TODO Should it start at i = 0 ? nope, cf state.hxx (crocoddyl file)
    localFrameExtForces[i] = d->pinocchio.oMi[i].actInv(extforces_[i]);
    
    //localFrameExtForces[i] = pinocchio::ForceTpl<Scalar>(testForce);
    //if (i == 1) {localFrameExtForces[i] = d->pinocchio.oMi[i].actInv(extforces_[i]);}
    std::cout << "Inside DAM extforces_[" << i << "] = " << extforces_[i] << "vs local = " << localFrameExtForces[i] << std::endl;
    std::cout << "oMi [" << i << "] = " << d->pinocchio.oMi[i] << std::endl;
  }
  //extforces_ = localFrameExtForces; // NOOOOO

  actuation_->calc(d->multibody.actuation, x, u);

  // Computing the dynamics using ABA or manually for armature case
  if (with_armature_) {
    d->xout = pinocchio::aba(pinocchio_, d->pinocchio, q, v, d->multibody.actuation->tau, localFrameExtForces); // previously extforces_
    //std::cout << "Inside DAM data.nle is" << d->pinocchio.nle << std::endl;
    pinocchio::updateGlobalPlacements(pinocchio_, d->pinocchio);
    /*std::cout << "Inside DAM data.nle is" << d->pinocchio.nle << std::endl; // added for quick test on b(q,dq) vector
    //std::cout << "Inside DAM data.f [0] = " << d->pinocchio.f[0] << std::endl;
    //std::cout << "Inside DAM data.f [1] = " << d->pinocchio.f[1] << std::endl;
    //std::cout << "Inside DAM data.f [2] = " << d->pinocchio.f[2] << std::endl;
    //std::cout << "Inside DAM data.f [3] = " << d->pinocchio.f[3] << std::endl;
    std::cout << "Inside DAM data.f [4] = " << d->pinocchio.f[4] << std::endl;
    std::cout << "Inside DAM data.f [5] = " << d->pinocchio.f[5] << std::endl;
    std::cout << "Inside DAM data.f [6] = " << d->pinocchio.f[6] << std::endl;
    std::cout << "Inside DAM data.f [7] = " << d->pinocchio.f[7] << std::endl;
    std::cout << "Inside DAM data.f [8] = " << d->pinocchio.f[8] << std::endl;
    std::cout << "Inside DAM data.f [9] = " << d->pinocchio.f[9] << std::endl;
    std::cout << "Inside DAM data.f [10] = " << d->pinocchio.f[10] << std::endl;
    std::cout << "Inside DAM data.f [11] = " << d->pinocchio.f[11] << std::endl;
    std::cout << "Inside DAM data.f [12] = " << d->pinocchio.f[12] << std::endl;
    std::cout << "Inside DAM data.f [13] = " << d->pinocchio.f[13] << std::endl;
    std::cout << "Inside DAM data.f [14] = " << d->pinocchio.f[14] << std::endl;
    std::cout << "Inside DAM data.f [15] = " << d->pinocchio.f[15] << std::endl;
    std::cout << "Inside DAM data.f [16] = " << d->pinocchio.f[16] << std::endl;
    std::cout << "Inside DAM data.f [17] = " << d->pinocchio.f[17] << std::endl;
    std::cout << "Inside DAM data.f [18] = " << d->pinocchio.f[18] << std::endl;
    std::cout << "Inside DAM data.f [19] = " << d->pinocchio.f[19] << std::endl; */

  } else {
    pinocchio::computeAllTerms(pinocchio_, d->pinocchio, q, v); //TODO: adapt this to take extforces_ into account
    d->pinocchio.M.diagonal() += armature_;
    pinocchio::cholesky::decompose(pinocchio_, d->pinocchio);
    d->Minv.setZero();
    pinocchio::cholesky::computeMinv(pinocchio_, d->pinocchio, d->Minv);
    d->u_drift = d->multibody.actuation->tau - d->pinocchio.nle;
    d->xout.noalias() = d->Minv * d->u_drift;
  }

  // Computing the cost value and residuals
  costs_->calc(d->costs, x, u);
  d->cost = d->costs->cost;
}

template <typename Scalar>
void DifferentialActionModelFreeFwdDynamicsExtForcesTpl<Scalar>::calcDiff(
    const boost::shared_ptr<DifferentialActionDataAbstract>& data, const Eigen::Ref<const VectorXs>& x,
    const Eigen::Ref<const VectorXs>& u) {
  if (static_cast<std::size_t>(x.size()) != state_->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "x has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
  }
  if (static_cast<std::size_t>(u.size()) != nu_) {
    throw_pretty("Invalid argument: "
                 << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
  }

  const long& nv = (long)state_->get_nv();
  const Eigen::VectorBlock<const Eigen::Ref<const VectorXs>, Eigen::Dynamic> q = x.head((long)state_->get_nq());
  const Eigen::VectorBlock<const Eigen::Ref<const VectorXs>, Eigen::Dynamic> v = x.tail(nv);

  Data* d = static_cast<Data*>(data.get());

  // Expressing extforces_ (expressed in world frame at joint center) in the local frame of the joints
  ForceAlignedVector localFrameExtForces = extforces_;
  VectorXs testForce(6);
  testForce(0) = 0.;
  testForce(1) = 400.;
  testForce(2) = 0.;
  testForce(3) = 0.;
  testForce(4) = 0.;
  testForce(5) = 0.;
  //pinocchio::ForceTpl<Scalar>(testForce, Vector3s::Zero())

  for (pinocchio::JointIndex i = 1; i < (pinocchio::JointIndex)pinocchio_.njoints; ++i) {
    localFrameExtForces[i] = d->pinocchio.oMi[i].actInv(extforces_[i]);
    //localFrameExtForces[i] = pinocchio::ForceTpl<Scalar>(testForce);
  }
  //extforces_ = localFrameExtForces;

  actuation_->calcDiff(d->multibody.actuation, x, u);

  // Computing the dynamics derivatives
  if (with_armature_) {
    pinocchio::computeABADerivatives(pinocchio_, d->pinocchio, q, v, d->multibody.actuation->tau,
                                     localFrameExtForces, d->Fx.leftCols(nv),
                                     d->Fx.rightCols(nv), d->pinocchio.Minv);
    d->Fx.noalias() += d->pinocchio.Minv * d->multibody.actuation->dtau_dx;
    d->Fu.noalias() = d->pinocchio.Minv * d->multibody.actuation->dtau_du;
  } else { //TODO: adapt this to take extforces_ into account
    pinocchio::computeRNEADerivatives(pinocchio_, d->pinocchio, q, v, d->xout);
    d->dtau_dx.leftCols(nv) = d->multibody.actuation->dtau_dx.leftCols(nv) - d->pinocchio.dtau_dq;
    d->dtau_dx.rightCols(nv) = d->multibody.actuation->dtau_dx.rightCols(nv) - d->pinocchio.dtau_dv;
    d->Fx.noalias() = d->Minv * d->dtau_dx;
    d->Fu.noalias() = d->Minv * d->multibody.actuation->dtau_du;
  }

  // Computing the cost derivatives
  costs_->calcDiff(d->costs, x, u);
}

template <typename Scalar>
boost::shared_ptr<DifferentialActionDataAbstractTpl<Scalar> >
DifferentialActionModelFreeFwdDynamicsExtForcesTpl<Scalar>::createData() {
  return boost::allocate_shared<Data>(Eigen::aligned_allocator<Data>(), this);
}

template <typename Scalar>
bool DifferentialActionModelFreeFwdDynamicsExtForcesTpl<Scalar>::checkData(
    const boost::shared_ptr<DifferentialActionDataAbstract>& data) {
  boost::shared_ptr<Data> d = boost::dynamic_pointer_cast<Data>(data);
  if (d != NULL) {
    return true;
  } else {
    return false;
  }
}
template <typename Scalar>
void DifferentialActionModelFreeFwdDynamicsExtForcesTpl<Scalar>::quasiStatic(
    const boost::shared_ptr<DifferentialActionDataAbstract>& data, Eigen::Ref<VectorXs> u,
    const Eigen::Ref<const VectorXs>& x, const std::size_t&, const Scalar&) {
  if (static_cast<std::size_t>(u.size()) != nu_) {
    throw_pretty("Invalid argument: "
                 << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
  }
  if (static_cast<std::size_t>(x.size()) != state_->get_nx()) {
    throw_pretty("Invalid argument: "
                 << "x has wrong dimension (it should be " + std::to_string(state_->get_nx()) + ")");
  }
  // Static casting the data
  Data* d = static_cast<Data*>(data.get());
  const Eigen::VectorBlock<const Eigen::Ref<const VectorXs>, Eigen::Dynamic> q = x.head((long)state_->get_nq());

  // Check the velocity input is zero
  assert_pretty(x.tail((long)state_->get_nv()).isZero(), "The velocity input should be zero for quasi-static to work.");

  // Expressing extforces_ (expressed in world frame at joint center) in the local frame of the joints
  ForceAlignedVector localFrameExtForces = extforces_;
  VectorXs testForce(6);
  testForce(0) = 0.;
  testForce(1) = 400.;
  testForce(2) = 0.;
  testForce(3) = 0.;
  testForce(4) = 0.;
  testForce(5) = 0.;
  for (pinocchio::JointIndex i = 1; i < (pinocchio::JointIndex)pinocchio_.njoints; ++i) {
    localFrameExtForces[i] = d->pinocchio.oMi[i].actInv(extforces_[i]);
    //localFrameExtForces[i] = pinocchio::ForceTpl<Scalar>(testForce);
  }
  //extforces_ = localFrameExtForces;

  d->pinocchio.tau =
      pinocchio::rnea(pinocchio_, d->pinocchio, q, VectorXs::Zero((long)state_->get_nv()),
                      VectorXs::Zero((long)state_->get_nv()), localFrameExtForces);

  d->tmp_xstatic.head((long)state_->get_nq()) = q;
  actuation_->calc(d->multibody.actuation, d->tmp_xstatic, VectorXs::Zero((long)nu_));
  actuation_->calcDiff(d->multibody.actuation, d->tmp_xstatic, VectorXs::Zero((long)nu_));

  u.noalias() = pseudoInverse(d->multibody.actuation->dtau_du) * d->pinocchio.tau;
  d->pinocchio.tau.setZero();
}

template <typename Scalar>
pinocchio::ModelTpl<Scalar>& DifferentialActionModelFreeFwdDynamicsExtForcesTpl<Scalar>::get_pinocchio() const {
  return pinocchio_;
}

template <typename Scalar>
const boost::shared_ptr<ActuationModelAbstractTpl<Scalar> >&
DifferentialActionModelFreeFwdDynamicsExtForcesTpl<Scalar>::get_actuation() const {
  return actuation_;
}

template <typename Scalar>
const boost::shared_ptr<CostModelSumTpl<Scalar> >& DifferentialActionModelFreeFwdDynamicsExtForcesTpl<Scalar>::get_costs()
    const {
  return costs_;
}

template <typename Scalar>
const typename MathBaseTpl<Scalar>::VectorXs& DifferentialActionModelFreeFwdDynamicsExtForcesTpl<Scalar>::get_armature() const {
  return armature_;
}

template <typename Scalar>
void DifferentialActionModelFreeFwdDynamicsExtForcesTpl<Scalar>::set_armature(const VectorXs& armature) {
  if (static_cast<std::size_t>(armature.size()) != state_->get_nv()) {
    throw_pretty("Invalid argument: "
                 << "The armature dimension is wrong (it should be " + std::to_string(state_->get_nv()) + ")");
  }

  armature_ = armature;
  with_armature_ = false;
}

}  // namespace example
}  // namespace gepetto
