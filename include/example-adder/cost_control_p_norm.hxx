///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2019-2020, LAAS-CNRS, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include "crocoddyl/core/utils/exception.hpp"

#include "example-adder/cost_control_p_norm.hpp"

namespace gepetto {
namespace example {
template <typename Scalar>
CostModelControlPNormTpl<Scalar>::CostModelControlPNormTpl(boost::shared_ptr<typename Base::StateAbstract> state,
                                                 const std::size_t& p)
    : Base(state, state->get_nv()), p_(p) {
  if (activation_->get_nr() != nu_) { // TODO unnecessary, right?
    throw_pretty("Invalid argument: "
                 << "nr is equals to " + std::to_string(nu_));
  }
  if (p < 2) {
    throw_pretty("Invalid argument: "
                 << "p must be >= 2; if p=2, prefer the CostModelControl with quadratic activation function");            
  
  }
}

template <typename Scalar>
void CostModelControlPNormTpl<Scalar>::calc(const boost::shared_ptr<CostDataAbstract>& data,
                                       const Eigen::Ref<const VectorXs>&, const Eigen::Ref<const VectorXs>& u) {
  if (nu_ == 0) {
    throw_pretty("Invalid argument: "
                 << "it seems to be an autonomous system, if so, don't add this cost function");
  }
  if (static_cast<std::size_t>(u.size()) != nu_) {
    throw_pretty("Invalid argument: "
                 << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
  }
  Scalar r_tmp = 0.;
  for (std::size_t i = 0; i < nu_; ++i){
    r_tmp += std::pow(std::abs(u[i]), p_);    
  }
  r_tmp = std::pow(r_tmp, 1/p_);
  
  data->r = u; // must be a vector, TODO think about putting: r = u - u_gravity_compensation
  // est-ce nécessaire pour activer les variables ?
  //activation_->calc(data->activation, data->r);
  data->cost = r_tmp;
}

template <typename Scalar>
void CostModelControlPNormTpl<Scalar>::calcDiff(const boost::shared_ptr<CostDataAbstract>& data,
                                           const Eigen::Ref<const VectorXs>&, const Eigen::Ref<const VectorXs>& u) {
  if (nu_ == 0) {
    throw_pretty("Invalid argument: "
                 << "it seems to be an autonomous system, if so, don't add this cost function");
  }
  if (static_cast<std::size_t>(u.size()) != nu_) {
    throw_pretty("Invalid argument: "
                 << "u has wrong dimension (it should be " + std::to_string(nu_) + ")");
  }

  VectorXs u_p_minus_two = VectorXs::Zero(nu_);
  for (std::size_t i = 0; i < nu_; ++i){
    u_p_minus_two[i] = std::pow(std::abs(u[i]), (p_ - 2))/(std::pow(data->cost, (p_ - 1)));
  }
  // cout << "Here is the matrix a:\n" << a << endl;

  //activation_->calcDiff(data->activation, data->r);
  data->Lu = u.cwiseProduct(u_p_minus_two);//data->activation->Ar;
  data->Luu = Scalar(2.) * data->Lu.transpose() * data->Lu; // TODO what do I do here? 2*Lu^T * Lu ?
}


}  // namespace example
}  // namespace gepetto
