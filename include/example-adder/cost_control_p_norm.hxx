///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2019-2020, LAAS-CNRS, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include <iostream> // WHY DO I NEED THAT ?!!

#include "crocoddyl/core/utils/exception.hpp"

#include "example-adder/cost_control_p_norm.hpp"

namespace gepetto {
namespace example {
template <typename Scalar>
CostModelControlPNormTpl<Scalar>::CostModelControlPNormTpl(boost::shared_ptr<typename Base::StateAbstract> state,
                                                 const std::size_t& nu,
                                                 const std::size_t& p)
    : Base(state, nu, nu), p_(p) { // specifying nu mandatory due to free-flyer state
  if (activation_->get_nr() != nu_) { // TODO unnecessary, right?
    throw_pretty("Invalid argument: "
                 << "nr is equals to " + std::to_string(nu_));
  }
  if (p < 2) {
    throw_pretty("Invalid argument: "
                 << "p must be >= 2; if p=2, prefer the CostModelControl with 2-norm smooth activation function -- OR NOT actually");            
  
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
  //std::cout << "r_tmp before + eps = " << r_tmp << std::endl;
  r_tmp += Scalar(1.); //TODO write epsilon param properly
  //std::cout << "r_tmp after plus eps = " << r_tmp << std::endl;
  //r_tmp = std::pow(r_tmp, (int(1)/p_)); // doesn't work apparently, returns 1.
  //std::cout << "r_tmp final = " << r_tmp << std::endl;
  data->r = u; // must be a vector, TODO think about putting: r = u - u_gravity_compensation
  // est-ce nécessaire pour activer les variables ?
  //activation_->calc(data->activation, data->r);
  //std::cout << "p = " << p_ << " 1/p_ = " << 1/p_ << std::endl;
  //std::cout << "1/float(p_) = " << 1/float(p_) << std::endl;
  data->cost = std::pow(r_tmp, 1/float(p_)); // float(p_) needed or integer part of quotient returned
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

  VectorXs abs_u_p_minus_two = VectorXs::Zero(nu_);
  for (std::size_t i = 0; i < nu_; ++i){
    abs_u_p_minus_two[i] = std::pow(std::abs(u[i]), (p_ - Scalar(2)))/(std::pow(data->cost, (p_ - Scalar(1))));
  }
  // cout << "Here is the matrix a:\n" << a << endl;

  MatrixXs hessian = MatrixXs::Zero(nu_, nu_); //TODO check calculus properly pasted from numdiff checking file
  for (std::size_t i = 0; i < nu_; ++i){
    for (std::size_t j = 0; j < nu_; ++j){
      if (i == j){
        //std::cout << "i = j = " << i << std::endl;
        Scalar long_bit = std::pow(std::abs(data->r[i]), int(2*p_ - 4));
        //std::cout << "int(2*p_ - 4) = " << int(2*p_ - 4) << std::endl;
        Scalar main_bit = (int(1 - p_))*std::pow(data->cost, int(1 - 2*p_))*std::pow(data->r[i], 2)*long_bit;
        //std::cout << "int(1 - p_) = " << int(1 - p_) << std::endl;
        Scalar second_term = std::pow(std::abs(data->r[i]), int(p_ - 2));
        //std::cout << "int(p_ - 2) = " << int(p_ - 2) << std::endl;
        Scalar third_term = (int(p_ - 2))*std::pow(std::abs(data->r[i]), int(p_ - 2));
        hessian(i,j) = main_bit + std::pow(data->cost, int(1 - p_))*(second_term + third_term);
        //std::cout << "INSIDE IF p_ = " << p_ << " AND 1 - 2*p_ = " << (1 - 2*p_) << std::endl;
        //std::cout << "int(1-2*p_) = " << int(1 - 2*p_) << std::endl; // THIS ONE GOOD syntax
        //std::cout << "int(1) - int(2*p_) = " << int(1) - int(2*p_) << std::endl; // this one good syntax as well
      }
      else{
        Scalar bit_in_i = data->r[i]*std::pow(std::abs(data->r[i]), int(p_ - 2));
        Scalar bit_in_j = data->r[j]*std::pow(std::abs(data->r[j]), int(p_ - 2));
        hessian(i,j) = (int(1 - p_))*std::pow(data->cost, int(1 - 2*p_))*bit_in_i*bit_in_j;
        //std::cout << "INSIDE ELSE int(1) - int(2)*p_ = " << int(1) - int(2)*p_ << std::endl;
        //std::cout << "i = " << i << " j = " << j << std::endl;
      }
    }
  }

  //activation_->calcDiff(data->activation, data->r);
  data->Lu = u.cwiseProduct(abs_u_p_minus_two);//data->activation->Ar;
  data->Luu = hessian; // TODO what do I do here? 2*Lu^T * Lu ?
}


}  // namespace example
}  // namespace gepetto
