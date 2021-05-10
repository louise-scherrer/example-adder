///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2019-2020, LAAS-CNRS, University of Edinburgh
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#ifndef __example_adder_cost_control_p_norm__
#define __example_adder_cost_control_p_norm__

#include "crocoddyl/core/fwd.hpp"
#include "crocoddyl/core/cost-base.hpp"
//#include "crocoddyl/core/utils/deprecate.hpp"

#include <pinocchio/spatial/force.hpp> // needed for example-adder/fwd.hpp

#include "example-adder/fwd.hpp"

namespace gepetto {
namespace example {

using namespace crocoddyl;

/**
 * @brief Cost on the p-norm of the control vector
 *
 * TODO write properly
 * This cost function defines a residual vector as \f$\mathbf{r}=\mathbf{u}-\mathbf{u}^*\f$, where
 * \f$\mathbf{u}\in~\mathbb{R}^{nu}\f$ is the current control input and \f$\mathbf{p}\f$ is the value of the norm.
 * This cost doesn't allow any activation function as a p-norm of the control vector is computed directly from it inside the calc
 * function.
 *
 * Both cost and residual derivatives are computed analytically. TODO do I compute the Hessian analytically as well?
 * For the computation of the cost Hessian, we use the Gauss-Newton approximation, e.g.
 * \f$\mathbf{l_{xx}} = \mathbf{l_{x}}^T \mathbf{l_{x}} \f$.
 *
 * As described in CostModelAbstractTpl(), the cost value and its derivatives are calculated by `calc` and `calcDiff`,
 * respectively.
 *
 * \sa `CostModelAbstractTpl`, `calc()`, `calcDiff()`, `createData()`
 */
template <typename _Scalar>
class CostModelControlPNormTpl : public CostModelAbstractTpl<_Scalar> {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef _Scalar Scalar;
  typedef MathBaseTpl<Scalar> MathBase;
  typedef CostModelAbstractTpl<Scalar> Base;
  typedef CostDataAbstractTpl<Scalar> CostDataAbstract;
  typedef ActivationModelAbstractTpl<Scalar> ActivationModelAbstract; // It is kept becaus cost-sum.hxx is using its nr
  // parameter, but not used in the computations, p-norm used instead
  typedef typename MathBase::VectorXs VectorXs;
  typedef typename MathBase::MatrixXs MatrixXs;

  /**
   * @brief Initialize the control cost model
   *
   * The default `nu` value is obtained from `StateAbstractTpl::get_nv()`.
   *
   * @param[in] state       State of the multibody system
   * @param[in] nu          Dimension of the control vector
   * @param[in] p           Value of the p-norm
   */
  CostModelControlPNormTpl(boost::shared_ptr<typename Base::StateAbstract> state, const std::size_t& nu,
                      const std::size_t& p);

  /**
   * @brief Compute the control cost
   *
   * @param[in] data  Control cost data
   * @param[in] x     State point \f$\mathbf{x}\in\mathbb{R}^{ndx}\f$
   * @param[in] u     Control input \f$\mathbf{u}\in\mathbb{R}^{nu}\f$
   */
  virtual void calc(const boost::shared_ptr<CostDataAbstract>& data, const Eigen::Ref<const VectorXs>& x,
                    const Eigen::Ref<const VectorXs>& u);

  /**
   * @brief Compute the derivatives of the control cost
   *
   * @param[in] data  Control cost data
   * @param[in] x     State point \f$\mathbf{x}\in\mathbb{R}^{ndx}\f$
   * @param[in] u     Control input \f$\mathbf{u}\in\mathbb{R}^{nu}\f$
   */
  virtual void calcDiff(const boost::shared_ptr<CostDataAbstract>& data, const Eigen::Ref<const VectorXs>& x,
                        const Eigen::Ref<const VectorXs>& u);

 protected:

  using Base::activation_;
  using Base::nu_;
  using Base::state_;
  using Base::unone_;

 private:
  const std::size_t p_;
};

}  // namespace example
}  // namespace gepetto

/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
#include "example-adder/cost_control_p_norm.hxx"

#endif  // CROCODDYL_CORE_COSTS_CONTROL_HPP_
