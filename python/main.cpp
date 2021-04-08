#include "example-adder/python.hpp"

BOOST_PYTHON_MODULE(example_adder) {
  gepetto::example::exposeExampleAdder();
  gepetto::example::exposeActionUniEx();
  gepetto::example::exposeDifferentialActionFreeFwdDynamicsExtForces();
  gepetto::example::exposeCostControlPNorm();
}
