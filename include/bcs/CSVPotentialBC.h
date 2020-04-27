#pragma once

#include "NodalBC.h"
#include "SplineInterpolation.h"

/**
 * Implements the Dirichlet boundary condition
 * c*u + u^2 + v^2 = _value
 * where "u" is the current variable, and "v" is a coupled variable.
 * Note: without the constant term, a zero initial guess gives you a
 * zero row in the Jacobian, which is a bad thing.
 */
class CSVPotentialBC : public NodalBC
{
public:
  static InputParameters validParams();

  CSVPotentialBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  //virtual Real computeQpJacobian();
  //virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  // The coupled variable
  SplineInterpolation _voltage;
};
