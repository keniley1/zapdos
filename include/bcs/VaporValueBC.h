//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADDirichletBCBase.h"

/**
 * Boundary condition of a Dirichlet type
 *
 * Sets the values of a nodal variable at nodes
 */
class VaporValueBC : public ADDirichletBCBase
{
public:
  static InputParameters validParams();

  VaporValueBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpValue() override;

  /// The value for this BC
  //const ADReal & _value;
  const ADVariableValue & _gas_temp;
  
  ADReal _vapor_pressure;
};
