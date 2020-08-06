//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADInterfaceKernel.h"

/**
 * DG kernel for interfacing diffusion between two variables on adjacent blocks,
 * using the automatic differentiation system to calculate the Jacobian.
 */
class ADInterfaceAdvection : public ADInterfaceKernel
{
public:
  static InputParameters validParams();

  ADInterfaceAdvection(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual(Moose::DGResidualType type) override;

  Real _r_units;
  Real _r_neighbor_units;

  //MooseVariable & _potential_neighbor_var;
  const ADVariableGradient & _grad_potential_neighbor;

  const ADMaterialProperty<Real> & _mu_neighbor;
  const MaterialProperty<Real> & _sgn_neighbor;
};
