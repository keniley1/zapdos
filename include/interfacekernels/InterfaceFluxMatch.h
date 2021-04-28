//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#pragma once

#include "ADInterfaceKernel.h"

/**
 * DG kernel for interfacing advection on adjacent blocks
 */
class InterfaceFluxMatch : public ADInterfaceKernel
{
public:
  static InputParameters validParams();

  InterfaceFluxMatch(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual(Moose::DGResidualType type) override;

  Real _r_units;
  Real _r_neighbor_units;

  const ADVariableGradient & _grad_potential;
  MooseVariable & _potential_neighbor_var;
  const ADVariableGradient & _grad_potential_neighbor;

  const ADMaterialProperty<Real> & _mu;
  const ADMaterialProperty<Real> & _diff;
  const MaterialProperty<Real> & _sgn;

  const ADMaterialProperty<Real> & _mu_neighbor;
  const ADMaterialProperty<Real> & _diff_neighbor;
  const MaterialProperty<Real> & _sgn_neighbor;
};
