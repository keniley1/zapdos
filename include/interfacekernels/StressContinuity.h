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
 * DG kernel for interfacing diffusion between two variables on adjacent blocks
 */
class StressContinuity : public ADInterfaceKernel
{
public:
  static InputParameters validParams();

  StressContinuity(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual(Moose::DGResidualType type);

  unsigned _component;

  const ADVariableGradient & _grad_vel_x; 
  const ADVariableGradient & _grad_vel_y; 

  MooseVariable & _neighbor_vel_x_var;
  MooseVariable & _neighbor_vel_y_var;
  const ADVariableGradient & _grad_neighbor_vel_x; 
  const ADVariableGradient & _grad_neighbor_vel_y; 

  const ADVariableValue & _p;
  const ADVariableValue & _p_neighbor;

  const MaterialProperty<Real> & _mu;
  const MaterialProperty<Real> & _mu_neighbor;

  RealVectorValue _I;

  ADRealVectorValue _stress;
  ADRealVectorValue _neighbor_stress;
};
