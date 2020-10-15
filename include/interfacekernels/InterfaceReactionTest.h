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
class InterfaceReactionTest : public ADInterfaceKernel
{
public:
  static InputParameters validParams();

  InterfaceReactionTest(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual(Moose::DGResidualType type);

  Real _r_units;
  Real _r_neighbor_units;
  Real _kf;
  Real _kb;

  const ADMaterialProperty<Real> & _diff;
  const ADMaterialProperty<Real> & _diff_neighbor;

  ADReal num;
  ADReal _diff_ij;

  //MooseVariable & _mean_en_neighbor_var;
  //const VariableValue & _mean_en_neighbor;
};
