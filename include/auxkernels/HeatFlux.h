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

#include "AuxKernel.h"

template <bool is_ad>
class HeatFluxTempl : public AuxKernel
{
public:
  HeatFluxTempl(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeValue();

  Real _r_units;

  // Coupled variables
  const VariableValue & _temp;
  const VariableGradient & _grad_temp;
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;

  // Material properties

  const GenericMaterialProperty<Real, is_ad> & _rho;
  const GenericMaterialProperty<Real, is_ad> & _k;
  const GenericMaterialProperty<Real, is_ad> & _cp;
};

typedef HeatFluxTempl<false> HeatFlux;
typedef HeatFluxTempl<true> ADHeatFlux;
