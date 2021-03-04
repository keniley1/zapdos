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

class HeatFlux : public AuxKernel
{
public:
  static InputParameters validParams();

  HeatFlux(const InputParameters & parameters);

  virtual ~HeatFlux() {}

protected:
  virtual Real computeValue();

  Real _r_units;
  const unsigned _component;

  // Coupled variables
  const VariableValue & _temp;
  const VariableGradient & _grad_temp;
  const VariableValue & _u_vel;
  const VariableValue & _v_vel;
  const VariableValue & _w_vel;

  // Material properties

  const MaterialProperty<Real> & _rho;
  const MaterialProperty<Real> & _k;
  const MaterialProperty<Real> & _cp;
};
