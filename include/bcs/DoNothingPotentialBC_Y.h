//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DONOTHINGPOTENTIALBC_Y_H
#define DONOTHINGPOTENTIALBC_Y_H

#include "IntegratedBC.h"

class DoNothingPotentialBC_Y;

template <>
InputParameters validParams<DoNothingPotentialBC_Y>();

// This diffusion kernel should only be used with species whose values are in the logarithmic form.

class DoNothingPotentialBC_Y: public IntegratedBC
{
public:
  DoNothingPotentialBC_Y(const InputParameters & parameters);
  virtual ~DoNothingPotentialBC_Y();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  Real _r_units;

  const MaterialProperty<Real> & _diffusivity;
};

#endif /* DONOTHINGPOTENTIALBC_Y_H */
