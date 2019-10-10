//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef DONOTHINGENERGYBC_H
#define DONOTHINGENERGYBC_H

#include "IntegratedBC.h"

class DoNothingEnergyBC;

template <>
InputParameters validParams<DoNothingEnergyBC>();

// This diffusion kernel should only be used with species whose values are in the logarithmic form.

class DoNothingEnergyBC: public IntegratedBC
{
public:
  DoNothingEnergyBC(const InputParameters & parameters);
  virtual ~DoNothingEnergyBC();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  Real _r_units;

  const MaterialProperty<Real> & _muel;
  const MaterialProperty<Real> & _d_muel_d_actual_mean_en;
  const MaterialProperty<Real> & _sign;
  const MaterialProperty<Real> & _diffel;
  const MaterialProperty<Real> & _d_diffel_d_actual_mean_en;

  // Coupled variables
  unsigned int _potential_id;
  const VariableGradient & _grad_potential;
  const VariableValue & _em;
  unsigned int _em_id;

  Real _d_actual_mean_en_d_em;
  Real _d_muel_d_em;
  Real _d_actual_mean_en_d_u;
  Real _d_muel_d_u;
  Real _d_diffel_d_u;
  Real _d_diffel_d_em;
};

#endif /* DONOTHINGENERGYBC_H */
