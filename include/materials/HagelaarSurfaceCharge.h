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

#include "Material.h"

class HagelaarSurfaceCharge;

template <>
InputParameters validParams<HagelaarSurfaceCharge>();

class HagelaarSurfaceCharge : public Material
{
public:
  HagelaarSurfaceCharge(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

  MaterialProperty<Real> & _sigma;
  const MaterialProperty<Real> & _sigma_old;
  //MaterialProperty<std::vector<Real>> & _d_ion_flux_d_ions;

  const MaterialProperty<Real> & _muem;
  //const MaterialProperty<Real> & _diffem;

  Real _r_units;
  Real _r_ion;
  Real _r_electron;
  const MaterialProperty<Real> & _kb;
  Real _v_thermal;
  Real _ve_thermal;
  Real _user_velocity;
  const MaterialProperty<Real> & _e;
  Real _se_coeff;
  const std::string _potential_units;
  Real _a;
  Real _b;

  const VariableGradient & _grad_potential;
  const VariableValue & _mean_en;
  const VariableValue & _em;
  const VariableGradient & _grad_em;

  unsigned int _num_ions;
  std::vector<const VariableValue *> _ions;
  std::vector<const VariableGradient *> _grad_ions;
  std::vector<const MaterialProperty<Real> *> _mu_ions;
  //std::vector<const MaterialProperty<Real> *> _diff_ions;
  std::vector<const MaterialProperty<Real> *> _sgn_ions;
  std::vector<const MaterialProperty<Real> *> _mass_ions;
  std::vector<const MaterialProperty<Real> *> _T_ions;

  Real _ion_flux;
  Real _n_gamma;

  Real _electron_flux;
  Real _voltage_scaling;
};

