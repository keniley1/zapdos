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

#include "ADMaterial.h"
#include "DerivativeMaterialPropertyNameInterface.h"

class SurfaceChargeNew : public ADMaterial, public DerivativeMaterialPropertyNameInterface
{
public:
  static InputParameters validParams();
  SurfaceChargeNew(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;
  virtual void initQpStatefulProperties() override;

  void computeChargeFlux();

  void computeElectronFlux();

  int _bc_type;
  ADMaterialProperty<Real> & _sigma;
  const MaterialProperty<Real> & _sigma_old;
  MaterialProperty<Real> & _plasma_current;
  const bool _include_electrons;
  const bool _include_secondary_electrons;
  
  // Electron parameters
  const ADVariableValue * _electron;
  const ADVariableValue * _mean_energy;
  const MaterialProperty<Real> * _mass_em;
  const ADMaterialProperty<Real> * _muem;
  const MaterialProperty<Real> * _se_coeff;
  const MaterialProperty<Real> * _e;
  const MaterialProperty<Real> * _kb;

  Real _r_units;
  const std::string _potential_units;

  const ADVariableGradient & _grad_potential;

  unsigned int _num_species;
  std::vector<const ADVariableValue *> _species;
  std::vector<const ADVariableGradient *> _grad_species;
  std::vector<const ADMaterialProperty<Real> *> _mu;
  std::vector<const ADMaterialProperty<Real> *> _diff;
  std::vector<const MaterialProperty<Real> *> _sgn;
  std::vector<const MaterialProperty<Real> *> _mass;
  std::vector<const ADMaterialProperty<Real> *> _Tion;

  // Recombination coefficient for Lymberopoulos-style BCs
  Real _ks;

  ADReal _charge_flux;

  ADReal _vi_thermal;
  ADReal _ve_thermal;
  ADReal _electron_flux;
  ADReal _n_gamma;

  Real _voltage_scaling;
  Real _q_times_NA;
  Real _r_ion;
  Real _r_electron;
  Real _a;
  Real _b;
};
