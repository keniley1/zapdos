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

template <ComputeStage compute_stage>
class ADHeavySpeciesMaterial : public ADMaterial<compute_stage>
{
public:
  static InputParameters validParams();
  ADHeavySpeciesMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  Real _user_massHeavy;
  Real _user_sgnHeavy;
  std::string _potential_units;
  Real _voltage_scaling;

  ADMaterialProperty(Real) & _massHeavy;        // Replaces _massArp
  ADMaterialProperty(Real) & _temperatureHeavy; // Replaces _tempArp
  ADMaterialProperty(Real) & _sgnHeavy;         // Replaces _sgnArp (unused though)
  ADMaterialProperty(Real) & _muHeavy;          // Replaces _muArp
  ADMaterialProperty(Real) & _diffHeavy;        // Replaces _diffArp

  const MaterialProperty<Real> & _T_gas;
  const MaterialProperty<Real> & _p_gas;

  Real _time_units;
  bool _calc_mobility;
  bool _calc_diffusivity;
  ADMaterialProperty(RealVectorValue) & _grad_mu;
  ADMaterialProperty(RealVectorValue) & _grad_diff;

  usingMaterialMembers;
};
