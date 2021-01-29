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

#include "ADIntegratedBC.h"
#include "ProvideMobility.h"

// Forward Declarations
class IonHeating;

declareADValidParams(IonHeating);

class IonHeating : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  IonHeating(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  Real _r_units;
  std::vector<const ADVariableValue *> _ion;
  std::vector<const ADVariableGradient *> _grad_ion;
  std::vector<unsigned int> _ion_id;
  //const ADVariableValue & _mean_en;
  //const ADVariableValue & _em;

  //const MaterialProperty<Real> & _se_coeff;
  std::vector<const ADMaterialProperty<Real> *> _muion;
  //const MaterialProperty<Real> & _eps;
  //const MaterialProperty<Real> & _N_A;
  std::vector<const MaterialProperty<Real> *> _sgnion;
  //std::vector<const ADMaterialProperty<Real> *> _Dion;
  //const ADMaterialProperty<Real> & _muem;
  //const MaterialProperty<Real> & _e;
  //const MaterialProperty<Real> & _massem;
  std::vector<const ADMaterialProperty<Real> *> _T_heavy;
  const MaterialProperty<Real> & _kb;
  const ADVariableGradient & _grad_potential;
  const MaterialProperty<Real> & _work_function;
  std::vector<const MaterialProperty<Real> *> _mass;

  std::string _potential_units;
  Real _r;

  //ADRealVectorValue _ion_flux;
  ADReal _ion_flux;
  //ADReal _n_gamma;
  //ADReal _actual_mean_en;
  //ADReal _v_e_th;
  ADReal _v_i_th;
  //Real _a;
  Real _b;
  //ADReal _numerator;
  //ADReal _denominator;

  Real _voltage_scaling;

  ADReal _ion_drift;
  //ADReal _secondary_ion;
  unsigned int _num_ions;
};
