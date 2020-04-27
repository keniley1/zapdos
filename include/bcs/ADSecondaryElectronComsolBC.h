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

// Forward Declarations
class ADSecondaryElectronComsolBC;

declareADValidParams(ADSecondaryElectronComsolBC);

class ADSecondaryElectronComsolBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  ADSecondaryElectronComsolBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  Real _r_units;
  Real _r;

  // Coupled variables

  const ADVariableGradient & _grad_potential;
  //const ADVariableValue & _mean_en;
  std::vector<const ADVariableValue *> _ip;
  std::vector<const ADVariableGradient *> _grad_ip;

  //const ADMaterialProperty<Real> & _muem;
  //const MaterialProperty<Real> & _massem;
  //const MaterialProperty<Real> & _e;
  std::vector<const MaterialProperty<Real> *> _sgnip;
  std::vector<const ADMaterialProperty<Real> *> _muip;
  std::vector<const ADMaterialProperty<Real> *> _Dip;
  //const MaterialProperty<Real> & _se_coeff;
  const Real _se_coeff;

  Real _a;
  //ADReal _v_thermal;
  ADRealVectorValue _ion_flux;
  //ADReal _n_gamma;

  unsigned int _num_ions;
};
