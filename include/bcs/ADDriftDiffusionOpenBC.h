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
class ADDriftDiffusionOpenBC;

declareADValidParams(ADDriftDiffusionOpenBC);

class ADDriftDiffusionOpenBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();
  ADDriftDiffusionOpenBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  Real _r_units;

  // Coupled variables

  const ADVariableGradient & _grad_potential;

  const MaterialProperty<Real> & _sign;
  const ADMaterialProperty<Real> & _mu;
  const ADMaterialProperty<Real> & _diffusivity;
  const MaterialProperty<Real> & _e;
};
