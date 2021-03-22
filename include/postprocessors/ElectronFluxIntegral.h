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

// MOOSE includes
#include "SideIntegralVariablePostprocessor.h"

// Forward Declarations
//class ElectronFluxIntegral;

/*
template <>
InputParameters validParams<ElectronFluxIntegral>();
*/
template <bool>
class ElectronFluxIntegralTempl;
typedef ElectronFluxIntegralTempl<false> ElectronFluxIntegral;
typedef ElectronFluxIntegralTempl<true> ADElectronFluxIntegral;

template <>
InputParameters validParams<ElectronFluxIntegral>();

/**
 * This postprocessor computes a side integral of the mass flux.
 */
template <bool is_ad>
class ElectronFluxIntegralTempl : public SideIntegralVariablePostprocessor
{
public:
  static InputParameters validParams();
  ElectronFluxIntegralTempl(const InputParameters & parameters);

protected:
  //virtual ~ElectronFluxIntegralTempl() {}
  virtual Real computeQpIntegral() override;

  Real _r_units;
  Real _r;
  const MaterialProperty<Real> & _kb;
  const MaterialProperty<Real> & _e;
  const MaterialProperty<Real> & _massem;
  Real _v_thermal;
  Real _a;

  const GenericMaterialProperty<Real, is_ad> & _muem;
  const VariableValue & _mean_en;
  const VariableGradient & _grad_potential;
};


