//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElectronFluxIntegral.h"

#include "metaphysicl/raw_type.h"

using MetaPhysicL::raw_value;

registerMooseObject("ZapdosApp", ElectronFluxIntegral);
registerMooseObject("ZapdosApp", ADElectronFluxIntegral);

defineLegacyParams(ElectronFluxIntegral);

/*
template <is_ad>
InputParameters
validParams<ElectronFluxIntegral>()
*/
template <bool is_ad>
InputParameters
ElectronFluxIntegralTempl<is_ad>::validParams()
{
  //InputParameters params = validParams<SideIntegralVariablePostprocessor>();
  InputParameters params = SideIntegralVariablePostprocessor::validParams();
  params.addRequiredCoupledVar("potential", "The potential that drives the advective flux.");
  params.addRequiredCoupledVar("mean_energy", "The mean electron energy.");
  params.addRequiredParam<Real>("r", "The reflection coefficient");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  params.addClassDescription("Returns the flux of the electrons at the specified boundary.");
  return params;
}

template <bool is_ad>
ElectronFluxIntegralTempl<is_ad>::ElectronFluxIntegralTempl(const InputParameters & parameters)
  : SideIntegralVariablePostprocessor(parameters),
    _r_units(1. / getParam<Real>("position_units")),
    _r(getParam<Real>("r")),
    _kb(getMaterialProperty<Real>("k_boltz")),
    _e(getMaterialProperty<Real>("e")),
    _massem(getMaterialProperty<Real>("massem")),
    _v_thermal(0),
    _a(0.5),

    _muem(getGenericMaterialProperty<Real, is_ad>("muem")),
    _mean_en(coupledValue("mean_energy")),
    _grad_potential(coupledGradient("potential"))
{
}

template <bool is_ad>
Real
ElectronFluxIntegralTempl<is_ad>::computeQpIntegral()
{
  // Output units for base case are: mol / (m^2 * s)

  _v_thermal =
      std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(raw_value(_mean_en[_qp] - _u[_qp])) / (M_PI * _massem[_qp]));

  if (_normals[_qp] * _grad_potential[_qp] > 0.0)
    _a = 1.0;
  else
    _a = 0.0;

  return raw_value((1. - _r) / (1. + _r) *
         (0.5 * _v_thermal +
          ((1 - 2 * _a) * _muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp])) *
         std::exp(_u[_qp]));
}

template class ElectronFluxIntegralTempl<false>;
template class ElectronFluxIntegralTempl<true>;
