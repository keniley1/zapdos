//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADSecondaryElectronComsolBC.h"

// MOOSE includes
#include "MooseVariable.h"

registerADMooseObject("ZapdosApp", ADSecondaryElectronComsolBC);

defineADLegacyParams(ADSecondaryElectronComsolBC);

InputParameters
ADSecondaryElectronComsolBC::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addRequiredParam<Real>("r", "The reflection coefficient");
  params.addRequiredCoupledVar("potential", "The electric potential");
  params.addRequiredCoupledVar("mean_en", "The electron density.");
  params.addRequiredCoupledVar("ip", "The ion density.");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  params.addParam<Real>("se_coeff", 1e-3, "The secondary electron coefficient. Default: 1e-3");
  return params;
}

ADSecondaryElectronComsolBC::ADSecondaryElectronComsolBC(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
    _r_units(1. / getParam<Real>("position_units")),
    _r(getParam<Real>("r")),

    // Coupled Variables
    _grad_potential(adCoupledGradient("potential")),
    //_mean_en(adCoupledValue("mean_en")),

    //_muem(getADMaterialProperty<Real>("muem")),
    //_massem(getMaterialProperty<Real>("massem")),
    //_e(getMaterialProperty<Real>("e"))
    //_se_coeff(getMaterialProperty<Real>("se_coeff"))
    _se_coeff(getParam<Real>("se_coeff"))
{
  _ion_flux.zero();
  _a = 0.5;

  _num_ions = coupledComponents("ip");

  // Resize the vectors to store _num_ions values:
  _ip.resize(_num_ions);
  _grad_ip.resize(_num_ions);
  _muip.resize(_num_ions);
  _Dip.resize(_num_ions);
  _sgnip.resize(_num_ions);

  // Retrieve the values for each ion and store in the relevant vectors.
  // Note that these need to be dereferenced to get the values inside the
  // main body of the code.
  // e.g. instead of "_ip[_qp]" it would be "(*_ip[i])[_qp]", where "i"
  // refers to a single ion species.
  for (unsigned int i = 0; i < _num_ions; ++i)
  {
    _ip[i] = &adCoupledValue("ip", i);
    _grad_ip[i] = &adCoupledGradient("ip", i);
    _muip[i] = &getADMaterialProperty<Real>("mu" + (*getVar("ip", i)).name());
    _Dip[i] = &getADMaterialProperty<Real>("diff" + (*getVar("ip", i)).name());
    _sgnip[i] = &getMaterialProperty<Real>("sgn" + (*getVar("ip", i)).name());
  }
}

ADReal
ADSecondaryElectronComsolBC::computeQpResidual()
{

  _ion_flux.zero();
  for (unsigned int i = 0; i < _num_ions; ++i)
  {
    if (_normals[_qp] * (*_sgnip[i])[_qp] * -_grad_potential[_qp] > 0.0)
    {
      _a = 1.0;
    }
    else
    {
      _a = 0.0;
    }

    _ion_flux += _a * (*_sgnip[i])[_qp] * (*_muip[i])[_qp] * -_grad_potential[_qp] * _r_units *
                     std::exp((*_ip[i])[_qp]) -
                 (*_Dip[i])[_qp] * std::exp((*_ip[i])[_qp]) * (*_grad_ip[i])[_qp] * _r_units;
  }

  return -_test[_i][_qp] * _se_coeff * _ion_flux * _normals[_qp];
}
