//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ChargeDensity.h"

registerMooseObject("ZapdosApp", ChargeDensity);

template <>
InputParameters
validParams<ChargeDensity>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("charged", "List of all charged particles.");
  params.addRequiredParam<std::string>("potential_units", "The potential units.");
  params.addClassDescription("Returns space charge density.");
  return params;
}

ChargeDensity::ChargeDensity(const InputParameters & parameters)
  : AuxKernel(parameters),

    _e(getMaterialProperty<Real>("e")),
    _N_A(getMaterialProperty<Real>("N_A")),
    _potential_units(getParam<std::string>("potential_units"))
{
  _num_ions = coupledComponents("charged");

  // Resize the vectors to store _num_ions values:
  _charged.resize(_num_ions);
  _sgn.resize(_num_ions);

  for (unsigned int i = 0; i < _num_ions; ++i)
  {
    _charged[i] = &coupledValue("charged", i);
    _sgn[i] = &getMaterialProperty<Real>("sgn" + (*getVar("charged", i)).name());
  }

  if (_potential_units.compare("V") == 0)
    _voltage_scaling = 1.;
  else if (_potential_units.compare("kV") == 0)
    _voltage_scaling = 1000;
}

Real
ChargeDensity::computeValue()
{
  //return -_test[_i][_qp] * _e[_qp] * _sgn[_qp] * _N_A[_qp] * std::exp(_charged[_qp]) /
  //       _voltage_scaling;
  _val = 0;

  for (unsigned int i = 0; i < _num_ions; ++i)
  {
    _val += std::exp((*_charged[i])[_qp]) * (*_sgn[i])[_qp];
  } 
  return _val * _e[_qp] * _voltage_scaling * _N_A[_qp];
}
