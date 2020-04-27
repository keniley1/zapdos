//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DoNothingPotentialBC_Y.h"

// MOOSE includes
#include "MooseVariable.h"

registerMooseObject("ZapdosApp", DoNothingPotentialBC_Y);

template <>
InputParameters
validParams<DoNothingPotentialBC_Y>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<Real>("position_units", "Units of position.");
  return params;
}

DoNothingPotentialBC_Y::DoNothingPotentialBC_Y(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _r_units(1. / getParam<Real>("position_units")),
    _diffusivity(getMaterialProperty<Real>("diff" + _var.name()))
{
}

DoNothingPotentialBC_Y::~DoNothingPotentialBC_Y() {}

Real
DoNothingPotentialBC_Y::computeQpResidual()
{
  return (_u[_qp] + _q_point[_qp](1) * _grad_u[_qp](1)) / _q_point[_qp](0) * _test[_i][_qp];
  // return _diffusivity[_qp] * _grad_u[_qp] * _r_units * -_normals[_qp] * _test[_i][_qp] *
  // _r_units;
}

Real
DoNothingPotentialBC_Y::computeQpJacobian()
{
  return (_phi[_j][_qp] + _q_point[_qp](1) * _grad_phi[_j][_qp](1)) / _q_point[_qp](0) *
         _test[_i][_qp];
  // return _diffusivity[_qp] * _grad_phi[_j][_qp] * _r_units * -_normals[_qp] * _test[_i][_qp] *
  // _r_units;
}
