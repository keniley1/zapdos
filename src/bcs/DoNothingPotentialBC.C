//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DoNothingPotentialBC.h"

// MOOSE includes
#include "MooseVariable.h"

registerMooseObject("ZapdosApp", DoNothingPotentialBC);

template <>
InputParameters
validParams<DoNothingPotentialBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<Real>("position_units", "Units of position.");
  return params;
}

DoNothingPotentialBC::DoNothingPotentialBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _r_units(1. / getParam<Real>("position_units")),
    _diffusivity(getMaterialProperty<Real>("diff" + _var.name()))
{
}

DoNothingPotentialBC::~DoNothingPotentialBC() {}

Real
DoNothingPotentialBC::computeQpResidual()
{
  return _diffusivity[_qp] * _grad_u[_qp] * _r_units * -_normals[_qp] * _test[_i][_qp] * _r_units; 
}

Real
DoNothingPotentialBC::computeQpJacobian()
{
  return _diffusivity[_qp] * _grad_phi[_j][_qp] * _r_units * -_normals[_qp] * _test[_i][_qp] * _r_units;
}
