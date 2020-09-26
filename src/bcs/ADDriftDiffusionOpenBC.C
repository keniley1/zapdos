//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADDriftDiffusionOpenBC.h"

registerADMooseObject("ZapdosApp", ADDriftDiffusionOpenBC);

defineADLegacyParams(ADDriftDiffusionOpenBC);

InputParameters
ADDriftDiffusionOpenBC::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addRequiredCoupledVar("potential", "The electric potential");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  params.addClassDescription("Do nothing boundary condition for electrons, ions, and mean electron "
                             "energy. Acts as an open outflow BC.");
  return params;
}

ADDriftDiffusionOpenBC::ADDriftDiffusionOpenBC(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
    _r_units(1. / getParam<Real>("position_units")),

    // Coupled Variables
    _grad_potential(adCoupledGradient("potential")),

    _sign(getMaterialProperty<Real>("sgn" + _var.name())),
    _mu(getADMaterialProperty<Real>("mu" + _var.name())),
    _diffusivity(getADMaterialProperty<Real>("diff" + _var.name())),
    _e(getMaterialProperty<Real>("e"))
{
}

ADReal
ADDriftDiffusionOpenBC::computeQpResidual()
{
  return _mu[_qp] * _sign[_qp] * std::exp(_u[_qp]) * -_grad_potential[_qp] * _r_units *
             _normals[_qp] * _test[_i][_qp] * _r_units -
         _diffusivity[_qp] * std::exp(_u[_qp]) * _grad_u[_qp] * _r_units * _normals[_qp] *
             _test[_i][_qp] * _r_units;
}
