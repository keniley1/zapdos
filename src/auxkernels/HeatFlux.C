//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeatFlux.h"
#include "MooseMesh.h"

registerMooseObject("ZapdosApp", HeatFlux);

registerMooseObject("ZapdosApp", HeatFlux);

InputParameters
HeatFlux::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", "y-velocity");
  params.addCoupledVar("w", "z-velocity");
  params.addRequiredCoupledVar("temperature", "The variable representing the log of the density.");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  params.addRequiredParam<unsigned>("component", "The direction.");
  params.addClassDescription("Returns the heat flux of defined species");

  return params;
}

HeatFlux::HeatFlux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _r_units(1. / getParam<Real>("position_units")),
    _component(getParam<unsigned>("component")),

    // Coupled variables
    _temp(coupledValue("temperature")),
    _grad_temp(coupledGradient("temperature")),

    _u_vel(coupledValue("u")),
    _v_vel(_mesh.dimension() >= 2 ? coupledValue("v") : _zero),
    _w_vel(_mesh.dimension() == 3 ? coupledValue("w") : _zero),

    // Material properties
    _rho(getMaterialProperty<Real>("rho")),
    _k(getMaterialProperty<Real>("k")),
    _cp(getMaterialProperty<Real>("cp"))
{
}

Real
HeatFlux::computeValue()
{
  /*
  return _rho[_qp] * _cp[_qp] *
             (_u_vel[_qp] * _grad_temp[_qp](0) + _v_vel[_qp] * _grad_temp[_qp](1) +
              _w_vel[_qp] * _grad_temp[_qp](2)) -
         _k[_qp] * _grad_temp[_qp];
         */
  return -_k[_qp] * _grad_temp[_qp](_component);
  /*
  return -raw_value(_diff[_qp]) * std::exp(_density_log[_qp]) * _grad_density_log[_qp](0) *
         _r_units * 6.02e23;
         */
}
