//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "StressContinuity.h"

// MOOSE includes
#include "MooseVariable.h"
#include "MooseMesh.h"

#include <cmath>

registerMooseObject("ZapdosApp", StressContinuity);

InputParameters
StressContinuity::validParams()
{
  InputParameters params = ADInterfaceKernel::validParams();
  params.addRequiredCoupledVar("vel_x", "x-component of the velocity.");
  params.addCoupledVar("vel_y", "y-component of the velocity.");
  params.addRequiredCoupledVar("neighbor_vel_x", "neighboring x-component of the velocity.");
  params.addCoupledVar("neighbor_vel_y", "neighboring y-component of the velocity.");
  params.addRequiredParam<unsigned>("component", "The velocity component that this is applied to.");
  params.addCoupledVar("pressure", 0, "The pressure variable.");
  params.addCoupledVar("pressure_neighbor", 0, "The neighbor pressure variable.");
  params.addClassDescription(
      "Used to include the diffusive flux of species into or out of a neighboring"
      "subdomain. Currently specific to electrons.");
  return params;
}

StressContinuity::StressContinuity(const InputParameters & parameters)
  : ADInterfaceKernel(parameters),
    _component(getParam<unsigned>("component")),
    _grad_vel_x(adCoupledGradient("vel_x")),
    _grad_vel_y(adCoupledGradient("vel_y")),
    _neighbor_vel_x_var(*getVar("neighbor_vel_x", 0)),
    _neighbor_vel_y_var(*getVar("neighbor_vel_y", 0)),
    _grad_neighbor_vel_x(_neighbor_vel_x_var.adGradSlnNeighbor()),
    _grad_neighbor_vel_y(_neighbor_vel_y_var.adGradSlnNeighbor()),
    _p(adCoupledValue("pressure")),
    _p_neighbor(adCoupledNeighborValue("pressure_neighbor")),
    _mu(getMaterialProperty<Real>("mu")),
    _mu_neighbor(getNeighborMaterialProperty<Real>("mu"))

//_grad_potential_neighbor(_potential_neighbor_var.adGradSlnNeighbor()),
{
  if (!parameters.isParamValid("boundary"))
  {
    mooseError("In order to use the StressContinuity dgkernel, you must specify a "
               "boundary where it will live.");
  }
  //RealVectorValue _I(0, 0, 0);
  _I.zero();
  _I(_component) = 1.0;
}

ADReal
StressContinuity::computeQpResidual(Moose::DGResidualType type)
{
  ADReal r;
  r = 0;

  if (_component == 0)
  {
    _stress(0) = _p[_qp] - 2.0 * _mu[_qp] * _grad_vel_x[_qp](0);
    _stress(1) = -_mu[_qp] * (_grad_vel_x[_qp](1) + _grad_vel_y[_qp](0));
    _stress(2) = 0;


    _neighbor_stress(0) = _p_neighbor[_qp] - 2.0 * _mu_neighbor[_qp] * _grad_neighbor_vel_x[_qp](0);
    _neighbor_stress(1) = -_mu_neighbor[_qp] * (_grad_neighbor_vel_x[_qp](1) + _grad_neighbor_vel_y[_qp](0));
    _neighbor_stress(2) = 0;
  }
  else if (_component == 1)
  {
    _stress(0) = -_mu[_qp] * (_grad_vel_x[_qp](1) + _grad_vel_y[_qp](0));
    _stress(1) = _p[_qp] - 2.0 * _mu[_qp] * _grad_vel_y[_qp](1);
    _stress(2) = 0;


    _neighbor_stress(0) = -_mu_neighbor[_qp] * (_grad_neighbor_vel_x[_qp](1) + _grad_neighbor_vel_y[_qp](0));
    _neighbor_stress(1) = _p_neighbor[_qp] - 2.0 * _mu_neighbor[_qp] * _grad_neighbor_vel_y[_qp](1);
    _neighbor_stress(2) = 0;
  }

  switch (type)
  {
    case Moose::Element:
      r = _test[_i][_qp] * _neighbor_stress * _normals[_qp];
      break;

    case Moose::Neighbor:
      r = -_test_neighbor[_i][_qp] * _stress * _normals[_qp];
      break;
  }

  return r;
}
