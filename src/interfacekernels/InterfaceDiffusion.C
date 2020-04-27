//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceDiffusion.h"

registerMooseObject("ZapdosApp", InterfaceDiffusion);

InputParameters
InterfaceDiffusion::validParams()
{
  InputParameters params = InterfaceKernel::validParams();
  params.addParam<MaterialPropertyName>("diff", "D", "The diffusion coefficient.");
  params.addParam<MaterialPropertyName>(
      "D_neighbor", "D_neighbor", "The neighboring diffusion coefficient.");
  params.addRequiredParam<Real>("position_units", "The units of position.");
  return params;
}

InterfaceDiffusion::InterfaceDiffusion(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _r_units(1. / getParam<Real>("position_units")),
    _D(getMaterialProperty<Real>("diff" + _var.name())),
    _D_neighbor(getNeighborMaterialProperty<Real>("diff" + _neighbor_var.name())),
    _sigma(getMaterialProperty<Real>("surface_charge"))
{
}

Real
InterfaceDiffusion::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;

  switch (type)
  {
    case Moose::Element:
      r = -_test[_i][_qp] *
          (_D_neighbor[_qp] * _grad_neighbor_value[_qp] * _r_units * _normals[_qp] +
           _sigma[_qp]);
      break;

    case Moose::Neighbor:
      r = _test_neighbor[_i][_qp] *
          (_D[_qp] * _grad_u[_qp] * _r_units * _normals[_qp] + _sigma[_qp]);
      break;
  }

  return r;
}

Real
InterfaceDiffusion::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;

  switch (type)
  {
    case Moose::ElementElement:
    case Moose::NeighborNeighbor:
      break;

    case Moose::NeighborElement:
      jac = _test_neighbor[_i][_qp] * _D[_qp] * _grad_phi[_j][_qp] * _r_units * _normals[_qp];
      break;

    case Moose::ElementNeighbor:
      jac = -_test[_i][_qp] * _D_neighbor[_qp] * _grad_phi_neighbor[_j][_qp] * _r_units *
            _normals[_qp];
      break;
  }

  return jac;
}

