//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "InterfaceFluxMatch.h"

// MOOSE includes
#include "MooseVariable.h"

#include <cmath>

registerMooseObject("ZapdosApp", InterfaceFluxMatch);

InputParameters
InterfaceFluxMatch::validParams()
{
  InputParameters params = ADInterfaceKernel::validParams();
  params.addRequiredCoupledVar("potential", "The potential on the primary side of the interface.");
  params.addRequiredCoupledVar("potential_neighbor",
                               "The potential on the secondary side of the interface.");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  params.addRequiredParam<Real>("neighbor_position_units",
                                "The units of position in the neighboring domain.");
  params.addClassDescription(
      "Used to include the electric field driven advective flux of species"
      "into or out of a neighboring subdomain. Currently this interface kernel"
      "is specific to electrons because the transport coefficients are assumed"
      "to be a function of the mean electron energy. A generic interface"
      "kernel with constant transport coefficients will have a much simpler Jacobian");
  return params;
}

InterfaceFluxMatch::InterfaceFluxMatch(const InputParameters & parameters)
  : ADInterfaceKernel(parameters),
    _r_units(1. / getParam<Real>("position_units")),
    _r_neighbor_units(1. / getParam<Real>("neighbor_position_units")),

    _grad_potential(adCoupledGradient("potential")),
    _potential_neighbor_var(*getVar("potential_neighbor", 0)),
    _grad_potential_neighbor(_potential_neighbor_var.adGradSlnNeighbor()),

    // Primary properties
    _mu(getADMaterialProperty<Real>("mu" + _var.name())),
    _diff(getADMaterialProperty<Real>("diff" + _var.name())),
    _sgn(getMaterialProperty<Real>("sgn" + _var.name())),

    // Neighbor properties
    _mu_neighbor(getNeighborADMaterialProperty<Real>("mu" + _neighbor_var.name())),
    _diff_neighbor(getNeighborADMaterialProperty<Real>("diff" + _neighbor_var.name())),
    _sgn_neighbor(getNeighborMaterialProperty<Real>("sgn" + _neighbor_var.name()))
{
}

ADReal
InterfaceFluxMatch::computeQpResidual(Moose::DGResidualType type)
{
  /*
  ADReal flux = std::exp(_u[_qp]) *
                (_sgn[_qp] * -_grad_potential[_qp] * _r_units - _diff[_qp] * _grad_u[_qp]) *
                _normals[_qp];
                */
  ADReal flux_neighbor = std::exp(_neighbor_value[_qp]) *
                         (_sgn_neighbor[_qp] * _mu_neighbor[_qp] * -_grad_potential_neighbor[_qp] -
                          _diff_neighbor[_qp] * _grad_neighbor_value[_qp]) *
                         _r_neighbor_units * _normals[_qp];
  ADReal r = 0;

  switch (type)
  {
    case Moose::Element:
      /*
      r = _mu_neighbor[_qp] * _sgn_neighbor[_qp] * -_grad_potential_neighbor[_qp] *
          _r_neighbor_units * std::exp(_neighbor_value[_qp]) * _normals[_qp] * _test[_i][_qp] *
          _r_units;
          */
      r = _test[_i][_qp] * _r_units * flux_neighbor;
      break;

    case Moose::Neighbor:
      r = 0.;
      // r = -_test_neighbor[_i][_qp] * _r_neighbor_units * flux;
      break;
  }

  return r;
}
