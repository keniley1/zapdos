#include "InterfaceFluxConservation.h"
#include "MooseVariable.h"
#include <cmath>

registerMooseObject("ZapdosApp", InterfaceFluxConservation);

template <>
InputParameters
validParams<InterfaceFluxConservation>()
{
  InputParameters params = validParams<InterfaceKernel>();
  params.addRequiredParam<std::string>("region_name", "The name of the primary region.");
  params.addRequiredParam<std::string>("neighbor_region_name",
                                       "The name of the neighboring region.");
  params.addRequiredParam<Real>("position_units",
                                "The position unit scaling of the primary region.");
  params.addRequiredParam<Real>("neighbor_position_units",
                                "The position unit scalin gof the neighboring region.");
  return params;
}

InterfaceFluxConservation::InterfaceFluxConservation(const InputParameters & parameters)
  : InterfaceKernel(parameters),
    _r_units(1. / getParam<Real>("position_units")),
    _r_neighbor_units(1. / getParam<Real>("neighbor_position_units")),
    _diffusivity_main(getMaterialProperty<Real>(getParam<std::string>("region_name"))),
    _diffusivity_neighbor(getMaterialProperty<Real>(getParam<std::string>("neighbor_region_name")))
{
  if (!parameters.isParamValid("boundary"))
    mooseError("In order to use the InterfaceFluxConservation DGKernel, you must specify a "
               "boundary on which it will live.");
}

Real
InterfaceFluxConservation::computeQpResidual(Moose::DGResidualType type)
{
  Real r = 0;
  switch (type)
  {
    case Moose::Element:
      r = -_diffusivity_neighbor[_qp] * _grad_neighbor_value[_qp] * _normals[_qp] * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      r = _diffusivity_main[_qp] * _grad_u[_qp] * _normals[_qp] * _test_neighbor[_i][_qp];
      break;
  }
  return r;
  //mooseError("Internal Error. (InterfaceFluxConservation)");
}

Real
InterfaceFluxConservation::computeQpJacobian(Moose::DGJacobianType type)
{
  Real jac = 0;

  switch (type)
  {
    case Moose::ElementElement:
    case Moose::NeighborNeighbor:
      break;

    case Moose::NeighborElement:
      jac = _test_neighbor[_i][_qp] * _diffusivity_main[_qp] * _grad_phi[_j][_qp] * _normals[_qp];
      break;

    case Moose::ElementNeighbor:
      jac = _test[_i][_qp] * -_diffusivity_neighbor[_qp] * _grad_phi_neighbor[_j][_qp] * _normals[_qp];
      break;
  }

  return jac;

  /*
  Real jac = 0;

  switch (type)
  {
    case Moose::ElementElement:
      jac -= _diffusivity_main[_qp] * _grad_phi[_j][_qp] * _normals[_qp] * _test[_i][_qp];
      return jac;

    case Moose::NeighborNeighbor:
      jac += _diffusivity_neighbor[_qp] * _grad_phi_neighbor[_j][_qp] * _normals[_qp] *
             _test_neighbor[_i][_qp];
      return jac;

    case Moose::NeighborElement:
      jac += _diffusivity_main[_qp] * _grad_phi[_j][_qp] * _normals[_qp] * _test_neighbor[_i][_qp];
      return jac;

    case Moose::ElementNeighbor:
      jac -=
          _diffusivity_neighbor[_qp] * _grad_phi_neighbor[_j][_qp] * _normals[_qp] * _test[_i][_qp];
      return jac;
  }

  mooseError("Internal Error. (InterfaceFluxConservation)");
  */
}

