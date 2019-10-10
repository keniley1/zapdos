//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "DriftDiffusionDoNothingElectronBC.h"

// MOOSE includes
#include "MooseVariable.h"

registerMooseObject("ZapdosApp", DriftDiffusionDoNothingElectronBC);

template <>
InputParameters
validParams<DriftDiffusionDoNothingElectronBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredCoupledVar("mean_en",
                               "The mean electron energy. Used for jacobian calculations.");
  params.addCoupledVar(
      "potential", "The gradient of the potential will be used to compute the advection velocity.");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  params.addParam<Real>("EField",
                        "Optionally can use a specified electric field for 1D "
                        "simulations in place of a potential variable");
  params.addParam<bool>("use_material_props", true, "Whether to use a material for properties.");
  return params;
}

DriftDiffusionDoNothingElectronBC::DriftDiffusionDoNothingElectronBC(
    const InputParameters & parameters)
  : IntegratedBC(parameters),

    _r_units(1. / getParam<Real>("position_units")),

    _muem(getMaterialProperty<Real>("muem")),
    _d_muem_d_actual_mean_en(getMaterialProperty<Real>("d_muem_d_actual_mean_en")),
    _sign(getMaterialProperty<Real>("sgnem")),
    _diffem(getMaterialProperty<Real>("diffem")),
    _d_diffem_d_actual_mean_en(getMaterialProperty<Real>("d_diffem_d_actual_mean_en")),

    _potential_id(coupled("potential")),
    _grad_potential(coupledGradient("potential")),
    _mean_en(coupledValue("mean_en")),
    _mean_en_id(coupled("mean_en")),

    _d_actual_mean_en_d_mean_en(0),
    _d_muem_d_mean_en(0),
    _d_actual_mean_en_d_u(0),
    _d_muem_d_u(0),
    _d_diffem_d_u(0),
    _d_diffem_d_mean_en(0)
{
  if (!(isCoupled("potential") || parameters.isParamSetByUser("EField")))
    mooseError("You must either couple in a potential variable or set an EField.");
}

DriftDiffusionDoNothingElectronBC::~DriftDiffusionDoNothingElectronBC() {}

Real
DriftDiffusionDoNothingElectronBC::computeQpResidual()
{
  return (_muem[_qp] * _sign[_qp] * std::exp(_u[_qp]) * -_grad_potential[_qp] * _r_units *
               -_normals[_qp] -
           _diffem[_qp] * std::exp(_u[_qp]) * _grad_u[_qp] * _r_units * -_normals[_qp]) *
         _test[_i][_qp] * _r_units;
}

Real
DriftDiffusionDoNothingElectronBC::computeQpJacobian()
{
  _d_actual_mean_en_d_u = std::exp(_mean_en[_qp] - _u[_qp]) * -_phi[_j][_qp];
  _d_muem_d_u = _d_muem_d_actual_mean_en[_qp] * _d_actual_mean_en_d_u;
  _d_diffem_d_u =
      _d_diffem_d_actual_mean_en[_qp] * std::exp(_mean_en[_qp] - _u[_qp]) * -_phi[_j][_qp];

  return ((_d_muem_d_u * _sign[_qp] * std::exp(_u[_qp]) * -_grad_potential[_qp] * _r_units *
                -_normals[_qp] +
            _muem[_qp] * _sign[_qp] * std::exp(_u[_qp]) * _phi[_j][_qp] * -_grad_potential[_qp] *
                -_normals[_qp] * _r_units) *
               _test[_i][_qp] * _r_units -
           _diffem[_qp] *
               (std::exp(_u[_qp]) * _grad_phi[_j][_qp] * -_normals[_qp] * _r_units +
                std::exp(_u[_qp]) * _phi[_j][_qp] * _grad_u[_qp] * _r_units * -_normals[_qp]) *
               _test[_i][_qp] * _r_units -
           _d_diffem_d_u * std::exp(_u[_qp]) * _grad_u[_qp] * _r_units * -_normals[_qp] *
               _test[_i][_qp] * _r_units);
}

Real
DriftDiffusionDoNothingElectronBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _potential_id)
  {
    return -_muem[_qp] * _sign[_qp] * std::exp(_u[_qp]) * -_grad_phi[_j][_qp] * _r_units *
           -_normals[_qp] * _test[_i][_qp] * _r_units;
  }
  else if (jvar == _mean_en_id)
  {
    _d_actual_mean_en_d_mean_en = std::exp(_mean_en[_qp] - _u[_qp]) * _phi[_j][_qp];
    _d_muem_d_mean_en = _d_muem_d_actual_mean_en[_qp] * _d_actual_mean_en_d_mean_en;
    _d_diffem_d_mean_en =
        _d_diffem_d_actual_mean_en[_qp] * std::exp(_mean_en[_qp] - _u[_qp]) * _phi[_j][_qp];

    return (_d_muem_d_mean_en * _sign[_qp] * std::exp(_u[_qp]) * -_grad_potential[_qp] * _r_units *
                 -_normals[_qp] * _test[_i][_qp] * _r_units -
             _d_diffem_d_mean_en * std::exp(_u[_qp]) * _grad_u[_qp] * _r_units * -_normals[_qp] *
                 _test[_i][_qp] * _r_units);
  }
  else
  {
    return 0.;
  }
}
