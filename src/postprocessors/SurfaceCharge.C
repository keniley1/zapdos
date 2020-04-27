//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SurfaceCharge.h"

// MOOSE includes
#include "MooseVariable.h"

registerMooseObject("ZapdosApp", SurfaceCharge);

template <>
InputParameters
validParams<SurfaceCharge>()
{
  InputParameters params = validParams<SideIntegralVariablePostprocessor>();
  // params.addRequiredParam<std::string>("diffusivity", "The name of the diffusivity material
  // property that will be used in the flux computation.");
  params.addRequiredCoupledVar("potential", "The potential that drives the advective flux.");
  params.addParam<Real>("r_ion", 0.0, "The ion reflection coefficient on this boundary.");
  params.addParam<Real>("r_electron", 0.0, "The electron reflection coefficient on this boundary.");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  params.addRequiredCoupledVar("mean_en", "Electron energy.");
  params.addRequiredCoupledVar("ions", "All of the ions that can interact with this boundary.");
  params.addParam<Real>("se_coeff", "The secondary electron emission coefficient.");
  return params;
}

SurfaceCharge::SurfaceCharge(const InputParameters & parameters)
  : SideIntegralVariablePostprocessor(parameters),
    _muem(getMaterialProperty<Real>("muem")),
    //_diffem(getMaterialProperty<Real>("diffem")),
    _r_units(1. / getParam<Real>("position_units")),
    _r_ion(getParam<Real>("r_ion")),
    _r_electron(getParam<Real>("r_electron")),
    _kb(getMaterialProperty<Real>("k_boltz")),
    _e(getMaterialProperty<Real>("e")),
    _se_coeff(isParamValid("se_coeff") ? getParam<Real>("se_coeff") : 0),
    _a(0.5),
    _grad_potential(coupledGradient("potential")),
    _mean_en(coupledValue("mean_en"))
{
  _num_ions = coupledComponents("ions");

  // Resize the vectors to store _num_ions pointers
  _ions.resize(_num_ions);
  _mu_ions.resize(_num_ions);
  //_diff_ions.resize(_num_ions);
  _sgn_ions.resize(_num_ions);
  _mass_ions.resize(_num_ions);
  _T_ions.resize(_num_ions);
  _grad_ions.resize(_num_ions);

  for (unsigned int i = 0; i < _ions.size(); ++i)
  {
    _ions[i] = &coupledValue("ions", i);
    _grad_ions[i] = &coupledGradient("ions", i);
    _mu_ions[i] = &getMaterialProperty<Real>("mu" + (*getVar("ions", i)).name());
    //_diff_ions[i] = &getMaterialProperty<Real>("diff" + (*getVar("ions", i)).name());
    _sgn_ions[i] = &getMaterialProperty<Real>("sgn" + (*getVar("ions", i)).name());
    _mass_ions[i] = &getMaterialProperty<Real>("mass" + (*getVar("ions", i)).name());
    _T_ions[i] = &getMaterialProperty<Real>("T" + (*getVar("ions", i)).name());
  }
}

Real
SurfaceCharge::computeQpIntegral()
{
  // Output units for base case are: mol / (m^2 * s)

  if (_normals[_qp] * -1 * -_grad_potential[_qp] > 0.0)
  {
    _b = 1.0;
  }
  else
  {
    _b = 0.0;
  }

  _ve_thermal =
      std::sqrt(8 * 1.602e-19 * 2.0 / 3 * std::exp(_mean_en[_qp] - _u[_qp]) / (M_PI * 9.11e-31));

  Real electron_flux;

  electron_flux = (1.0 - _r_electron) / (1.0 + _r_electron) *
                  (-(2 * _b - 1) * _muem[_qp] * -_grad_potential[_qp] * _r_units *
                       std::exp(_u[_qp]) * _normals[_qp] +
                   0.5 * _ve_thermal * std::exp(_u[_qp]));

  _ion_flux = 0.0;
  for (unsigned int i = 0; i < _num_ions; ++i)
  {
    if (_normals[_qp] * (*_sgn_ions[i])[_qp] * -_grad_potential[_qp] > 0.0)
    {
      _a = 1.0;
    }
    else
    {
      _a = 0.0;
    }
    _ion_flux +=
        std::exp((*_ions[i])[_qp]) *
        (0.5 * std::sqrt(8 * _kb[_qp] * (*_T_ions[i])[_qp] / (M_PI * (*_mass_ions[i])[_qp])) +
         (2 * _a - 1) * (*_sgn_ions[i])[_qp] * (*_mu_ions[i])[_qp] * -_grad_potential[_qp] *
             _r_units * _normals[_qp]);
  }
  _ion_flux *= (1. - _r_ion) / (1. + _r_ion);

  _n_gamma = (1. - _b) * _se_coeff * _ion_flux /
             (_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp] +
              std::numeric_limits<double>::epsilon());

  // Add secondary electrons
  electron_flux += ((1. - _r_electron) / (1. + _r_electron) * (-0.5 * _ve_thermal * _n_gamma) -
                    (2.0 / (1 + _r_electron) * (1. - _b) * _se_coeff * _ion_flux));

  /*
  return _test[_i][_qp] * _r_units * (1. - _r) / (1. + _r) * (-0.5 * _ve_thermal * _n_gamma) -
         _test[_i][_qp] * _r_units * 2. / (1. + _r) * (1. - _a) * _se_coeff[_qp] * _ion_flux *
             _normals[_qp];
             */

  return (_ion_flux - electron_flux) * 6.022e23 * 1.602e-19;
}
