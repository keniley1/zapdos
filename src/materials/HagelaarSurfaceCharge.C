//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HagelaarSurfaceCharge.h"

registerMooseObject("ZapdosApp", HagelaarSurfaceCharge);

template <>
InputParameters
validParams<HagelaarSurfaceCharge>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("potential", "The potential that drives the advective flux.");
  params.addParam<Real>("r_ion", 0.0, "The ion reflection coefficient on this boundary.");
  params.addParam<Real>("r_electron", 0.0, "The electron reflection coefficient on this boundary.");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  params.addRequiredCoupledVar("mean_en", "Electron energy.");
  params.addRequiredCoupledVar("em", "Electron density.");
  params.addRequiredCoupledVar("ions", "All of the ions that can interact with this boundary.");
  params.addParam<Real>("se_coeff", "The secondary electron emission coefficient.");
  params.addRequiredParam<std::string>("potential_units", "The potential units.");
  return params;
}

HagelaarSurfaceCharge::HagelaarSurfaceCharge(const InputParameters & parameters)
  : Material(parameters),

    // Declare material properties
    _sigma(declareProperty<Real>("surface_charge")),
    _sigma_old(getMaterialPropertyOld<Real>("surface_charge")),
    //_d_ion_flux_d_ions(declareProperty<std::vector<Real>>("d_ion_flux_d_ions")),

    // Coupled Variables
    _muem(getMaterialProperty<Real>("muem")),
    //_diffem(getMaterialProperty<Real>("diffem")),
    _r_units(1. / getParam<Real>("position_units")),
    _r_ion(getParam<Real>("r_ion")),
    _r_electron(getParam<Real>("r_electron")),
    _kb(getMaterialProperty<Real>("k_boltz")),
    _e(getMaterialProperty<Real>("e")),
    _se_coeff(isParamValid("se_coeff") ? getParam<Real>("se_coeff") : 0),
    _potential_units(getParam<std::string>("potential_units")),
    _a(0.5),
    _b(0.5),
    _grad_potential(coupledGradient("potential")),
    _mean_en(coupledValue("mean_en")),
    _em(coupledValue("em")),
    _grad_em(coupledGradient("em"))
{
  // Define voltage scaling parameter
  if (_potential_units.compare("V") == 0)
    _voltage_scaling = 1.;
  else if (_potential_units.compare("kV") == 0)
    _voltage_scaling = 1000;

  // Initialize the vectors to store ion density and transport values
  _num_ions = coupledComponents("ions");

  //_d_ion_flux_d_ions.resize(_num_ions);

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

void
HagelaarSurfaceCharge::initQpStatefulProperties()
{
  _sigma[_qp] = 0;
}

void
HagelaarSurfaceCharge::computeQpProperties()
{
  if (_material_data_type == Moose::FACE_MATERIAL_DATA || boundaryRestricted())
  {
    if (_normals[_qp] * -1 * -_grad_potential[_qp] > 0.0)
    {
      _b = 1.0;
    }
    else
    {
      _b = 0.0;
    }

    _ve_thermal =
        std::sqrt(8 * 1.602e-19 * 2.0 / 3 * std::exp(_mean_en[_qp] - _em[_qp]) / (M_PI * 9.11e-31));

    // Real electron_flux;
/*
    _electron_flux =
        std::exp(_em[_qp]) *
        (_b * -_muem[_qp] * -_grad_potential[_qp] * _r_units - _diffem[_qp] * _grad_em[_qp] * _r_units) *
        _normals[_qp];
        */
    _electron_flux = (1.0 - _r_electron) / (1.0 + _r_electron) *
                     (-(2 * _b - 1) * _muem[_qp] * -_grad_potential[_qp] * _r_units *
                          std::exp(_em[_qp]) * _normals[_qp] +
                      0.5 * _ve_thermal * std::exp(_em[_qp]));

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

      /*
      _ion_flux += (*_sgn_ions[i])[_qp] * std::exp((*_ions[i])[_qp]) *
                   (_a * (*_sgn_ions[i])[_qp] * (*_mu_ions[i])[_qp] * -_grad_potential[_qp] * _r_units -
                    (*_diff_ions[i])[_qp] * (*_grad_ions[i])[_qp] * _r_units) *
                   _normals[_qp];
                   */
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

    // Subtract secondary electrons from the electron flux
    _electron_flux += ((1. - _r_electron) / (1. + _r_electron) * (-0.5 * _ve_thermal * _n_gamma) -
                       (2.0 / (1 + _r_electron) * (1. - _b) * _se_coeff * _ion_flux));
    /*
    _electron_flux += (1. - _b) * _se_coeff * _ion_flux;
    */

    _sigma[_qp] =
        _sigma_old[_qp] + (_ion_flux - _electron_flux) * 6.022e23 * 1.602e-19 * _dt / _voltage_scaling;


  }
  else
    _sigma[_qp] = 0.;
}
