//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HagelaarEnergyBC.h"

// MOOSE includes
#include "MooseVariable.h"

registerMooseObject("ZapdosApp", HagelaarEnergyBC);

template <>
InputParameters
validParams<HagelaarEnergyBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<Real>("r", "The reflection coefficient");
  params.addRequiredCoupledVar("potential", "The electric potential");
  params.addRequiredCoupledVar("em", "The electron density.");
  params.addRequiredCoupledVar("ip", "The ion density.");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  return params;
}

HagelaarEnergyBC::HagelaarEnergyBC(const InputParameters & parameters)
  : IntegratedBC(parameters),

    _r_units(1. / getParam<Real>("position_units")),
    _r(getParam<Real>("r")),

    // Coupled Variables
    _grad_potential(coupledGradient("potential")),
    _potential_id(coupled("potential")),
    _em(coupledValue("em")),
    _em_id(coupled("em")),
    _ip_var(*getVar("ip", 0)),
    //_ip(coupledValue("ip")),
    //_grad_ip(coupledGradient("ip")),
    //_ip_id(coupled("ip")),

    _muem(getMaterialProperty<Real>("muem")),
    _d_muem_d_actual_mean_en(getMaterialProperty<Real>("d_muem_d_actual_mean_en")),
    _massem(getMaterialProperty<Real>("massem")),
    _e(getMaterialProperty<Real>("e")),
    _sgnip(getMaterialProperty<Real>("sgn" + _ip_var.name())),
    _muip(getMaterialProperty<Real>("mu" + _ip_var.name())),
    _Dip(getMaterialProperty<Real>("diff" + _ip_var.name())),
    _se_coeff(getMaterialProperty<Real>("se_coeff")),
    _se_energy(getMaterialProperty<Real>("se_energy")),
    _mumean_en(getMaterialProperty<Real>("mumean_en")),
    _d_mumean_en_d_actual_mean_en(getMaterialProperty<Real>("d_mumean_en_d_actual_mean_en")),
    _a(0.5),
    _v_thermal(0),
    _ion_flux(0, 0, 0),
    _n_gamma(0),
    _d_v_thermal_d_u(0),
    _d_v_thermal_d_em(0),
    _d_ion_flux_d_potential(0, 0, 0),
    _d_ion_flux_d_ip(0, 0, 0),
    _d_n_gamma_d_potential(0),
    _d_n_gamma_d_ip(0),
    _d_n_gamma_d_u(0),
    _d_n_gamma_d_em(0),
    _actual_mean_en(0)
{
  int n = coupledComponents("ip");

  //_ip_var.resize(n);
  _ip.resize(n);
  _grad_ip.resize(n);
  _ip_id.resize(n);
  //_sgnip.resize(n);
  //_muip.resize(n);
  //_Dip.resize(n);
  for (unsigned int i = 0; i < _ip.size(); ++i)
  {
    //_ip_var[i] = getVar("ip", i);
    _ip[i] = &coupledValue("ip", i);
    _grad_ip[i] = &coupledGradient("ip", i);
    _ip_id[i] = coupled("ip", i);
    //_sgnip[i] = &getMaterialProperty<Real>("sgn" + _ip[i]
  }
}

Real
HagelaarEnergyBC::computeQpResidual()
{
  if (_normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0)
  {
    _a = 1.0;
  }
  else
  {
    _a = 0.0;
  }
  // Real total_ion_flux;
  _ion_flux = 0.0;
  for (unsigned int i = 0; i < _ip.size(); ++i)
  {
    _ion_flux +=
        _sgnip[_qp] * _muip[_qp] * -_grad_potential[_qp] * _r_units * std::exp((*_ip[i])[_qp]) -
        _Dip[_qp] * std::exp((*_ip[i])[_qp]) * (*_grad_ip[i])[_qp] * _r_units;
  }

  //_ion_flux = _sgnip[_qp] * _muip[_qp] * -_grad_potential[_qp] * _r_units * std::exp(_ip[_qp]) -
  //            _Dip[_qp] * std::exp(_ip[_qp]) * _grad_ip[_qp] * _r_units;
  // _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (_muem[_qp] *
  // -_grad_potential[_qp] * _r_units * _normals[_qp] + std::numeric_limits<double>::epsilon());
  _v_thermal =
      std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_u[_qp] - _em[_qp]) / (M_PI * _massem[_qp]));

  return _test[_i][_qp] * _r_units / (6. * (_r + 1.)) *
         (10. * _ion_flux * _normals[_qp] * _se_energy[_qp] * _se_coeff[_qp] * (_a - 1.) *
              (_r + 1.) +
          (_r - 1.) * (std::exp(_u[_qp]) - _se_energy[_qp] * _n_gamma) *
              (6. * -_grad_potential[_qp] * _r_units * _normals[_qp] * _mumean_en[_qp] *
                   (2. * _a - 1.) -
               5. * _v_thermal));
}

Real
HagelaarEnergyBC::computeQpJacobian()
{
  if (_normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0)
  {
    _a = 1.0;
  }
  else
  {
    _a = 0.0;
  }

  _ion_flux = 0.0;
  for (unsigned int i = 0; i < _ip.size(); ++i)
  {
    _ion_flux +=
        _sgnip[_qp] * _muip[_qp] * -_grad_potential[_qp] * _r_units * std::exp((*_ip[i])[_qp]) -
        _Dip[_qp] * std::exp((*_ip[i])[_qp]) * (*_grad_ip[i])[_qp] * _r_units;
  }
  //_ion_flux = _sgnip[_qp] * _muip[_qp] * -_grad_potential[_qp] * _r_units * std::exp(_ip[_qp]) -
  //            _Dip[_qp] * std::exp(_ip[_qp]) * _grad_ip[_qp] * _r_units;
  // _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (_muem[_qp] *
  // -_grad_potential[_qp] * _r_units * _normals[_qp] + std::numeric_limits<double>::epsilon());
  _actual_mean_en = std::exp(_u[_qp] - _em[_qp]);
  // _d_n_gamma_d_u = -1. * (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] /
  // (std::pow(_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp], 2.) +
  // std::numeric_limits<double>::epsilon()) * -_grad_potential[_qp] * _r_units * _normals[_qp] *
  // _d_muem_d_actual_mean_en[_qp] * _actual_mean_en * _phi[_j][_qp];
  _v_thermal =
      std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_u[_qp] - _em[_qp]) / (M_PI * _massem[_qp]));
  _d_v_thermal_d_u = 0.5 / _v_thermal * 8 * _e[_qp] * 2.0 / 3 * std::exp(_u[_qp] - _em[_qp]) /
                     (M_PI * _massem[_qp]) * _phi[_j][_qp];

  return _test[_i][_qp] * _r_units / (6. * (_r + 1.)) * (_r - 1.) *
         ((std::exp(_u[_qp]) * _phi[_j][_qp] - _se_energy[_qp] * _d_n_gamma_d_u) *
              (6. * -_grad_potential[_qp] * _r_units * _normals[_qp] * _mumean_en[_qp] *
                   (2. * _a - 1.) -
               5. * _v_thermal) +
          (std::exp(_u[_qp]) - _se_energy[_qp] * _n_gamma) *
              (6. * -_grad_potential[_qp] * _r_units * _normals[_qp] *
                   _d_mumean_en_d_actual_mean_en[_qp] * _actual_mean_en * _phi[_j][_qp] *
                   (2. * _a - 1.) -
               5. * _d_v_thermal_d_u));
}

Real
HagelaarEnergyBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _potential_id)
  {
    if (_normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0)
      _a = 1.0;
    else
      _a = 0.0;

    _ion_flux = 0.0;
    _d_ion_flux_d_potential = 0.0;
    for (unsigned int i = 0; i < _ip.size(); ++i)
    {
      _ion_flux +=
          _sgnip[_qp] * _muip[_qp] * -_grad_potential[_qp] * _r_units * std::exp((*_ip[i])[_qp]) -
          _Dip[_qp] * std::exp((*_ip[i])[_qp]) * (*_grad_ip[i])[_qp] * _r_units;
      _d_ion_flux_d_potential +=
          _sgnip[_qp] * _muip[_qp] * -_grad_phi[_j][_qp] * _r_units * std::exp((*_ip[i])[_qp]);
    }
    //_ion_flux = _sgnip[_qp] * _muip[_qp] * -_grad_potential[_qp] * _r_units * std::exp(_ip[_qp]) -
    //            _Dip[_qp] * std::exp(_ip[_qp]) * _grad_ip[_qp] * _r_units;
    //_d_ion_flux_d_potential =
    //    _sgnip[_qp] * _muip[_qp] * -_grad_phi[_j][_qp] * _r_units * std::exp(_ip[_qp]);
    // _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (_muem[_qp] *
    // -_grad_potential[_qp] * _r_units * _normals[_qp] + std::numeric_limits<double>::epsilon());
    // _d_n_gamma_d_potential = (1. - _a) * _se_coeff[_qp] / _muem[_qp] * (_d_ion_flux_d_potential *
    // _normals[_qp] / (-_grad_potential[_qp] * _r_units * _normals[_qp] +
    // std::numeric_limits<double>::epsilon()) - _ion_flux * _normals[_qp] /
    // (std::pow(-_grad_potential[_qp] * _r_units * _normals[_qp], 2.) +
    // std::numeric_limits<double>::epsilon()) * -_grad_phi[_j][_qp] * _r_units * _normals[_qp]);
    _v_thermal =
        std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_u[_qp] - _em[_qp]) / (M_PI * _massem[_qp]));

    return _test[_i][_qp] * _r_units / (6. * (_r + 1.)) *
           (10. * _d_ion_flux_d_potential * _normals[_qp] * _se_energy[_qp] * _se_coeff[_qp] *
                (_a - 1.) * (_r + 1.) +
            (_r - 1.) * ((-_se_energy[_qp] * _d_n_gamma_d_potential) *
                             (6. * -_grad_potential[_qp] * _r_units * _normals[_qp] *
                                  _mumean_en[_qp] * (2. * _a - 1.) -
                              5. * _v_thermal) +
                         (std::exp(_u[_qp]) - _se_energy[_qp] * _n_gamma) *
                             (6. * -_grad_phi[_j][_qp] * _r_units * _normals[_qp] *
                              _mumean_en[_qp] * (2. * _a - 1.))));
  }

  else if (jvar == _em_id)
  {
    if (_normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0)
    {
      _a = 1.0;
    }
    else
    {
      _a = 0.0;
    }
    _v_thermal =
        std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_u[_qp] - _em[_qp]) / (M_PI * _massem[_qp]));
    _d_v_thermal_d_em = 0.5 / _v_thermal * 8 * _e[_qp] * 2.0 / 3 * std::exp(_u[_qp] - _em[_qp]) /
                        (M_PI * _massem[_qp]) * -_phi[_j][_qp];
    // _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] / (_muem[_qp] *
    // -_grad_potential[_qp] * _r_units * _normals[_qp] + std::numeric_limits<double>::epsilon());
    _actual_mean_en = std::exp(_u[_qp] - _em[_qp]);
    // _d_n_gamma_d_em = -1. * (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] /
    // (std::pow(_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp], 2.) +
    // std::numeric_limits<double>::epsilon()) * -_grad_potential[_qp] * _r_units * _normals[_qp] *
    // _d_muem_d_actual_mean_en[_qp] * _actual_mean_en * -_phi[_j][_qp];

    return _test[_i][_qp] * _r_units / (6. * (_r + 1.)) *
           ((_r - 1.) * (std::exp(_u[_qp]) - _se_energy[_qp] * _n_gamma) *
                (6. * -_grad_potential[_qp] * _r_units * _normals[_qp] *
                     _d_mumean_en_d_actual_mean_en[_qp] * _actual_mean_en * -_phi[_j][_qp] *
                     (2. * _a - 1.) -
                 5. * _d_v_thermal_d_em) +
            (_r - 1.) *
                (6. * -_grad_potential[_qp] * _r_units * _normals[_qp] * _mumean_en[_qp] *
                     (2. * _a - 1.) -
                 5. * _v_thermal) *
                -_se_energy[_qp] * _d_n_gamma_d_em);
  }

  // else if (jvar == _ip_id)
  else
  {
    // First find if jvar refers to one of the ions by checking the _ip_id numbers.
    _it = std::find(_ip_id.begin(), _ip_id.end(), jvar);
    if (_it != _ip_id.end())
    {
      // If it does, store the index of that ion (_ip_num) and use to compute the relevant Jacobian
      // terms.
      _ip_num = std::distance(_ip_id.begin(), _it);
      if (_normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0)
      {
        _a = 1.0;
      }
      else
      {
        _a = 0.0;
      }
      _v_thermal =
          std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_u[_qp] - _em[_qp]) / (M_PI * _massem[_qp]));

      _d_ion_flux_d_ip =
          _sgnip[_qp] * _muip[_qp] * -_grad_potential[_qp] * _r_units *
              std::exp((*_ip[_ip_num])[_qp]) * _phi[_j][_qp] -
          _Dip[_qp] * std::exp((*_ip[_ip_num])[_qp]) * _grad_phi[_j][_qp] * _r_units -
          _Dip[_qp] * std::exp((*_ip[_ip_num])[_qp]) * _phi[_j][_qp] * (*_grad_ip[_ip_num])[_qp] * _r_units;
      //_d_ion_flux_d_ip = _sgnip[_qp] * _muip[_qp] * -_grad_potential[_qp] * _r_units *
      //                       std::exp(_ip[_qp]) * _phi[_j][_qp] -
      //                   _Dip[_qp] * std::exp(_ip[_qp]) * _grad_phi[_j][_qp] * _r_units -
      //                   _Dip[_qp] * std::exp(_ip[_qp]) * _phi[_j][_qp] * _grad_ip[_qp] *
      //                   _r_units;
      // _d_n_gamma_d_ip = (1. - _a) * _se_coeff[_qp] * _d_ion_flux_d_ip * _normals[_qp] /
      // (_muem[_qp]
      // * -_grad_potential[_qp] * _r_units * _normals[_qp] +
      // std::numeric_limits<double>::epsilon());

      return _test[_i][_qp] * _r_units / (6. * (_r + 1.)) *
             (10. * _d_ion_flux_d_ip * _normals[_qp] * _se_energy[_qp] * _se_coeff[_qp] *
                  (_a - 1.) * (_r + 1.) +
              (_r - 1.) * (-_se_energy[_qp] * _d_n_gamma_d_ip) *
                  (6. * -_grad_potential[_qp] * _r_units * _normals[_qp] * _mumean_en[_qp] *
                       (2. * _a - 1.) -
                   5. * _v_thermal));
    }

    else
      return 0.0;
  }
}
