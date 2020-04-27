//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComsolElectronBC.h"

// MOOSE includes
#include "MooseVariable.h"

registerMooseObject("ZapdosApp", ComsolElectronBC);

template <>
InputParameters
validParams<ComsolElectronBC>()
{
  InputParameters params = validParams<IntegratedBC>();
  params.addRequiredParam<Real>("r", "The reflection coefficient");
  params.addRequiredCoupledVar("potential", "The electric potential");
  params.addRequiredCoupledVar("mean_en", "The mean energy.");
  params.addRequiredCoupledVar("ip", "The ion density.");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  params.addParam<Real>("se_coeff", 0.01, "The secondary electron emission coefficient. Default: 0.01.");
  return params;
}

ComsolElectronBC::ComsolElectronBC(const InputParameters & parameters)
  : IntegratedBC(parameters),

    _r_units(1. / getParam<Real>("position_units")),
    _r(getParam<Real>("r")),

    // Coupled Variables
    _grad_potential(coupledGradient("potential")),
    _potential_id(coupled("potential")),
    _mean_en(coupledValue("mean_en")),
    _mean_en_id(coupled("mean_en")),

    _muem(getMaterialProperty<Real>("muem")),
    _d_muem_d_actual_mean_en(getMaterialProperty<Real>("d_muem_d_actual_mean_en")),
    _massem(getMaterialProperty<Real>("massem")),
    _e(getMaterialProperty<Real>("e")),
    //_se_coeff(getMaterialProperty<Real>("se_coeff")),
    _se_coeff(getParam<Real>("se_coeff")),
    _a(0.5),
    _v_thermal(0),
    _ion_flux(0, 0, 0),
    _n_gamma(0),
    _d_v_thermal_d_u(0),
    _d_v_thermal_d_mean_en(0),
    _d_ion_flux_d_potential(0, 0, 0),
    _d_ion_flux_d_ip(0, 0, 0),
    _d_n_gamma_d_potential(0),
    _d_n_gamma_d_ip(0),
    _d_n_gamma_d_u(0),
    _d_n_gamma_d_mean_en(0),
    _actual_mean_en(0)
{
  _num_ions = coupledComponents("ip");

  // Resize the vectors to store _num_ions values:
  _ip.resize(_num_ions);
  _grad_ip.resize(_num_ions);
  _ip_var.resize(_num_ions);
  _ion_id.resize(_num_ions);
  _muip.resize(_num_ions);
  _Dip.resize(_num_ions);
  _sgnip.resize(_num_ions);

  // Retrieve the values for each ion and store in the relevant vectors.
  // Note that these need to be dereferenced to get the values inside the
  // main body of the code.
  // e.g. instead of "_ip[_qp]" it would be "(*_ip[i])[_qp]", where "i"
  // refers to a single ion species.
  for (unsigned int i = 0; i < _ip.size(); ++i)
  {
    _ip_var[i] = getVar("ip", i);
    _ip[i] = &coupledValue("ip", i);
    _grad_ip[i] = &coupledGradient("ip", i);
    _ion_id[i] = _ip_var[i]->number();
    _muip[i] = &getMaterialProperty<Real>("mu" + (*getVar("ip", i)).name());
    _Dip[i] = &getMaterialProperty<Real>("diff" + (*getVar("ip", i)).name());
    _sgnip[i] = &getMaterialProperty<Real>("sgn" + (*getVar("ip", i)).name());
  }
}

Real
ComsolElectronBC::computeQpResidual()
{
  if (_normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0)
  {
    _a = 1.0;
  }
  else
  {
    _a = 0.0;
  }

  _ion_flux.zero();
  for (unsigned int i = 0; i < _num_ions; ++i)
  {
    _ion_flux += (*_sgnip[i])[_qp] * (*_muip[i])[_qp] * -_grad_potential[_qp] * _r_units *
                     std::exp((*_ip[i])[_qp]) -
                 (*_Dip[i])[_qp] * std::exp((*_ip[i])[_qp]) * (*_grad_ip[i])[_qp] * _r_units;
  }
  /*
  _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] /
             (_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp] +
              std::numeric_limits<double>::epsilon());
              */
  _n_gamma = (1. - _a) * _se_coeff * _ion_flux * _normals[_qp] /
             (_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp] +
              std::numeric_limits<double>::epsilon());
  
  _v_thermal =
      std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_mean_en[_qp] - _u[_qp]) / (M_PI * _massem[_qp]));

  /*
  return _test[_i][_qp] * _r_units * (1. - _r) / (1. + _r) * (-0.5 * _v_thermal * _n_gamma) -
         _test[_i][_qp] * _r_units * 2. / (1. + _r) * (1. - _a) * _se_coeff[_qp] * _ion_flux *
             _normals[_qp];
             */
  return _test[_i][_qp] * _r_units * (1. - _r) / (1. + _r) * (-0.5 * _v_thermal * _n_gamma) -
         _test[_i][_qp] * _r_units * 2. / (1. + _r) * (1. - _a) * _se_coeff * _ion_flux *
             _normals[_qp];
}

Real
ComsolElectronBC::computeQpJacobian()
{
  if (_normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0)
  {
    _a = 1.0;
  }
  else
  {
    _a = 0.0;
  }

  _ion_flux.zero();
  for (unsigned int i = 0; i < _num_ions; ++i)
  {
    _ion_flux += (*_sgnip[i])[_qp] * (*_muip[i])[_qp] * -_grad_potential[_qp] * _r_units *
                     std::exp((*_ip[i])[_qp]) -
                 (*_Dip[i])[_qp] * std::exp((*_ip[i])[_qp]) * (*_grad_ip[i])[_qp] * _r_units;
  }

  /*
  _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] /
             (_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp] +
              std::numeric_limits<double>::epsilon());
              */
  _n_gamma = (1. - _a) * _se_coeff * _ion_flux * _normals[_qp] /
             (_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp] +
              std::numeric_limits<double>::epsilon());

  _actual_mean_en = std::exp(_mean_en[_qp] - _u[_qp]);
  /*
  _d_n_gamma_d_u = -1. * (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] /
                   (std::pow(_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp], 2.) +
                    std::numeric_limits<double>::epsilon()) *
                   -_grad_potential[_qp] * _r_units * _normals[_qp] *
                   _d_muem_d_actual_mean_en[_qp] * _actual_mean_en * -_phi[_j][_qp];
                   */
  _d_n_gamma_d_u = -1. * (1. - _a) * _se_coeff * _ion_flux * _normals[_qp] /
                   (std::pow(_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp], 2.) +
                    std::numeric_limits<double>::epsilon()) *
                   -_grad_potential[_qp] * _r_units * _normals[_qp] *
                   _d_muem_d_actual_mean_en[_qp] * _actual_mean_en * -_phi[_j][_qp];
  _v_thermal =
      std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_mean_en[_qp] - _u[_qp]) / (M_PI * _massem[_qp]));
  _d_v_thermal_d_u = 0.5 / _v_thermal * 8 * _e[_qp] * 2.0 / 3 * std::exp(_mean_en[_qp] - _u[_qp]) /
                     (M_PI * _massem[_qp]) * -_phi[_j][_qp];

  return _test[_i][_qp] * _r_units * (1. - _r) / (1. + _r) * -0.5 *
         (_d_v_thermal_d_u * _n_gamma + _v_thermal * _d_n_gamma_d_u);
  // return 0.;
}

Real
ComsolElectronBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  _iter = std::find(_ion_id.begin(), _ion_id.end(), jvar);
  if (jvar == _potential_id)
  {
    if (_normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0)
      _a = 1.0;
    else
      _a = 0.0;

    _ion_flux.zero();
    _d_ion_flux_d_potential.zero();
    for (unsigned int i = 0; i < _num_ions; ++i)
    {
      _ion_flux += (*_sgnip[i])[_qp] * (*_muip[i])[_qp] * -_grad_potential[_qp] * _r_units *
                       std::exp((*_ip[i])[_qp]) -
                   (*_Dip[i])[_qp] * std::exp((*_ip[i])[_qp]) * (*_grad_ip[i])[_qp] * _r_units;

      _d_ion_flux_d_potential += (*_sgnip[i])[_qp] * (*_muip[i])[_qp] * -_grad_phi[_j][_qp] *
                                 _r_units * std::exp((*_ip[i])[_qp]);
    }

    /*
    _d_n_gamma_d_potential = (1. - _a) * _se_coeff[_qp] / _muem[_qp] *
                             (_d_ion_flux_d_potential * _normals[_qp] /
                                  (-_grad_potential[_qp] * _r_units * _normals[_qp] +
                                   std::numeric_limits<double>::epsilon()) -
                              _ion_flux * _normals[_qp] /
                                  (std::pow(-_grad_potential[_qp] * _r_units * _normals[_qp], 2.) +
                                   std::numeric_limits<double>::epsilon()) *
                                  -_grad_phi[_j][_qp] * _r_units * _normals[_qp]);
                                  */
    _d_n_gamma_d_potential = (1. - _a) * _se_coeff / _muem[_qp] *
                             (_d_ion_flux_d_potential * _normals[_qp] /
                                  (-_grad_potential[_qp] * _r_units * _normals[_qp] +
                                   std::numeric_limits<double>::epsilon()) -
                              _ion_flux * _normals[_qp] /
                                  (std::pow(-_grad_potential[_qp] * _r_units * _normals[_qp], 2.) +
                                   std::numeric_limits<double>::epsilon()) *
                                  -_grad_phi[_j][_qp] * _r_units * _normals[_qp]);

    _v_thermal = std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_mean_en[_qp] - _u[_qp]) /
                           (M_PI * _massem[_qp]));

    /*
    return _test[_i][_qp] * _r_units * (1. / (1. + _r)) *
           ((1. - _r) * (-0.5 * _v_thermal * _d_n_gamma_d_potential) -
            2. * (1. - _a) * _se_coeff[_qp] * _d_ion_flux_d_potential * _normals[_qp]);
            */
    return _test[_i][_qp] * _r_units * (1. / (1. + _r)) *
           ((1. - _r) * (-0.5 * _v_thermal * _d_n_gamma_d_potential) -
            2. * (1. - _a) * _se_coeff * _d_ion_flux_d_potential * _normals[_qp]);
  }

  else if (jvar == _mean_en_id)
  {
    if (_normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0)
    {
      _a = 1.0;
    }
    else
    {
      _a = 0.0;
    }

    _ion_flux.zero();
    for (unsigned int i = 0; i < _num_ions; ++i)
    {
      _ion_flux += (*_sgnip[i])[_qp] * (*_muip[i])[_qp] * -_grad_potential[_qp] * _r_units *
                       std::exp((*_ip[i])[_qp]) -
                   (*_Dip[i])[_qp] * std::exp((*_ip[i])[_qp]) * (*_grad_ip[i])[_qp] * _r_units;
    }

    _v_thermal = std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_mean_en[_qp] - _u[_qp]) /
                           (M_PI * _massem[_qp]));

    _d_v_thermal_d_mean_en = 0.5 / _v_thermal * 8 * _e[_qp] * 2.0 / 3 *
                             std::exp(_mean_en[_qp] - _u[_qp]) / (M_PI * _massem[_qp]) *
                             _phi[_j][_qp];

    /*
    _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] /
               (_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp] +
                std::numeric_limits<double>::epsilon());
                */
    _n_gamma = (1. - _a) * _se_coeff * _ion_flux * _normals[_qp] /
               (_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp] +
                std::numeric_limits<double>::epsilon());

    _actual_mean_en = std::exp(_mean_en[_qp] - _u[_qp]);

    /*
    _d_n_gamma_d_mean_en =
        -1. * (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] /
        (std::pow(_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp], 2.) +
         std::numeric_limits<double>::epsilon()) *
        -_grad_potential[_qp] * _r_units * _normals[_qp] * _d_muem_d_actual_mean_en[_qp] *
        _actual_mean_en * _phi[_j][_qp];
        */
    _d_n_gamma_d_mean_en =
        -1. * (1. - _a) * _se_coeff * _ion_flux * _normals[_qp] /
        (std::pow(_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp], 2.) +
         std::numeric_limits<double>::epsilon()) *
        -_grad_potential[_qp] * _r_units * _normals[_qp] * _d_muem_d_actual_mean_en[_qp] *
        _actual_mean_en * _phi[_j][_qp];

    return _test[_i][_qp] * _r_units * (1. - _r) / (1. + _r) * -0.5 *
           (_d_v_thermal_d_mean_en * _n_gamma + _v_thermal * _d_n_gamma_d_mean_en);
  }
  else if (_iter != _ion_id.end())
  {
    _ip_index = std::distance(_ion_id.begin(), _iter);
    if (_normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0)
    {
      _a = 1.0;
    }
    else
    {
      _a = 0.0;
    }

    _v_thermal = std::sqrt(8 * _e[_qp] * 2.0 / 3 * std::exp(_mean_en[_qp] - _u[_qp]) /
                           (M_PI * _massem[_qp]));

    _d_ion_flux_d_ip =
        std::exp((*_ip[_ip_index])[_qp]) *
        ((*_sgnip[_ip_index])[_qp] * (*_muip[_ip_index])[_qp] * -_grad_potential[_qp] * _r_units *
             _phi[_j][_qp] -
         (*_Dip[_ip_index])[_qp] * _grad_phi[_j][_qp] * _r_units -
         (*_Dip[_ip_index])[_qp] * _phi[_j][_qp] * (*_grad_ip[_ip_index])[_qp] * _r_units);

    /*
    _d_n_gamma_d_ip = (1. - _a) * _se_coeff[_qp] * _d_ion_flux_d_ip * _normals[_qp] /
                      (_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp] +
                       std::numeric_limits<double>::epsilon());
                       */
    _d_n_gamma_d_ip = (1. - _a) * _se_coeff * _d_ion_flux_d_ip * _normals[_qp] /
                      (_muem[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp] +
                       std::numeric_limits<double>::epsilon());

    /*
    return _test[_i][_qp] * _r_units * (1. - _r) / (1. + _r) *
               (-0.5 * _v_thermal * _d_n_gamma_d_ip) -
           _test[_i][_qp] * _r_units * 2. / (1. + _r) * (1. - _a) * _se_coeff[_qp] *
               _d_ion_flux_d_ip * _normals[_qp];
               */
    return _test[_i][_qp] * _r_units * (1. - _r) / (1. + _r) *
               (-0.5 * _v_thermal * _d_n_gamma_d_ip) -
           _test[_i][_qp] * _r_units * 2. / (1. + _r) * (1. - _a) * _se_coeff *
               _d_ion_flux_d_ip * _normals[_qp];
  }

  else
    return 0.0;
}
