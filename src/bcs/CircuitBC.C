//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CircuitBC.h"

// MOOSE includes
#include "MooseVariable.h"
#include "Function.h"

registerADMooseObject("ZapdosApp", CircuitBC);

defineADLegacyParams(CircuitBC);

InputParameters
CircuitBC::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addRequiredParam<FunctionName>("function", "The function.");
  params.addRequiredParam<UserObjectName>(
      "data_provider",
      "The name of the UserObject that can provide some data to materials, bcs, etc.");
  params.addRequiredCoupledVar("ip", "The ion density.");
  params.addRequiredCoupledVar("em", "The log of the electron density.");
  params.addRequiredCoupledVar(
      "mean_en", "The log of the product of the mean energy and the electron density.");
  params.addRequiredParam<std::string>("potential_units", "The potential units.");
  params.addRequiredParam<Real>("r",
                                "The reflection coefficient applied to both electrons and ions");
  params.addRequiredParam<Real>("position_units", "Units of position.");

  return params;
}

CircuitBC::CircuitBC(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
    _r_units(1. / getParam<Real>("position_units")),
    _V_bat(getFunction("function")),
    _data(getUserObject<ProvideMobility>("data_provider")),
    _mean_en(adCoupledValue("mean_en")),
    _em(adCoupledValue("em")),
    _se_coeff(getMaterialProperty<Real>("se_coeff")),
    _eps(getMaterialProperty<Real>("eps")),
    _N_A(getMaterialProperty<Real>("N_A")),
    _muem(getADMaterialProperty<Real>("muem")),
    _e(getMaterialProperty<Real>("e")),
    _massem(getMaterialProperty<Real>("massem")),
    _kb(getMaterialProperty<Real>("k_boltz")),

    _potential_units(getParam<std::string>("potential_units")),
    _r(getParam<Real>("r"))
{
  _ion_flux.zero();
  _n_gamma = 0.0;
  _actual_mean_en = 0.0;
  _v_e_th = 0.0;
  _v_i_th = 0.0;
  _a = 0.0;
  _b = 0.0;
  _numerator = 0.0;
  _denominator = 0.0;
  if (_potential_units.compare("V") == 0)
    _voltage_scaling = 1.;
  else if (_potential_units.compare("kV") == 0)
    _voltage_scaling = 1000;

  // First we need to initialize all of the ion densities and material properties.
  // Find the number of ions coupled into this BC:
  _num_ions = coupledComponents("ip");

  // Resize the vectors to store _num_ions values:
  _ip.resize(_num_ions);
  _grad_ip.resize(_num_ions);
  _T_heavy.resize(_num_ions);
  _muip.resize(_num_ions);
  _Dip.resize(_num_ions);
  _sgnip.resize(_num_ions);
  _mass.resize(_num_ions);

  // Retrieve the values for each ion and store in the relevant vectors.
  // Note that these need to be dereferenced to get the values inside the
  // main body of the code.
  // e.g. instead of "_ip[_qp]" it would be "(*_ip[i])[_qp]", where "i"
  // refers to a single ion species.
  for (unsigned int i = 0; i < _ip.size(); ++i)
  {
    _ip[i] = &adCoupledValue("ip", i);
    _grad_ip[i] = &adCoupledGradient("ip", i);
    _T_heavy[i] = &getADMaterialProperty<Real>("T" + (*getVar("ip", i)).name());
    _muip[i] = &getADMaterialProperty<Real>("mu" + (*getVar("ip", i)).name());
    _Dip[i] = &getADMaterialProperty<Real>("diff" + (*getVar("ip", i)).name());
    _sgnip[i] = &getMaterialProperty<Real>("sgn" + (*getVar("ip", i)).name());
    _mass[i] = &getMaterialProperty<Real>("mass" + (*getVar("ip", i)).name());
  }
}

ADReal
CircuitBC::computeQpResidual()
{
  if (_normals[_qp] * -1.0 * -_grad_u[_qp] > 0.0)
  {
    _a = 1.0;
    _b = 0.;
  }
  else
  {
    _a = 0.0;
    _b = 1.;
  }

  _ion_flux.zero();
  _n_gamma = 0;
  _secondary_ion = 0;
  _ion_drift = 0;

  ADReal _v_i_th = 0;
  for (unsigned int i = 0; i < _num_ions; ++i)
  {
    _ion_flux +=
        (*_sgnip[i])[_qp] * (*_muip[i])[_qp] * -_grad_u[_qp] * _r_units * std::exp((*_ip[i])[_qp]) -
        (*_Dip[i])[_qp] * std::exp((*_ip[i])[_qp]) * (*_grad_ip[i])[_qp] * _r_units;

    _secondary_ion += std::exp((*_ip[i])[_qp]) * (*_muip[i])[_qp];

    _ion_drift += std::sqrt(8 * _kb[_qp] * (*_T_heavy[i])[_qp] / (M_PI * (*_mass[i])[_qp])) *
                  std::exp((*_ip[i])[_qp]);

    _v_i_th += std::sqrt(8 * _kb[_qp] * (*_T_heavy[i])[_qp] / (M_PI * (*_mass[i])[_qp]));
  }
  _n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] /
             (_muem[_qp] * -_grad_u[_qp] * _r_units * _normals[_qp]);

  _v_e_th = std::sqrt(8 * _data.coulomb_charge() * 2.0 / 3 * std::exp(_mean_en[_qp] - _em[_qp]) /
                      (M_PI * _massem[_qp]));

  Real C, D, q;
  C = (1 - _r) / (1 + _r);
  D = 2.0 / (1 + _r);
  q = 1.602e-19;

  /*
  return -_test[_i][_qp] * _r_units * _eps[_qp] *
         ((-1.0 / (C * _data.electrode_area() * _data.ballast_resist() * 1.602e-19) *
               (_u[_qp] + _V_bat.value(_t, _q_point[_qp])) +
           0.5 * _v_i_th * std::exp((*_ip[0])[_qp]) * _N_A[_qp] -
           0.5 * _v_e_th * (std::exp(_em[_qp]) - _n_gamma) * _N_A[_qp] +
           D * (1 - _a) * _se_coeff[_qp] * _v_i_th * std::exp((*_ip[0])[_qp]) * _N_A[_qp]) /
          ((2 * _b - 1) * (*_muip[0])[_qp] * std::exp((*_ip[0])[_qp]) * _N_A[_qp] /
               _voltage_scaling * (1.0 + D * (1.0 - _a)) +
           (2 * _a - 1) * _muem[_qp] * std::exp(_em[_qp])) * _N_A[_qp]/
          _voltage_scaling);
          */

  /*
  return -_test[_i][_qp] * _r_units * _eps[_qp] *
         (-(_u[_qp] + _V_bat.value(_t, _q_point[_qp])) /
              (C * q * _data.electrode_area() * _data.ballast_resist()) +
          (_v_i_th * std::exp((*_ip[0])[_qp]) * (1.0 + D * (1 - _a) * _se_coeff[_qp]) -
           _v_e_th * (std::exp(_em[_qp]) - _n_gamma) * 0.5 * _N_A[_qp]) /
              ((*_muip[0])[_qp] * std::exp((*_ip[0])[_qp]) *
                   ((2 * _b - 1) + D * (1 - _a) * _se_coeff[_qp] * (2 * _b - 1)) +
               (2 * _a - 1) * _muem[_qp] * std::exp(_em[_qp])) / _voltage_scaling * _N_A[_qp]);
               */
  /*
  return -_test[_i][_qp] * _r_units * _eps[_qp] *
         ((-(_u[_qp] - _V_bat.value(_t, _q_point[_qp])) /
               (C * q * _data.electrode_area() * _data.ballast_resist()) +
           (_v_i_th * std::exp((*_ip[0])[_qp]) * (1.0 + D * (1 - _a) * _se_coeff[_qp]) -
            _v_e_th * (std::exp(_em[_qp]) - _n_gamma)) *
               0.5 * _N_A[_qp]) /
          (((*_muip[0])[_qp] * std::exp((*_ip[0])[_qp]) *
                ((2 * _b - 1) + D * (1 - _a) * _se_coeff[_qp] * (2 * _b - 1)) +
            (2 * _a - 1) * _muem[_qp] * std::exp(_em[_qp])) /
           _voltage_scaling * _N_A[_qp]));
           */

  Real num0, num1, den0, den1, test0, test1;

  num0 = MetaPhysicL::raw_value(
      -(_u[_qp] - _V_bat.value(_t, _q_point[_qp])) /
          (C * q * _data.electrode_area() * _data.ballast_resist() / _voltage_scaling) +
      (_v_i_th * std::exp((*_ip[0])[_qp]) * (1.0 + D * (1 - _a) * _se_coeff[_qp]) -
       _v_e_th * (std::exp(_em[_qp]) - _n_gamma)) *
          0.5 * _N_A[_qp]);

  num1 = MetaPhysicL::raw_value((-2. * (1. + _r) * _u[_qp] -
                                 2. * (1. + _r) * -_V_bat.value(_t, _q_point[_qp]) +
                                 _data.electrode_area() * _data.coulomb_charge() *
                                     _data.ballast_resist() / _voltage_scaling * (-1. + _r) *
                                     ((-1. + (-1. + _a) * _se_coeff[_qp]) * _N_A[_qp] * _ion_drift +
                                      _N_A[_qp] * (std::exp(_em[_qp]) - _n_gamma) * _v_e_th)));

  den0 = MetaPhysicL::raw_value(((*_muip[0])[_qp] * std::exp((*_ip[0])[_qp]) *
                                     ((2 * _b - 1) + D * (1 - _a) * _se_coeff[_qp] * (2 * _b - 1)) +
                                 (2 * _a - 1) * _muem[_qp] * std::exp(_em[_qp])) /
                                _voltage_scaling * _N_A[_qp]);

  den1 = MetaPhysicL::raw_value((2. * _data.electrode_area() * _data.coulomb_charge() *
                                 ((-1. + 2. * _a) * _muem[_qp] / _voltage_scaling * _N_A[_qp] *
                                      (std::exp(_em[_qp]) - _n_gamma) -
                                  (-1. + 2. * _b) * (-1. + (-1. + _a) * _se_coeff[_qp]) *
                                      _secondary_ion / _voltage_scaling * _N_A[_qp]) *
                                 _data.ballast_resist() * (-1. + _r)));

  test0 = MetaPhysicL::raw_value(
      (((_u[_qp] - _V_bat.value(_t, _q_point[_qp])) /
            (C * q * _data.electrode_area() * _data.ballast_resist()) +
        (_v_i_th * std::exp((*_ip[0])[_qp]) * (1.0 + D * (1 - _a) * _se_coeff[_qp]) -
         _v_e_th * (std::exp(_em[_qp]) - _n_gamma)) *
            0.5 * _N_A[_qp]) /
       (((*_muip[0])[_qp] * std::exp((*_ip[0])[_qp]) *
             ((2 * _b - 1) + D * (1 - _a) * _se_coeff[_qp] * (2 * _b - 1)) +
         (2 * _a - 1) * _muem[_qp] * std::exp(_em[_qp])) /
        _voltage_scaling * _N_A[_qp])));

  test1 = MetaPhysicL::raw_value(
      (-2. * (1. + _r) * _u[_qp] - 2. * (1. + _r) * -_V_bat.value(_t, _q_point[_qp]) +
       _data.electrode_area() * _data.coulomb_charge() * _data.ballast_resist() / _voltage_scaling *
           (-1. + _r) *
           ((-1. + (-1. + _a) * _se_coeff[_qp]) * _N_A[_qp] * _ion_drift +
            _N_A[_qp] * (std::exp(_em[_qp]) - _n_gamma) * _v_e_th)) /
      (2. * _data.electrode_area() * _data.coulomb_charge() *
       ((-1. + 2. * _a) * _muem[_qp] / _voltage_scaling * _N_A[_qp] *
            (std::exp(_em[_qp]) - _n_gamma) -
        (-1. + 2. * _b) * (-1. + (-1. + _a) * _se_coeff[_qp]) * _secondary_ion / _voltage_scaling *
            _N_A[_qp]) *
       _data.ballast_resist() * (-1. + _r)));

  // std::cout << den0 << ", " << den1 << std::endl;
  std::cout << test0 << ", " << test1 << std::endl;

  return _test[_i][_qp] * _r_units * _eps[_qp] *
         (-2. * (1. + _r) * _u[_qp] - 2. * (1. + _r) * -_V_bat.value(_t, _q_point[_qp]) +
          _data.electrode_area() * _data.coulomb_charge() * _data.ballast_resist() /
              _voltage_scaling * (-1. + _r) *
              ((-1. + (-1. + _a) * _se_coeff[_qp]) * _N_A[_qp] * _ion_drift +
               _N_A[_qp] * (std::exp(_em[_qp]) - _n_gamma) * _v_e_th)) /
         (2. * _data.electrode_area() * _data.coulomb_charge() *
          ((-1. + 2. * _a) * _muem[_qp] / _voltage_scaling * _N_A[_qp] *
               (std::exp(_em[_qp]) - _n_gamma) -
           (-1. + 2. * _b) * (-1. + (-1. + _a) * _se_coeff[_qp]) * _secondary_ion /
               _voltage_scaling * _N_A[_qp]) *
          _data.ballast_resist() * (-1. + _r));
}

