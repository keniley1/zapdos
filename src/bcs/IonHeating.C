//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "IonHeating.h"

// MOOSE includes
#include "MooseVariable.h"
#include "Function.h"

registerADMooseObject("ZapdosApp", IonHeating);

defineADLegacyParams(IonHeating);

InputParameters
IonHeating::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addRequiredCoupledVar("ions", "The ion density variabels.");
  //params.addRequiredCoupledVar("em", "The log of the electron density.");
  //params.addRequiredCoupledVar(
  //    "mean_en", "The log of the product of the mean energy and the electron density.");
  params.addRequiredCoupledVar("potential", "The electric potential variable.");
  params.addRequiredParam<std::string>("potential_units", "The potential units.");
  params.addRequiredParam<Real>("r",
                                "The reflection coefficient applied to both electrons and ions");
  params.addRequiredParam<Real>("position_units", "Units of position.");

  return params;
}

IonHeating::IonHeating(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
    _r_units(1. / getParam<Real>("position_units")),
    //_mean_en(adCoupledValue("mean_en")),
    //_em(adCoupledValue("em")),
    //_se_coeff(getMaterialProperty<Real>("se_coeff")),
    //_eps(getMaterialProperty<Real>("eps")),
    //_N_A(getMaterialProperty<Real>("N_A")),
    //_muem(getADMaterialProperty<Real>("muem")),
    //_e(getMaterialProperty<Real>("e")),
    //_massem(getMaterialProperty<Real>("massem")),
    _grad_potential(adCoupledGradient("potential")),
    _kb(getMaterialProperty<Real>("k_boltz")),

    _work_function(getMaterialProperty<Real>("work_function")),
    _potential_units(getParam<std::string>("potential_units")),
    _r(getParam<Real>("r"))
{
  //_ion_flux.zero();
  //_n_gamma = 0.0;
  //_actual_mean_en = 0.0;
  //_v_e_th = 0.0;
  _v_i_th = 0.0;
  //_a = 0.0;
  _b = 0.0;
  //_numerator = 0.0;
  //_denominator = 0.0;
  if (_potential_units.compare("V") == 0)
    _voltage_scaling = 1.;
  else if (_potential_units.compare("kV") == 0)
    _voltage_scaling = 1000;

  // First we need to initialize all of the ion densities and material properties.
  // Find the number of ions coupled into this BC:
  _num_ions = coupledComponents("ions");

  // Resize the vectors to store _num_ions values:
  _ion.resize(_num_ions);
  _grad_ion.resize(_num_ions);
  _T_heavy.resize(_num_ions);
  _muion.resize(_num_ions);
  //_Dip.resize(_num_ions);
  _sgnion.resize(_num_ions);
  _mass.resize(_num_ions);

  // Retrieve the values for each ion and store in the relevant vectors.
  // Note that these need to be dereferenced to get the values inside the
  // main body of the code.
  // e.g. instead of "_ion[_qp]" it would be "(*_ion[i])[_qp]", where "i"
  // refers to a single ion species.
  for (unsigned int i = 0; i < _ion.size(); ++i)
  {
    _ion[i] = &adCoupledValue("ions", i);
    _grad_ion[i] = &adCoupledGradient("ions", i);
    _T_heavy[i] = &getADMaterialProperty<Real>("T" + (*getVar("ions", i)).name());
    _muion[i] = &getADMaterialProperty<Real>("mu" + (*getVar("ions", i)).name());
    //_Dip[i] = &getADMaterialProperty<Real>("diff" + (*getVar("ions", i)).name());
    _sgnion[i] = &getMaterialProperty<Real>("sgn" + (*getVar("ions", i)).name());
    _mass[i] = &getMaterialProperty<Real>("mass" + (*getVar("ions", i)).name());
  }
}

ADReal
IonHeating::computeQpResidual()
{
  /*
  if (_normals[_qp] * -1.0 * -_grad_potential[_qp] > 0.0)
  {
    _a = 1.0;
  }
  else
  {
    _a = 0.0;
  }
  */

  //_ion_flux.zero();
  _ion_flux = 0;
  //_n_gamma = 0;
  //_secondary_ion = 0;
  //_ion_drift = 0;
  for (unsigned int i = 0; i < _num_ions; ++i)
  {
    /*
    _ion_flux +=
        (*_sgnion[i])[_qp] * (*_muion[i])[_qp] * -_grad_potential[_qp] * _r_units * std::exp((*_ion[i])[_qp]) +
        0.5 * _vi_thermal;
  }
  */
    if (_normals[_qp] * (*_sgnion[i])[_qp] * -_grad_potential[_qp] > 0.0)
      _b = 1.0;
    else
      _b = 0.0;
    _ion_flux += std::exp((*_ion[i])[_qp]) *
                 (0.5 * std::sqrt(8 * _kb[_qp] * (*_T_heavy[i])[_qp] / (M_PI * (*_mass[i])[_qp])) +
                  (2 * _b - 1) * (*_sgnion[i])[_qp] * (*_muion[i])[_qp] * -_grad_potential[_qp] *
                      _r_units * _normals[_qp]);
  }
  _ion_flux *= (1.0 - _r) / (1.0 + _r) * 1.602e-19 * 6.022e23;
  //_n_gamma = (1. - _a) * _se_coeff[_qp] * _ion_flux * _normals[_qp] /
  //           (_muem[_qp] * -_grad_u[_qp] * _r_units * _normals[_qp]);

  //_v_e_th = std::sqrt(8 * _data.coulomb_charge() * 2.0 / 3 * std::exp(_mean_en[_qp] - _em[_qp]) /
  //                    (M_PI * _massem[_qp]));

  //_electron_flux = - 2./(1. - _r) * _se_coeff * _ion_flux;

  // 15.76 is argon ionziation energy. Needs to be replaced with material property or just real input?
  // Latter would be more optimal, but it would require more input for the user. Actions could take care of it though.
  // 2 is the ion acceleration through the sheath. Just an estimate. Actually a nonlocal phenomenon that would be very hard to capture here.
  return _test[_i][_qp] * _r_units * _ion_flux * (15.76 * 2. - _work_function[_qp]);

  /*
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
          */
}
