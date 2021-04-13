//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeavySpeciesMaterial.h"
#include "MooseUtils.h"

// MOOSE includes
#include "MooseVariable.h"

registerMooseObject("ZapdosApp", HeavySpeciesMaterial);
registerMooseObject("ZapdosApp", ADHeavySpeciesMaterial);


/*****
 */
template <bool is_ad>
InputParameters
HeavySpeciesMaterialTempl<is_ad>::validParams()
{
  InputParameters params = Material::validParams();

  params.addRequiredParam<std::string>("heavy_species_name", "The name of the heavy species");
  params.addRequiredParam<Real>("heavy_species_mass", "Mass of the heavy species");
  params.addRequiredParam<std::string>("potential_units", "The potential units.");
  params.addRequiredParam<Real>("heavy_species_charge", "Charge of heavy species.");
  params.addParam<Real>("time_units", 1, "Units of time");
  params.addParam<Real>("mobility", "The species mobility (if applicable).");
  params.addParam<Real>("diffusivity", "The species diffusivity (if applicable).");
  params.addParam<bool>(
      "species_temperature",
      false,
      "Whether this individual species temperature is being tracked as a nonlinear variable. "
      "Default is false. (Note that this feature is not yet supported, and an error will return if "
      "true.");
  params.addCoupledVar(
      "gas_temperature",
      "The gas temperature of the system, in Kelvin. Note that this is assumed to be the total gas "
      "temperature, not the temperature of an individual species.");
  params.addParam<Real>(
      "temperature_scale",
      1,
      "Whether or not to scale diffusivity/mobility based on temperature. This "
      "will take the user-provided diffusivity and/or mobility and scale it by "
      "this value^(-1/3). Transport coefficients typically have a cube root "
      "dependence on temperature (see: Chapman-Enskog theory), but they are frequently reported in "
      "literature as constant values at 300 K. If the user "
      "already took this into account, this option may be set to 1.");
  return params;
}

template <bool is_ad>
HeavySpeciesMaterialTempl<is_ad>::HeavySpeciesMaterialTempl(const InputParameters & parameters)
  : Material(parameters),
    _user_massHeavy(getParam<Real>("heavy_species_mass")),
    _user_sgnHeavy(getParam<Real>("heavy_species_charge")),
    _potential_units(getParam<std::string>("potential_units")),
    _massHeavy(declareProperty<Real>("mass" + getParam<std::string>("heavy_species_name"))),
    _temperatureHeavy(
        declareGenericProperty<Real, is_ad>("T" + getParam<std::string>("heavy_species_name"))),
    _sgnHeavy(declareProperty<Real>("sgn" + getParam<std::string>("heavy_species_name"))),
    _muHeavy(
        declareGenericProperty<Real, is_ad>("mu" + getParam<std::string>("heavy_species_name"))),
    _diffHeavy(
        declareGenericProperty<Real, is_ad>("diff" + getParam<std::string>("heavy_species_name"))),
    _T_gas(getMaterialProperty<Real>("T_gas")),
    _p_gas(getMaterialProperty<Real>("p_gas")),
    _gas_temp(isCoupled("gas_temperature") ? coupledValue("gas_temperature") : _zero),
    //_temp_scale(getParam<Real>("temperature_scale")),
    _time_units(getParam<Real>("time_units"))

{
  if (isParamValid("mobility") && isParamValid("diffusivity"))
  {
    _calc_mobility = false;
    _calc_diffusivity = false;
  }
  else if (isParamValid("mobility") && !isParamValid("diffusivity"))
  {
    _calc_mobility = false;
    _calc_diffusivity = true;
  }
  else if (!isParamValid("mobility") && isParamValid("diffusivity"))
  {
    _calc_mobility = true;
    _calc_diffusivity = false;
  }
  else
  {
    _calc_mobility = true;
    _calc_diffusivity = true;
  }

  // Several options are available for temperature dependence.
  // First, each species may have its own temperature and associated energy equation. (note that
  // this is not yet implemented as of March 2021.)
  // Second, gas temperature is coupled. (These two options are mutually exclusive.)
  //
  // In either case, temperature scaling is provided as an option.
  // temperature_scale will scale values for diffusivity and/or mobility by the indicated
  // value to the -1/3 power (transport coefficients depend on
  // the cube root of the temperature based on Chapman-Enskog theory).
  //
  // This option is included because transport coefficients are frequently provided as constant
  // values at a specific temperature. This option allows the user to scale these values to the
  // temperature used in this simulation. 
  //
  // Note: temperature_scale defaults to 1, e.g. no scaling.
  _temp_scale = std::pow(getParam<Real>("temperature_scale"), -3. / 2.);
  /*
  if (getParam<bool>("species_temperature"))
    mooseError(
        "The ability to track an individual species temperature is not yet implemented. Please set "
        "species_temperature to false for all HeavySpeciesMaterials added to this input file.");
  else if (isCoupled("gas_temperature"))
  {
  }
  else
  {
  }
  */
}

template <bool is_ad>
void
HeavySpeciesMaterialTempl<is_ad>::computeQpProperties()
{
  if (_potential_units.compare("V") == 0)
    _voltage_scaling = 1.;
  else if (_potential_units.compare("kV") == 0)
    _voltage_scaling = 1000;

  _massHeavy[_qp] = _user_massHeavy;
  _sgnHeavy[_qp] = _user_sgnHeavy;

  // _T_gas[_qp] = _user_T_gas;
  // _p_gas[_qp] = _user_p_gas;

  //_temperatureHeavy[_qp] = _T_gas[_qp]; // Needs to be changed.
  if (isCoupled("gas_temperature"))
    _temperatureHeavy[_qp] = _gas_temp[_qp];
  else
    _temperatureHeavy[_qp] = _T_gas[_qp]; // Needs to be changed.

  // _n_gas[_qp] = _p_gas[_qp] / (8.3145 * _T_gas[_qp]);

  if (!_calc_mobility && !_calc_diffusivity)
  {
    _diffHeavy[_qp] = getParam<Real>("diffusivity") * _time_units;
    _muHeavy[_qp] = getParam<Real>("mobility") * _voltage_scaling * _time_units;
  }
  else if (_calc_mobility && !_calc_diffusivity)
  {
    //_diffHeavy[_qp] = getParam<Real>("diffusivity") * _time_units / _temp_scale * std::pow(_temperatureHeavy[_qp], 3./2.);
    _diffHeavy[_qp] = getParam<Real>("diffusivity") * _time_units;
    _muHeavy[_qp] =
        _diffHeavy[_qp] * 1.602e-19 / (_temperatureHeavy[_qp] * 1.3807e-23) * _voltage_scaling;
  }
  else if (!_calc_mobility && _calc_diffusivity)
  {
    _muHeavy[_qp] = getParam<Real>("mobility") * _voltage_scaling * _time_units;
    _diffHeavy[_qp] = _muHeavy[_qp] * _temperatureHeavy[_qp] * 1.3807e-23 /
                      (_user_sgnHeavy * 1.602e-19) / _voltage_scaling;
  }
  else
  {
    // If no mobility or diffusivity values are given, values are computed for Argon based on
    // Richards and Sawin paper
    /*
    _muHeavy[_qp] =
        1444. * _voltage_scaling * _time_units /
        (10000. * 760. * _p_gas[_qp] / 1.01E5); // units of m^2/(kV*s) if _voltage_scaling = 1000
        */

    _diffHeavy[_qp] =
        0.004 * _time_units / (760. * _p_gas[_qp] / 1.01E5); // covert to m^2 and include press

    _muHeavy[_qp] = _diffHeavy[_qp] * 1.602e-19 / (_temperatureHeavy[_qp] * 1.3807e-23) * _voltage_scaling;
    /*
    _diffHeavy[_qp] =
        0.004 * _time_units / (760. * _p_gas[_qp] / 1.01E5) * _temp_scale * std::pow(_temperatureHeavy[_qp], 3./2.); // covert to m^2 and include press

    _muHeavy[_qp] = _diffHeavy[_qp] * 1.602e-19 / (_temperatureHeavy[_qp] * 1.3807e-23) * _voltage_scaling;
    */
  }
}

template class HeavySpeciesMaterialTempl<false>;
template class HeavySpeciesMaterialTempl<true>;
