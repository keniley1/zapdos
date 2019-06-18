//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "HeavySpeciesMod.h"
#include "MooseUtils.h"

// MOOSE includes
#include "MooseVariable.h"

registerMooseObject("ZapdosApp", HeavySpeciesMod);

template <>
InputParameters
validParams<HeavySpeciesMod>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredParam<std::string>("heavy_species_name", "The name of the heavy species");
  params.addRequiredParam<Real>("heavy_species_mass", "Mass of the heavy species");
  params.addRequiredParam<std::string>("potential_units", "The potential units.");
  params.addRequiredParam<Real>("heavy_species_charge", "Charge of heavy species.");
  params.addParam<Real>("time_units", 1, "Units of time");
  params.addParam<Real>("mobility", "The species mobility (if applicable).");
  params.addParam<Real>("diffusivity", "The species diffusivity (if applicable).");
  // params.addRequiredParam<FileName>(
      // "reactions_file", "The file containing interpolation tables for material properties.");
  // params.addParam<Real>("user_T_gas", 300, "The gas temperature in Kelvin.");
  // params.addParam<Real>("user_p_gas", 1.01e5, "The gas pressure in Pascals.");
  return params;
}

HeavySpeciesMod::HeavySpeciesMod(const InputParameters & parameters)
  : Material(parameters),
    _user_massHeavy(getParam<Real>("heavy_species_mass")),
    _user_sgnHeavy(getParam<Real>("heavy_species_charge")),
    _potential_units(getParam<std::string>("potential_units")),
    _massHeavy(declareProperty<Real>("mass" + getParam<std::string>("heavy_species_name"))),
    _temperatureHeavy(declareProperty<Real>("T" + getParam<std::string>("heavy_species_name"))),
    _sgnHeavy(declareProperty<Real>("sgn" + getParam<std::string>("heavy_species_name"))),
    _muHeavy(declareProperty<Real>("mu" + getParam<std::string>("heavy_species_name"))),
    _diffHeavy(declareProperty<Real>("diff" + getParam<std::string>("heavy_species_name"))),
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
    mooseError("User must input AT LEAST ONE of the parameters 'mobility' and 'diffusivity' to HeavySpeciesMod!"); 
  }

}

void
HeavySpeciesMod::computeQpProperties()
{
  if (_potential_units.compare("V") == 0)
    _voltage_scaling = 1.;
  else if (_potential_units.compare("kV") == 0)
    _voltage_scaling = 1000;

  _massHeavy[_qp] = _user_massHeavy;
  _sgnHeavy[_qp] = _user_sgnHeavy;

  _temperatureHeavy[_qp] = 300;  // Kelvin - this is temporary

  if (!_calc_mobility && !_calc_diffusivity)
  {
    _diffHeavy[_qp] = getParam<Real>("diffusivity") * _time_units;
    _muHeavy[_qp] = getParam<Real>("mobility") * _voltage_scaling * _time_units;
  }
  else if (_calc_mobility && !_calc_diffusivity)
  {
    _diffHeavy[_qp] = getParam<Real>("diffusivity") * _time_units;
    _muHeavy[_qp] = _diffHeavy[_qp] * 1.602e-19 / (_temperatureHeavy[_qp] * 1.3807e-23) * _voltage_scaling; 
  } 
  else if (!_calc_mobility && _calc_diffusivity)
  {
    _muHeavy[_qp] = getParam<Real>("mobility") * _voltage_scaling * _time_units;
    _diffHeavy[_qp] = _muHeavy[_qp] * _temperatureHeavy[_qp] * 1.3807e-23 / (_user_sgnHeavy * 1.602e-19);
  }
  else
  {
    mooseError("This shouldn't happen!");
  }



  /*
  if (isParamValid("mobility"))
  {
    _muHeavy[_qp] = getParam<Real>("mobility") * _voltage_scaling * _time_units;
  }
  else
  {
    _muHeavy[_qp] =
        1444. * _voltage_scaling * _time_units /
        (10000. * 760. * _p_gas[_qp] / 1.01E5); // units of m^2/(kV*s) if _voltage_scaling = 1000
  }

  if (isParamValid("diffusivity"))
  {
    _diffHeavy[_qp] = getParam<Real>("diffusivity") * _time_units;
  }
  else
  {
    _diffHeavy[_qp] = 0.004 * _time_units / (760. * _p_gas[_qp] / 1.01E5); // covert to m^2 and include press
  }
  */
}
