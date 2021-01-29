//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VaporPressureBC.h"

#include "libmesh/node.h"
#include "Function.h"
#include <iostream>

registerADMooseObject("ZapdosApp", VaporPressureBC);
InputParameters
VaporPressureBC::validParams(){
  InputParameters params = ADNodalBC::validParams();
  params.addRequiredCoupledVar("gas_temp", "The temperature of the gas.");
  //^assign Tg to this variable above
  return params;
}
VaporPressureBC::VaporPressureBC(const InputParameters & parameters)
  : ADNodalBC(parameters),
  _gas_temp(adCoupledValue("gas_temp"))
{
}

ADReal
VaporPressureBC::computeQpResidual(){
  //we want to return the log-mole density of the water
  //using gasTemp, we can estimate vapor pressure, then use ideal gas law
  std::cout << "gas temp: " << MetaPhysicL::raw_value(_gas_temp[_qp]) << " Kelvin" << std::endl;
  _vapor_pressure = 0.61094 * std::exp((17.625 * (_gas_temp[_qp] - 273.15)) / ((_gas_temp[_qp] - 273.15) + 243.04));
  std::cout << "vapor pressure: " <<MetaPhysicL::raw_value( _vapor_pressure) << "kPa" << std::endl;
  //the ARM equation uses Celsius and KiloPascals
  _atoms_H2O = (_vapor_pressure * 1000) / (1.38 * std::pow(10, -23) * _gas_temp[_qp]);
  std::cout << "atoms of H2O: " << MetaPhysicL::raw_value(_atoms_H2O) << std::endl;
  _units_H2O = std::log(_atoms_H2O / (6.022 * std::pow(10, 23)));
  std::cout << "log mole of H2O: " << MetaPhysicL::raw_value(_units_H2O) << std::endl;


  return _u -_units_H2O;
}
