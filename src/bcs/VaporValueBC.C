//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VaporValueBC.h"

registerMooseObject("ZapdosApp", VaporValueBC);

InputParameters
VaporValueBC::validParams()
{
  InputParameters params = ADDirichletBCBase::validParams();
  params.addRequiredCoupledVar("gas_temperature",
                               "The gas temperature in the system (in K). This value "
                               "is expected to be a total gas temperature, "
                               "not the temperature of an individual species.");
  // params.declareControllable("value");
  params.addClassDescription(
      "Imposes a Dirichlet BC based on the vapor pressure of water, "
      "described by August-Roche-Magnus equation. (Temperature is in KELVIN, not CELSIUS.)");
  return params;
}

VaporValueBC::VaporValueBC(const InputParameters & parameters)
  : ADDirichletBCBase(parameters), _gas_temp(adCoupledValue("gas_temperature"))
{
}

ADReal
VaporValueBC::computeQpValue()
{
  _vapor_pressure = 0.61094 * std::exp((17.625 * (_gas_temp[_qp] - 273.15)) /
                                       ((_gas_temp[_qp] - 273.15) + 243.04));

  return std::log((_vapor_pressure * 1000.) / (1.38e-23 * _gas_temp[_qp]) / 6.022e23);
}
