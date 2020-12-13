//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FunctionFlux.h"
#include "Function.h"

registerADMooseObject("ZapdosApp", FunctionFlux);

defineADLegacyParams(FunctionFlux);

InputParameters
FunctionFlux::validParams()
{
  InputParameters params = ADIntegratedBC::validParams();
  params.addRequiredParam<FunctionName>("function", "The value of the input flux.");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  params.addClassDescription("Implements a constant flux boundary condition.");
  return params;
}

FunctionFlux::FunctionFlux(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
    _r_units(1. / getParam<Real>("position_units")),
    _func(getFunction("function"))
{
}

ADReal
FunctionFlux::computeQpResidual()
{
  return _test[_i][_qp] * _r_units * _func.value(_t, _q_point[_qp]);
}
