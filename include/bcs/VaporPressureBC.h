//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADNodalBC.h"
#include "MooseVariable.h"


class VaporPressureBC : public ADNodalBC
{
public:
  static InputParameters validParams();

  VaporPressureBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
  //ADReal computeQpValue();
  const ADVariableValue & _gas_temp;
  ADReal _vapor_pressure;
  ADReal _atoms_H2O;
  ADReal _units_H2O;



};
