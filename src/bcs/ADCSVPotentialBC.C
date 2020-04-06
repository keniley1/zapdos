//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADCSVPotentialBC.h"
#include "MooseUtils.h"

registerADMooseObject("ZapdosApp", ADCSVPotentialBC);

defineADLegacyParams(ADCSVPotentialBC);

template <ComputeStage compute_stage>
InputParameters
ADCSVPotentialBC<compute_stage>::validParams()
{
  InputParameters params = ADDirichletBCBase<compute_stage>::validParams();
  params.addClassDescription("Imposes the essential boundary condition $u=g$, where $g$ "
                             "is calculated by a function.");
  /*
  params.addParam<Real>(
      "recurrence_time",
      0,
      "The recurrence time of the voltage data. Example: If recurrence_time = 1e-3 and timesteps "
      "are in units of seconds, the CSV data will be repeated every 1e-3 seconds.");
      */
  // params.addParam<FunctionName>("function", 0, "The function describing the Dirichlet
  // condition");
  //
  params.addRequiredParam<FileName>("file_name", "The file containing voltage versus time data.");
  return params;
}

template <ComputeStage compute_stage>
ADCSVPotentialBC<compute_stage>::ADCSVPotentialBC(const InputParameters & parameters)
  : ADDirichletBCBase<compute_stage>(parameters)
{
  std::vector<Real> x_val;
  std::vector<Real> y_val;

  std::string file_name = getParam<FileName>("file_name");
  MooseUtils::checkFileReadable(file_name);
  const char * charPath = file_name.c_str();
  std::ifstream myfile(charPath);
  Real value;

  if (myfile.is_open())
  {
    while (myfile >> value)
    {
      x_val.push_back(value);
      myfile >> value;
      y_val.push_back(value);
    }
    myfile.close();
  }

  else
    mooseError("Unable to open file");

  //_voltage.setData(x_val, y_val);
  _voltage = libmesh_make_unique<LinearInterpolation>(x_val, y_val);
}

template <ComputeStage compute_stage>
ADReal
ADCSVPotentialBC<compute_stage>::computeQpValue()
{
  return _voltage->sample(_t);
}
