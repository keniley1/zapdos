/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "InterpolateAlongAxis.h"
#include "MooseMesh.h"

registerMooseObject("ZapdosApp", InterpolateAlongAxis);

template <>
InputParameters
validParams<InterpolateAlongAxis>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addParam<unsigned>("axis", 0, "The axis to interpolate along.");
  params.addParam<FileName>("property_file", "", "The file containing interpolation table.");
  params.addParam<std::string>(
      "file_location", "", "The name of the file that stores the reaction rate tables.");
  params.addParam<std::string>("sampling_format",
                               "reduced_field",
                               "The format that the rate constant files are in. Options: "
                               "reduced_field and electron_energy.");
  return params;
}

InterpolateAlongAxis::InterpolateAlongAxis(const InputParameters & parameters)
  : AuxKernel(parameters), _component(getParam<unsigned>("axis"))
{
  std::vector<Real> x_val;
  std::vector<Real> y_val;
  std::string file_name =
      getParam<std::string>("file_location") + "/" + getParam<FileName>("property_file");
  MooseUtils::checkFileReadable(file_name);
  const char * charPath = file_name.c_str();
  std::ifstream myfile(charPath);
  Real value;

  int i = 0;
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

  _coefficient_interpolation.setData(x_val, y_val);
}

Real
InterpolateAlongAxis::computeValue()
{
  Real val = 0;

  // std::cout << (*_current_node)(_component) << std::endl;
  // return _coefficient_interpolation.sample(_q_point[_qp](_component));
  return _coefficient_interpolation.sample((*_current_node)(_component));
}
