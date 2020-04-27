#include "BoundaryFlux.h"

registerMooseObject("ZapdosApp", BoundaryFlux);

template <>
InputParameters
validParams<BoundaryFlux>()
{
  InputParameters params = validParams<ODEKernel>();
  params.addRequiredParam<PostprocessorName>(
      "current",
      "The postprocessor response for calculating the current passing through the needle surface.");
  return params;
}

BoundaryFlux::BoundaryFlux(const InputParameters & parameters)
  : ODEKernel(parameters), _current(getPostprocessorValue("current"))
{
}

Real
BoundaryFlux::computeQpResidual()
{
  return -_current*0.001;
}

Real
BoundaryFlux::computeQpJacobian()
{
  // return -_stoichiometric_coeff * _rate_coefficient[_i];
  return 0.0;
}
