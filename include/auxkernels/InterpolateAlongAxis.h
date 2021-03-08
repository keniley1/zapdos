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

#pragma once

#include "AuxKernel.h"
#include "LinearInterpolation.h"

class InterpolateAlongAxis;

template <>
InputParameters validParams<InterpolateAlongAxis>();

class InterpolateAlongAxis : public AuxKernel
{
public:
  InterpolateAlongAxis(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  LinearInterpolation _coefficient_interpolation;

  unsigned _component;
};
