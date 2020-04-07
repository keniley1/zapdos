//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADCoeffDiffusion.h"

registerADMooseObject("ZapdosApp", ADCoeffDiffusion);

defineADLegacyParams(ADCoeffDiffusion);

template <ComputeStage compute_stage>
InputParameters
ADCoeffDiffusion<compute_stage>::validParams()
{
  InputParameters params = ADKernelGrad<compute_stage>::validParams();
  params.addRequiredParam<Real>("position_units", "Units of position.");
  return params;
}

template <ComputeStage compute_stage>
ADCoeffDiffusion<compute_stage>::ADCoeffDiffusion(const InputParameters & parameters)
  : ADKernelGrad<compute_stage>(parameters),
    _r_units(1. / getParam<Real>("position_units")),
    _diffusivity(getADMaterialProperty<Real>("diff" + _var.name()))
{
}

// ADRealVectorValue
// ADCoeffDiffusion<compute_stage>::precomputeQpResidual()
template <ComputeStage compute_stage>
ADRealVectorValue
ADCoeffDiffusion<compute_stage>::precomputeQpResidual()
{
  return _diffusivity[_qp] * std::exp(_u[_qp]) * _grad_u[_qp] * _r_units * _r_units;
}
