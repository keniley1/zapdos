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

#include "ADKernelGrad.h"

template <ComputeStage>
class ADCoeffDiffusion;

declareADValidParams(ADCoeffDiffusion);

template <ComputeStage compute_stage>
class ADCoeffDiffusion : public ADKernelGrad<compute_stage>
{
public:
  ADCoeffDiffusion(const InputParameters & parameters);
  static InputParameters validParams();

protected:
  virtual ADRealVectorValue precomputeQpResidual() override;

  usingKernelGradMembers;

private:
  const Real _r_units;

  const ADMaterialProperty(Real) & _diffusivity;
};
/*
template <ComputeStage compute_stage>

class ADCoeffDiffusion : public ADKernelGrad<compute_stage>
{
public:
  static InputParameters validParams();

  ADCoeffDiffusion(const InputParameters & parameters);

protected:
  virtual ADRealVectorValue precomputeQpResidual() override;
  // virtual ADReal computeQpResidual();

  usingKernelGradMembers;
  // using ADKernelGrad<compute_stage>::getPostprocessorValue;
  //usingKernelMembers;

private:
  /// Position units
  const Real _r_units;

  /// The diffusion coefficient (either constant or mixture-averaged)
  const ADMaterialProperty(Real) & _diffusivity;
};
*/
