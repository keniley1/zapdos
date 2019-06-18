#pragma once

#include "InterfaceKernel.h"

class InterfaceFluxConservation;

template <>
InputParameters validParams<InterfaceFluxConservation>();

class InterfaceFluxConservation : public InterfaceKernel
{
public:
  InterfaceFluxConservation(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual(Moose::DGResidualType type);
  virtual Real computeQpJacobian(Moose::DGJacobianType type);

  Real _r_units;
  Real _r_neighbor_units;

  const MaterialProperty<Real> & _diffusivity_main;
  const MaterialProperty<Real> & _diffusivity_neighbor;
  // const MaterialProperties<Real> & sigma;
};

