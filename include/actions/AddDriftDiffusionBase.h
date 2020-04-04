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

#include "AddVariableAction.h"
#include "Action.h"

class AddDriftDiffusionBase;

template <>
InputParameters validParams<AddDriftDiffusionBase>();

class AddDriftDiffusionBase : public Action
{
protected:
  // virtual std::vector<bool> skip_variable(std::vector<VariableName> variables_to_skip);
  std::set<SubdomainID> getSubdomainIDs();
  virtual void addBC(const std::string & kernel_name,
                     const std::string & var_name,
                     const Real & reflection_value,
                     const std::string & bc_type,
                     const BoundaryName & boundary);
  virtual void addLogStab(const std::string & kernel_name,
                          const std::string & var_name,
                          const Real & offset_value);
  virtual void addChargeSource(const std::string & kernel_name, const std::string & var_name);
  virtual void addAdvection(const std::string & kernel_name, const std::string & var_name);
  virtual void addDiffusion(const std::string & kernel_name, const std::string & var_name);
  virtual void addTimeDerivative(const std::string & kernel_name, const std::string & var_name);

public:
  AddDriftDiffusionBase(InputParameters params);

  const std::vector<NonlinearVariableName> _variables;
  const std::vector<Real> _charge;
  const std::vector<Real> _mass;
  const std::vector<VariableName> _potential;
  const bool _add_bc;
  const std::vector<SubdomainName> _block;
  const std::string _advection_bc;
  const std::string _diffusion_bc;
  const bool _skip_advection;
  const bool _skip_diffusion;
  const std::vector<BoundaryName> _boundaries;

  std::vector<Real> _offset;
  std::vector<Real> _reflection;
  std::vector<bool> _log_ignore;
  std::vector<bool> _bc_ignore;

  virtual void act();
};

