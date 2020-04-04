//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AddDriftDiffusionBase.h"
#include "Parser.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "AddVariableAction.h"
#include "Conversion.h"
#include "DirichletBC.h"
#include "ActionFactory.h"
#include "MooseObjectAction.h"
#include "MooseApp.h"

#include "libmesh/vector_value.h"

#include <sstream>
#include <stdexcept>

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/fe.h"

registerMooseAction("ZapdosApp", AddDriftDiffusionBase, "add_variable");
registerMooseAction("ZapdosApp", AddDriftDiffusionBase, "add_kernel");
registerMooseAction("ZapdosApp", AddDriftDiffusionBase, "add_bc");
registerMooseAction("ZapdosApp", AddDriftDiffusionBase, "add_material");

template <>
InputParameters
validParams<AddDriftDiffusionBase>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

  InputParameters params = validParams<Action>();
  params.addRequiredParam<std::vector<NonlinearVariableName>>(
      "species", "The names of the species that should be added.");
  params.addRequiredParam<std::vector<Real>>(
      "mass", "The mass of each species in kg or amu. (If amu is used, set mass_units = amu).");
  params.addParam<std::string>(
      "mass_units", "kg", "The units of mass to be used, either kg or amu. Default: kg.");
  params.addRequiredParam<std::vector<Real>>("charge", "The charge of each species.");
  params.addRequiredParam<std::vector<VariableName>>(
      "potential", "A dummy vector that holds the potential to couple in for advection");
  params.addRequiredParam<std::string>("potential_units", "The units of potential.");
  params.addParam<Real>("position_units", 1.0, "The units of position.");
  params.addParam<Real>("offset", 31, "The offset for log stabilization.");
  params.addParam<std::vector<Real>>(
      "offset_vals",
      "A vector of offset values. Use if you want different "
      "variables to have different log stabilization terms. Note that this overrides the 'offset' "
      "input parameter if both are set.");
  params.addParam<Real>("reflection",
                        0,
                        "The reflection coefficient. Note that this assumes that there is only one "
                        "reflection coefficient shared between all species. If each species has a "
                        "different coefficient, use reflection_vals input parameter instead.");
  params.addParam<std::vector<Real>>(
      "reflection_vals",
      "A vector of reflection coefficients. Use if you want different variables to have different "
      "reflection coefficients. Note that this overrides the 'reflection' input parameter if both "
      "are set.");
  params.addParam<std::vector<VariableName>>(
      "no_log_stabilization",
      "A list of variables for which the log stabilization term will be excluded.");
  params.addParam<bool>(
      "add_bcs", true, "Whether or not to add BCs. Set to false if you want to add them manually.");
  params.addParam<std::vector<NonlinearVariableName>>(
      "no_bcs",
      "A list of variables for which bc terms will be excluded. Use if you want to manually add "
      "BCs for specific species.");
  params.addParam<std::string>(
      "advection_bc",
      "ADHagelaarIonAdvectionBC",
      "The advection BC to add. Note that this BC will be added to every variable included under "
      "the species parameter. It is recommended to group variables that require the different BCs "
      "in separate blocks.");
  params.addParam<std::string>(
      "diffusion_bc",
      "ADHagelaarIonDiffusionBC",
      "The diffusion BC to add. Note that this BC will be added to every variable included under "
      "the species parameter. It is recommended to group variables that require the different BCs "
      "in separate blocks.");
  params.addParam<bool>("skip_advection", false, "Whether or not to skip advection BCs.");
  params.addParam<bool>("skip_diffusion", false, "Whether or not to skip diffusion BCs.");
  params.addParam<std::vector<BoundaryName>>("boundaries",
                                             "The boundaries that the BCs should act on.");
  params.addParam<std::vector<SubdomainName>>("block",
                                              "The subdomain that this action applies to.");

  return params;
}

AddDriftDiffusionBase::AddDriftDiffusionBase(InputParameters params)
  : Action(params),
    _variables(getParam<std::vector<NonlinearVariableName>>("species")),
    _charge(getParam<std::vector<Real>>("charge")),
    _mass(getParam<std::vector<Real>>("mass")),
    _potential(getParam<std::vector<VariableName>>("potential")),
    _add_bc(getParam<bool>("add_bcs")),
    _block(getParam<std::vector<SubdomainName>>("block")),
    _advection_bc(getParam<std::string>("advection_bc")),
    _diffusion_bc(getParam<std::string>("diffusion_bc")),
    _skip_advection(getParam<bool>("skip_advection")),
    _skip_diffusion(getParam<bool>("skip_diffusion")),
    _boundaries(getParam<std::vector<BoundaryName>>("boundaries"))
{
  if (_variables.size() != _charge.size())
    mooseError("Species list and charge list are not of equal size! Each species needs to have a "
               "charge associated with it.");
  if (_variables.size() != _mass.size())
    mooseError("Species list and mass list are not of equal size! Each species needs to have a "
               "mass associated with it.");

  if (isParamValid("offset_vals"))
  {
    _offset = getParam<std::vector<Real>>("offset_vals");
  }
  else
  {
    _offset.resize(_variables.size(), getParam<Real>("offset"));
  }

  if (isParamValid("reflection_vals"))
  {
    _reflection = getParam<std::vector<Real>>("reflection_vals");
  }
  else
  {
    _reflection.resize(_variables.size(), getParam<Real>("reflection"));
  }

  // Here, we check to see if the variable is included in no_stabilization parameter.
  // If it is, the log stabilization term for that variable is skipped.
  if (isParamValid("no_log_stabilization"))
  {
    //_log_ignore = skip_variable(getParam<std::vector<VariableName>>("no_log_stabilization"));
    auto skip_log = getParam<std::vector<VariableName>>("no_log_stabilization");

    std::vector<VariableName>::iterator iter;
    for (unsigned int i = 0; i < _variables.size(); ++i)
    {
      iter = std::find(skip_log.begin(), skip_log.end(), _variables[i]);
      if (iter != skip_log.end())
        _log_ignore.push_back(true);
      else
        _log_ignore.push_back(false);
    }
  }
  else
  {

    _log_ignore.resize(_variables.size(), false);
  }

  // Do the same for BCs, as long as add_bc is true.
  // (It is set to true by default.
  if (_add_bc)
  {
    if (!isParamValid("boundaries"))
      mooseError(
          "The boundaries that the BCs should be applied to are not included in the DriftDiffusion "
          "action's input parameters! Please add the boundaries that require BCs. Set "
          "add_bcs = false if you want to add BCs manually.");
    if (isParamValid("no_bcs"))
    {
      auto skip_bc = getParam<std::vector<VariableName>>("no_bcs");

      std::vector<VariableName>::iterator iter;
      for (unsigned int i = 0; i < _variables.size(); ++i)
      {
        iter = std::find(skip_bc.begin(), skip_bc.end(), _variables[i]);
        if (iter != skip_bc.end())
          _bc_ignore.push_back(true);
        else
          _bc_ignore.push_back(false);
      }
    }
    else
    {
      _bc_ignore.resize(_variables.size(), false);
    }
  }
}

std::set<SubdomainID>
AddDriftDiffusionBase::getSubdomainIDs()
{
  // Extract and return the block ids supplied in the input
  std::set<SubdomainID> blocks;
  std::vector<SubdomainName> block_param = getParam<std::vector<SubdomainName>>("block");
  for (const auto & subdomain_name : block_param)
  {
    SubdomainID blk_id = _problem->mesh().getSubdomainID(subdomain_name);
    blocks.insert(blk_id);
  }
  return blocks;
}

void
AddDriftDiffusionBase::act()
{
  // MooseSharedPointer<Action> action;
  // MooseSharedPointer<MooseObjectAction> moose_object_action;

  unsigned int number = _variables.size();

  if (_current_task == "add_kernel")
  {
    std::cout << "Adding Kernels..." << std::endl;
    for (unsigned int cur_num = 0; cur_num < number; cur_num++)
    {
      std::string var_name = _variables[cur_num];

      // First we add EField Advection term (only for charged species)
      if (_charge[cur_num] != 0)
      {
        addAdvection("ADEFieldAdvection", _variables[cur_num]);

        addChargeSource("ChargeSourceMoles_KV", _variables[cur_num]);
      }

      // Now we add diffusion term (all species)
      addDiffusion("ADCoeffDiffusion", _variables[cur_num]);

      // Add time derivative for everything
      addTimeDerivative("ADTimeDerivativeLog", _variables[cur_num]);

      // Add log stabilization terms
      if (!_log_ignore[cur_num])
      {
        addLogStab("LogStabilizationMoles", _variables[cur_num], _offset[cur_num]);
      }
    }
  }

  if (_current_task == "add_bc")
  {
    if (_add_bc)
    {
      std::cout << "Adding BCs..." << std::endl;
      for (unsigned int cur_num = 0; cur_num < number; cur_num++)
      {
        for (unsigned int boundary_id = 0; boundary_id < _boundaries.size(); ++boundary_id)
        {
          // First we add advection BCs
          if (!_skip_advection && _charge[cur_num] != 0 && !_bc_ignore[cur_num])
            addBC(_advection_bc,
                  _variables[cur_num],
                  _reflection[cur_num],
                  "advection",
                  _boundaries[boundary_id]);

          // Now we add diffusion BCs
          if (!_skip_diffusion && !_bc_ignore[cur_num])
            addBC(_diffusion_bc,
                  _variables[cur_num],
                  _reflection[cur_num],
                  "diffusion",
                  _boundaries[boundary_id]);
        }
      }
    }
  }

  // mass and charge to be used here
  if (_current_task == "add_material")
    std::cout << "Materials coming soon to an action near you..." << std::endl;
}

/*
std::vector<bool>
AddDriftDiffusionBase::skip_variable(std::vector<VariableName> variables_to_skip)
{
  std::vector<bool> var_ignore;
  std::vector<VariableName>::iterator iter;

  var_ignore.resize(_variables.size());
  for (unsigned int i = 0; i < _variables.size(); ++i)
  {
    iter = std::find(variables_to_skip.begin(), variables_to_skip.end(), _variables[i]);
    if (iter != variables_to_skip.end())
      var_ignore.push_back(true);
    else
      var_ignore.push_back(false);
  }

  return var_ignore;
}
*/

void
AddDriftDiffusionBase::addBC(const std::string & kernel_name,
                             const std::string & var_name,
                             const Real & reflection_value,
                             const std::string & bc_type,
                             const BoundaryName & boundary)
{

  std::string kernel_id = var_name + "_" + bc_type + "_bc";

  auto params = _factory.getValidParams(kernel_name + "<RESIDUAL>");
  params.set<NonlinearVariableName>("variable") = var_name;
  params.set<Real>("position_units") = getParam<Real>("position_units");
  params.set<std::vector<BoundaryName>>("boundary") = {boundary};
  if (params.isParamRequired("r"))
    params.set<Real>("r") = reflection_value;
  if (bc_type == "advection")
    params.set<std::vector<VariableName>>("potential") = _potential;

  _problem->addBoundaryCondition(
      kernel_name + "<RESIDUAL>", kernel_id + boundary + "_residual", params);
  _problem->addBoundaryCondition(
      kernel_name + "<JACOBIAN>", kernel_id + boundary + "_jacobian", params);
  _problem->haveADObjects(true);
}

void
AddDriftDiffusionBase::addLogStab(const std::string & kernel_name,
                                  const std::string & var_name,
                                  const Real & offset_value)
{
  std::string kernel_id = var_name + "_log_stabilization";

  auto params = _factory.getValidParams(kernel_name);
  params.set<NonlinearVariableName>("variable") = var_name;
  params.set<Real>("offset") = offset_value;
  if (isParamValid("block"))
    params.set<std::vector<SubdomainName>>("block") = _block;

  // No AD version of charged source
  _problem->addKernel(kernel_name, kernel_id, params);
}

void
AddDriftDiffusionBase::addChargeSource(const std::string & kernel_name,
                                       const std::string & var_name)
{
  std::string kernel_id = var_name + "_charge_source";

  auto params = _factory.getValidParams(kernel_name);
  params.set<NonlinearVariableName>("variable") = _potential[0];
  params.set<std::vector<VariableName>>("charged") = {var_name};
  if (isParamValid("block"))
    params.set<std::vector<SubdomainName>>("block") = _block;
  params.set<std::string>("potential_units") = getParam<std::string>("potential_units");

  // No AD version of charged source
  _problem->addKernel(kernel_name, kernel_id, params);
}

void
AddDriftDiffusionBase::addAdvection(const std::string & kernel_name, const std::string & var_name)
{
  std::string kernel_id = var_name + "_efield_advection";

  auto params = _factory.getValidParams(kernel_name + "<RESIDUAL>");
  params.set<NonlinearVariableName>("variable") = var_name;
  params.set<std::vector<VariableName>>("potential") = _potential;
  if (isParamValid("block"))
    params.set<std::vector<SubdomainName>>("block") = _block;
  params.set<Real>("position_units") = getParam<Real>("position_units");

  _problem->addKernel(kernel_name + "<RESIDUAL>", kernel_id + "_residual", params);
  _problem->addKernel(kernel_name + "<JACOBIAN>", kernel_id + "_jacobian", params);
  _problem->haveADObjects(true);
}

void
AddDriftDiffusionBase::addDiffusion(const std::string & kernel_name, const std::string & var_name)
{
  std::string kernel_id = var_name + "_diffusion";

  auto params = _factory.getValidParams(kernel_name + "<RESIDUAL>");
  params.set<NonlinearVariableName>("variable") = var_name;
  if (isParamValid("block"))
    params.set<std::vector<SubdomainName>>("block") = _block;
  params.set<Real>("position_units") = getParam<Real>("position_units");

  _problem->addKernel(kernel_name + "<RESIDUAL>", kernel_id + "_residual", params);
  _problem->addKernel(kernel_name + "<JACOBIAN>", kernel_id + "_jacobian", params);
  _problem->haveADObjects(true);
}

void
AddDriftDiffusionBase::addTimeDerivative(const std::string & kernel_name,
                                         const std::string & var_name)
{
  std::string kernel_id = var_name + "_time";

  auto params = _factory.getValidParams(kernel_name + "<RESIDUAL>");
  params.set<NonlinearVariableName>("variable") = var_name;
  if (isParamValid("block"))
    params.set<std::vector<SubdomainName>>("block") = _block;

  _problem->addKernel(kernel_name + "<RESIDUAL>", kernel_id + "_residual", params);
  _problem->addKernel(kernel_name + "<JACOBIAN>", kernel_id + "_jacobian", params);
  _problem->haveADObjects(true);
}
