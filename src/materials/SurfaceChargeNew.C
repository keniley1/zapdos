//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SurfaceChargeNew.h"
#include "MooseUtils.h"

registerADMooseObject("ZapdosApp", SurfaceChargeNew);

InputParameters
SurfaceChargeNew::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addParam<std::string>("bc_type",
                               "Hagelaar",
                               "The name of the family of BCs being used in this model. "
                               "Options: Hagelaar (DEFAULT), Sakiyama, Lymberopoulos.");
  params.addRequiredCoupledVar("potential", "The potential that drives the advective flux.");
  params.addParam<Real>("r_ion", 0.0, "The ion reflection coefficient on this boundary.");
  params.addParam<Real>("r_electron", 0.0, "The electron reflection coefficient on this boundary.");
  params.addRequiredParam<Real>("position_units", "Units of position.");
  params.addRequiredCoupledVar("species",
                               "All of the charged species that interact with this boundary.");
  params.addCoupledVar("em", "The electron density.");
  params.addCoupledVar("mean_en", "The electron density.");
  params.addParam<bool>("include_electrons", true, "Whether or not electrons are included.");
  params.addParam<bool>(
      "include_secondary_electrons", true, "Whether or not secondary electrons are included.");
  params.addParam<Real>("ks", "The recombination coefficient (for Lymberopoulos-type BC)");
  params.addParam<bool>(
      "secondary_electrons",
      true,
      "Whether or not to include secondary electron emission in the surface charge calculation. "
      "Note that this should be consistent with the selected boundary conditions; if a secondary "
      "electron BC is used on this boundary, this should be true. DEFAULT: true.");
  params.addRequiredParam<std::string>("potential_units", "The potential units.");
  params.addClassDescription(
      "Adds a surface charge material property based on the rate of change of the total charged "
      "flux to a boundary. (NOTE: this material is meant to be boundary-restricted.)");
  return params;
}

SurfaceChargeNew::SurfaceChargeNew(const InputParameters & parameters)
  : ADMaterial(parameters),
    _sigma(declareADProperty<Real>("surface_charge")),
    _sigma_old(getMaterialPropertyOld<Real>("surface_charge")),
    _plasma_current(declareProperty<Real>("plasma_current")),
    _include_electrons(getParam<bool>("include_electrons")),
    _include_secondary_electrons(getParam<bool>("include_secondary_electrons")),
    _electron(isParamValid("em") ? &adCoupledValue("em") : nullptr),
    _mean_energy(isParamValid("mean_en") ? &adCoupledValue("mean_en") : nullptr),
    _mass_em(_include_electrons ? &getMaterialProperty<Real>("massem") : nullptr),
    _muem(_include_electrons ? &getADMaterialProperty<Real>("muem") : nullptr),
    _se_coeff(_include_secondary_electrons ? &getMaterialProperty<Real>("se_coeff") : nullptr),
    _e(_include_electrons ? &getMaterialProperty<Real>("e") : nullptr),
    _kb(_include_electrons ? &getMaterialProperty<Real>("k_boltz") : nullptr),
    _r_units(1. / getParam<Real>("position_units")),
    _r_ion(getParam<Real>("r_ion")),
    _r_electron(getParam<Real>("r_electron")),
    _potential_units(getParam<std::string>("potential_units")),
    _grad_potential(adCoupledGradient("potential"))
{
  if (_include_electrons and !_include_secondary_electrons)
    mooseError("SurfaceChargeNew: secondary_electrons is set to true, but include electrons is set "
               "to false.");

  // Define voltage scaling parameter
  if (_potential_units.compare("V") == 0)
    _voltage_scaling = 1.;
  else if (_potential_units.compare("kV") == 0)
    _voltage_scaling = 1000;

  _num_species = coupledComponents("species");

  // Resize the vectors to store _num_species pointers
  _species.resize(_num_species);
  _grad_species.resize(_num_species);
  _mu.resize(_num_species);
  _diff.resize(_num_species);
  _sgn.resize(_num_species);
  _mass.resize(_num_species);
  _Tion.resize(_num_species);

  for (unsigned int i = 0; i < _num_species; ++i)
  {
    _species[i] = &adCoupledValue("species", i);
    _grad_species[i] = &adCoupledGradient("species", i);
    _mu[i] = &getADMaterialProperty<Real>("mu" + (*getVar("species", i)).name());
    _diff[i] = &getADMaterialProperty<Real>("diff" + (*getVar("species", i)).name());
    _sgn[i] = &getMaterialProperty<Real>("sgn" + (*getVar("species", i)).name());
    _mass[i] = &getMaterialProperty<Real>("mass" + (*getVar("species", i)).name());
    _Tion[i] = &getADMaterialProperty<Real>("T" + (*getVar("species", i)).name());
  }

  // Precalculate constant values
  _q_times_NA = 1.602e-19 * 6.022e23 / _voltage_scaling;
}

void
SurfaceChargeNew::initQpStatefulProperties()
{
  _sigma[_qp] = 0;
}

void
SurfaceChargeNew::computeElectronFlux()
{
  // _normals * -1 * -_grad_potential -- negative signs cancel!
  if (_normals[_qp] * _grad_potential[_qp] > 0.0)
    _a = 1.0;
  else
    _a = 0.0;
  _ve_thermal =
      std::sqrt(8 * (*_e)[_qp] * 2.0 / 3 * std::exp((*_mean_energy)[_qp] - (*_electron)[_qp]) /
                (M_PI * (*_mass_em)[_qp]));

  _n_gamma = (1. - _a) * (*_se_coeff)[_qp] * _charge_flux /
             ((*_muem)[_qp] * -_grad_potential[_qp] * _r_units * _normals[_qp] +
              std::numeric_limits<double>::epsilon());

  _electron_flux = ((1. - _r_electron) / (1. + _r_electron) *
                        (-(2 * _a - 1) * (*_muem)[_qp] * -_grad_potential[_qp] * _normals[_qp] *
                             std::exp((*_electron)[_qp]) +
                         (0.5 * _ve_thermal * (std::exp((*_electron)[_qp]) - _n_gamma))) -
                    2. / (1. + _r_electron) * (1. - _a) * (*_se_coeff)[_qp] * _charge_flux);
}

void
SurfaceChargeNew::computeQpProperties()
{
  if (_material_data_type == Moose::FACE_MATERIAL_DATA || boundaryRestricted())
  {

    // COMPUTE CHARGED FLUX
    _charge_flux = 0.0;
    for (unsigned int i = 0; i < _num_species; ++i)
    {
      if (_normals[_qp] * (*_sgn[i])[_qp] * -_grad_potential[_qp] > 0.0)
        _b = 1.0;
      else
        _b = 0.0;
      //_charge_flux += 0;
      _charge_flux +=
          (*_sgn[i])[_qp] * std::exp((*_species[i])[_qp]) *
          (0.5 * std::sqrt(8 * (*_kb)[_qp] * (*_Tion[i])[_qp] / (M_PI * (*_mass[i])[_qp])) +
           (2 * _b - 1) * (*_sgn[i])[_qp] * (*_mu[i])[_qp] * -_grad_potential[_qp] * _r_units *
               _normals[_qp]);
      /*
      _charge_flux += (*_sgn[i])[_qp] * std::exp((*_species[i])[_qp]) *
                      ((*_sgn[i])[_qp] * (*_mu[i])[_qp] * -_grad_potential[_qp] * _r_units -
                       (*_diff[i])[_qp] * (*_grad_species[i])[_qp] * _r_units) *
                      _normals[_qp];
                      */
    }
    if (_include_electrons)
    {
      computeElectronFlux();
      //std::cout << MetaPhysicL::raw_value(_charge_flux) << ", "
      //          << MetaPhysicL::raw_value(_electron_flux) << std::endl;
      _charge_flux -= _electron_flux;
    }

    _plasma_current[_qp] = MetaPhysicL::raw_value(_charge_flux) * 1.602e-19 * 6.022e23;

    _sigma[_qp] = _sigma_old[_qp] + _charge_flux * _q_times_NA * _dt;
    //_sigma[_qp] = 0;
  }
  else
    _sigma[_qp] = 0.;
}
