#include "CSVPotentialBC.h"

registerMooseObject("ZapdosApp", CSVPotentialBC);

InputParameters
CSVPotentialBC::validParams()
{
  InputParameters params = NodalBC::validParams();
  params.addRequiredParam<FileName>("file_name", "The file name.");
  params.declareControllable("file_name");
  return params;
}

CSVPotentialBC::CSVPotentialBC(const InputParameters & parameters) : NodalBC(parameters)
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
    mooseError("Unable to open file " + file_name + " for CSVPotentialBC.");

  _voltage.setData(x_val, y_val);
}

Real
CSVPotentialBC::computeQpResidual()
{
  std::cout << _u[_qp] - _voltage.sample(_t) << std::endl;
  return _u[_qp] - _voltage.sample(_t);
}

/*
Real
CSVPotentialBC::computeQpJacobian()
{
  return 0.0;
}
*/

/*
Real
CSVPotentialBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _v_num)
    return 2. * _v[_qp];
  else
    return 0.;
}
*/
