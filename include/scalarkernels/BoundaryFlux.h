#pragma once                                                                   
                                                                               
#include "ODEKernel.h"                                                         
// #include "RateCoefficientProvider.h"                                        
                                                                               
class BoundaryFlux;                                                     
                                                                               
template <>                                                                    
InputParameters validParams<BoundaryFlux>();                            
                                                                               
class BoundaryFlux : public ODEKernel                                   
{                                                                              
public:                                                                        
  BoundaryFlux(const InputParameters & parameters);                     
                                                                               
protected:                                                                     
  virtual Real computeQpResidual() override;                                   
  virtual Real computeQpJacobian() override;                                   
  //virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;         

  const PostprocessorValue & _current;
  // const RateCoefficientProvider & _data;                                    
};                              
