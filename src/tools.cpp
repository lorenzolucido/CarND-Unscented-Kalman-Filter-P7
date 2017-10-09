#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() == 0 || estimations.size() != ground_truth.size())
  {
      cout << "Invalid inputs" << endl;
  }
  //accumulate squared residuals
  VectorXd residuals;
  for(unsigned int i=0; i < estimations.size(); ++i){
    residuals = estimations[i] - ground_truth[i];
    residuals = residuals.array() * residuals.array();
    rmse += residuals;
  }
  rmse = rmse/estimations.size();
  rmse = rmse.array().sqrt();
  
  return rmse;
}