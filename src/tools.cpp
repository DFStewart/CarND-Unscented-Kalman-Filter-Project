#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		std::cout << "Invalid estimation or ground_truth data" << std::endl;
		return rmse;
	}
	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

float Tools::CalculateNIS(VectorXd &z,
                          VectorXd &z_pred,
                          MatrixXd &S_inverse)
{
	VectorXd z_diff = z - z_pred;
	float NIS       = z_diff.transpose() * S_inverse * z_diff;
	//return the result
	return NIS;
}

float Tools::WraptoPi(float angle)
{
	float ang_wrap = angle;

  if (ang_wrap>M_PI)
  {
    while (ang_wrap > M_PI)
    {
    	ang_wrap = ang_wrap - 2.0*M_PI;
    }
  }
  else if (ang_wrap < -M_PI)
  {
    while (ang_wrap < -M_PI)
    {
    	ang_wrap = ang_wrap + 2.0*M_PI;
    }
  }
  return ang_wrap;
}
