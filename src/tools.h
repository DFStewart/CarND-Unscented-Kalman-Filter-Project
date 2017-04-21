#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);
  float CalculateNIS(Eigen::VectorXd &z, Eigen::VectorXd &z_pred, Eigen::MatrixXd &S_inverse);
  float WraptoPi(float angle);
};

#endif /* TOOLS_H_ */
