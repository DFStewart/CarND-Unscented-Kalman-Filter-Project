#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_;

  ///*Previous timestamp
  long long prev_time_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  float NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  Tools tools;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  /**
   * Perform the measurement update and update the state and covariance update with the measurement
   * @param n_z    measurement state size
   * @param Zsig Z sigma points measurement state vector
   * @param z_pred predicted measurement state vector
   * @param S      measurement covariance matrix "S" output
   * @param z      raw measurement vector
   */
  void UpdateState(int n_z, MatrixXd &Zsig, VectorXd &z_pred, MatrixXd &S, VectorXd &z);

  /**
   * Computes the augmented sigma points
   * @param Xsig_out  augmented state
   */
  void AugmentedSigmaPoints(MatrixXd* Xsig_out);

  /**
   * Compute the predicted state mean and covariance using the sigma points
   * @param return_x predicted state mean output
   * @param return_P predicted covariance matrix output
   */
  void PredictMeanAndCovariance(VectorXd* return_x, MatrixXd* return_P);

  /**
   * Inserts every augmented sigma point into process model
   * @param Xsig_    aug augmented state vector input from AugmentedSigmaPoints()
   * @param delta_t  delta time
   * @param Xsig_out state vector output
   */
  void SigmaPointPrediction(MatrixXd &Xsig_aug, double delta_t, MatrixXd* Xsig_out);

  /**
   * Compute the mean measurement state vector z and measurement covariance matrix s
   * @param n_z      measurement state size
   * @param Zsig     Z sigma points state vector input
   * @param R        measurement noise matrix input
   * @param return_z mean measurement state vector "z" output
   * @param return_S measurement covariance matrix "S" output
   */
  void PredictMeas(int n_z, MatrixXd &Zsig, MatrixXd &R, VectorXd* z_out, MatrixXd* S_out);

};

#endif /* UKF_H */
