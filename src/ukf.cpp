#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_      = 3.0; // maximum acceleration of a bike is about 1/4 the sports car example 12 m/s/s from Leture slides

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_  = 0.7;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_  = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_  = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_   = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_  = 0.3;

  /**
  Complete the initialization. See ukf.h for other member properties.
  */

  // Initialize the filter
  is_initialized_ = false;

  // previous timestamp
  prev_time_ = 0;

  // State dimension
  n_x_   = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Lambda parameters
  lambda_ = 3.0 - n_aug_;

  // Weights
  weights_ = VectorXd(2*n_aug_+1);
  weights_.fill(0.0);
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i=1; i<2*n_aug_+1; i++)
  {
    weights_(i) = 0.5/(lambda_+n_aug_);
  }

  // Initialize State Vector
  x_ = VectorXd(n_x_);

  // Initialize Covariance Matrix
  P_ = MatrixXd(n_x_, n_x_);

  std::cout << "Initialized Variables" << endl;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  Complete this function! Make sure you switch between LIDAR and radar
  measurements.
  */
  // Initialize the filter on first pass
  if(!is_initialized_)
  {
	    // Initialize State Vector for RADAR
	    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
	    {
	      float rho    = meas_package.raw_measurements_(0);
	      float phi    = meas_package.raw_measurements_(1);
	      float rhodot = meas_package.raw_measurements_(2);
	      x_  <<  rho * cos(phi),
	              rho * sin(phi),
	              0.0,
	              0.0,
	              0.0;
	    }
	    // Initialize State Vector for LIDAR
	    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
	    {
	      float px = meas_package.raw_measurements_(0);
	      float py = meas_package.raw_measurements_(1);
	      x_  <<  meas_package.raw_measurements_(0),
	              meas_package.raw_measurements_(1),
	              0.0,
	              0.0,
	              0.0;
	    }

	    // Initialize Covariance Matrix
	    P_ << std_radr_*std_radr_, 0.0,                 0.0,           0.0,                     0.0,
	            0.0,               std_radr_*std_radr_, 0.0,           0.0,                     0.0,
	            0.0,               0.0,                 std_a_*std_a_, 0.0,                     0.0,
	            0.0,               0.0,                 0.0,           std_radphi_*std_radphi_, 0.0,
	            0.0,               0.0,                 0.0,           0.0,                     1.0;

	    prev_time_      = meas_package.timestamp_;
	    is_initialized_ = true;

	    return;
  }

  /*****************************************************************************
   *  Functionality to remove LIDAR or RADAR
   ****************************************************************************/

  // skip measurement if not supposed to use this type of sensor
  if (!use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    cout << "RADAR Not Enabled" << endl;
    return;
  }
  if (!use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
	cout << "LIDAR Not Enabled" << endl;
    return;
  }

  /*****************************************************************************
   *  Compute delta time
   ****************************************************************************/

  double dt   = (meas_package.timestamp_ - prev_time_) / 1e6;// usec to sec
  cout << "Timestamp          = " << meas_package.timestamp_ << endl;
  cout << "Previous Timestamp = " << prev_time_ << endl;
  cout << "dt                 = " << dt << endl;
  prev_time_  = meas_package.timestamp_;
  // dt is not constant in the dataset (sometimes larger than 0.1s)
  //perform more predictions to smooth data
  while (dt > 0.2)
  {
    Prediction(0.2);
    dt = dt - 0.2;
  }

  /*****************************************************************************
   *  Time Update / Prediction Update
   ****************************************************************************/
  std::cout << "Start Time Update\n";
  Prediction(dt);
  std::cout << "Finish Time Update\n";
  /*****************************************************************************
   *  Measurement Update
   ****************************************************************************/
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
	std::cout << "Start RADAR Update\n";
    UpdateRadar(meas_package);
    std::cout << "Finish RADAR Update\n";
  }
  else
  {
	std::cout << "Start LIDAR Update\n";
    UpdateLidar(meas_package);
    std::cout << "Finish LIDAR Update\n";
  }
  std::cout << "Completed Prediction + Measurement Update\n";
  std::cout << "\\-----------------------------------------\\" << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t Time between k and k+1 in s
 */
void UKF::Prediction(double delta_t) {
  /**

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  std::cout << "Augment Sigma Points\n";
  AugmentedSigmaPoints(&Xsig_aug);
  std::cout << "Sigma Point Prediction\n";
  SigmaPointPrediction(Xsig_aug, delta_t, &Xsig_pred_); // use our instance variable to later reference it in measurement update

  // Predict the Mean and Covariance
  std::cout << "Predict Mean and Covariance\n";
  PredictMeanAndCovariance(&x_, &P_);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param meas_package The measurement at k+1
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int n_z         = 2;
  MatrixXd Zsig   = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S      = MatrixXd(n_z,n_z);

  // Predict LIDAR measurement: Calculate Zsig, z_pred and S

  //transform sigma points into LIDAR measurement space
  for (int i=0; i<2*n_aug_+1; i++)
  {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    Zsig(0,i)  = p_x;
    Zsig(1,i)  = p_y;
  }

  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0.,
          0.,                    std_laspy_*std_laspy_;

  PredictMeas(n_z, Zsig, R, &z_pred, &S);

  VectorXd z = meas_package.raw_measurements_;
  UpdateState(n_z, Zsig, z_pred, S, z);

  // calculate NIS
  MatrixXd Sinv = S.inverse();
  NIS_laser_ = tools.CalculateNIS(z, z_pred, Sinv);
  cout << "NIS LIDAR          = " << NIS_laser_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param meas_package The measurement at k+1
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int n_z = 3;

  MatrixXd Zsig   = MatrixXd(n_z, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z);
  MatrixXd S      = MatrixXd(n_z,n_z);

  // Predict RADAR Measurement:
  //    Transform sigma points into RADAR measurement space
  for (int i=0; i<2*n_aug_+1; i++)
  {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double v1  = cos(yaw)*v;
    double v2  = sin(yaw)*v;
    Zsig(0,i)  = sqrt(p_x*p_x+p_y*p_y); //rho
    Zsig(1,i)  = atan2(p_y, p_x);       //phi
    //Prevent Divide by Zero
    if(Zsig(0,i) == 0.0)
    {
    	std::cout << "RADAR Zsig(0,i) == 0.0" << std::endl;
    	return;
    }
    else
    	Zsig(2,i)  = (p_x*v1 + p_y*v2) / Zsig(0,i);
  }

  // RADAR measurement noise matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0.,                      0.,
          0.,                  std_radphi_*std_radphi_, 0.,
          0.,                  0.,                      std_radrd_*std_radrd_;

  PredictMeas(n_z, Zsig, R, &z_pred, &S);

  VectorXd z = meas_package.raw_measurements_;
  UpdateState(n_z, Zsig, z_pred, S, z);

  // Calculate NIS
  MatrixXd Sinv = S.inverse();
  NIS_radar_ =  tools.CalculateNIS(z, z_pred, Sinv);
  cout << "NIS RADAR          = " << NIS_radar_ << endl;
}

/**
 * Computes the augmented sigma points
 * @param Xsig_out  augmented state
 */
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out)
{
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0.;
  x_aug(6) = 0.;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_,2*n_aug_+1);

  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1)          = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug(3,i+1)              = tools.WraptoPi(Xsig_aug(3,i+1));
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug(3,i+1+n_aug_)       = tools.WraptoPi(Xsig_aug(3,i+1+n_aug_));
  }

  *Xsig_out = Xsig_aug;
}

/**
 * Inserts every augmented sigma point into process model
 * @param Xsig_    aug augmented state vector input from AugmentedSigmaPoints()
 * @param delta_t  delta time
 * @param Xsig_out state vector output
 */
void UKF::SigmaPointPrediction(MatrixXd &Xsig_aug, double delta_t, MatrixXd* Xsig_out)
{
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //predict sigma points
  for (int i=0; i<2*n_aug_+1; i++)
  {
    float p_x       = Xsig_aug(0,i);
    float p_y       = Xsig_aug(1,i);
    float v         = Xsig_aug(2,i);
    float yaw       = Xsig_aug(3,i);
    float yawd      = Xsig_aug(4,i);
    float nu_a      = Xsig_aug(5,i);
    float nu_yawdd  = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd)>0.001)
    {
    	px_p = p_x + v/yawd * ( sin(yaw+yawd*delta_t) - sin(yaw)) ;
    	py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t));
    }
    else
    {
      // no change in yaw, going on straight line
    	px_p = p_x + v*delta_t * cos(yaw);
    	px_p = p_y + v*delta_t * sin(yaw);
    }

    double v_p    = v;
    double yaw_p  = yaw + yawd*delta_t;
    double yawd_p = yawd;

    // add noise
    px_p   = px_p + 0.5*nu_a*delta_t*delta_t*cos(yaw);
    py_p   = py_p + 0.5*nu_a*delta_t*delta_t*sin(yaw);
    v_p    = v + nu_a*delta_t;
    yaw_p  = yaw_p  + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  //write result
  *Xsig_out = Xsig_pred;
}

/**
 * Compute the predicted state mean and covariance using the sigma points
 * @param return_x predicted state mean output
 * @param return_P predicted covariance matrix output
 */
void UKF::PredictMeanAndCovariance(VectorXd* return_x, MatrixXd* return_P)
{
  //Note: Weights are set on initialization

  VectorXd x = VectorXd(n_x_);      //predicted state
  MatrixXd P = MatrixXd(n_x_, n_x_);//covariance matrix for prediction

  //predict state mean
  x.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++) // iterate over sigma points
  {
    x = x + weights_(i)*Xsig_pred_.col(i);
  }

  std::cout << "Pred x = " << x << endl;

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i=1; i<2*n_aug_+1; i++) // iterate over sigma points
  {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    x_diff(3) = tools.WraptoPi(x_diff(3));
    P         = P + weights_(i)*x_diff*x_diff.transpose();
  }

  std::cout << "Pred P = " << P << endl;

  *return_P = P;
  *return_x = x;
}

/**
 * Compute the mean measurement state vector z and measurement covariance matrix s
 * @param n_z      measurement state size
 * @param Zsig     Z sigma points state vector input
 * @param R        measurement noise matrix input
 * @param return_z mean measurement state vector "z" output
 * @param return_S measurement covariance matrix "S" output
 */
void UKF::PredictMeas(int n_z, MatrixXd &Zsig, MatrixXd &R, VectorXd* return_z, MatrixXd* return_S)
{
  //Compute mean measurement state vector z
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++) //2n+1 sigma points
  {
    z_pred = z_pred + weights_(i)*Zsig.col(i);
  }

  //Compute measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++) //2n+1 sigma points
  {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (n_z > 2) // RADAR has more states than the LIDAR
    {
      z_diff(1) = tools.WraptoPi(z_diff(1));
    }
    S = S + weights_(i)*z_diff*z_diff.transpose();
  }

  //Sum the measurement noise and measurement covariance matricies
  S = S + R;

  *return_z = z_pred;
  *return_S = S;
}

/**
 * Perform the measurement update and update the state and covariance update with the measurement
 * @param n_z    measurement state size
 * @param Zsig Z sigma points measurement state vector
 * @param z_pred predicted measurement state vector
 * @param S      measurement covariance matrix "S" output
 * @param z      raw measurement vector
 */
void UKF::UpdateState(int n_z, MatrixXd &Zsig, VectorXd &z_pred, MatrixXd &S, VectorXd &z)
{
  MatrixXd Tc = MatrixXd(n_x_, n_z);//Cross Correlation Tc

  //Calculate cross correlation matrix
  Tc.fill(0.0);

  for (int i=0; i<2*n_aug_+1; i++)//2n+1 sigma points
  {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    if (n_z == 3) // awkward check if its radar measurement
      z_diff(1) = tools.WraptoPi(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = tools.WraptoPi(x_diff(3));

    Tc = Tc + weights_(i)*x_diff*z_diff.transpose();
  }

  //Kalman gain k
  MatrixXd K = Tc*S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;
  if (n_z > 2) // RADAR has more states than the LIDAR
  {
    z_diff(1) = tools.WraptoPi(z_diff(1));
  }

  //update state mean and covariance matrix
  x_    = x_ + K*z_diff;
  P_    = P_ - K*S*K.transpose();
  std::cout << "x = " << x_ << endl;
}
