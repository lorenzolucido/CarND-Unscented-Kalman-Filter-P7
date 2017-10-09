#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

# define M_PI           3.14159265358979323846  /* pi */

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  ///* State dimension
  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Zero(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.25;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.25;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  ///* initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;
  
  ///* time when the state is true, in us
  time_us_ = 0;

  ///* Augmented state dimension
  n_aug_ = 7;

  ///* Number of generated sigma points
  n_sigma_ = 2 * n_aug_ + 1;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  ///* Weights of sigma points
  weights_ = VectorXd(n_sigma_);
  weights_.fill(1 / (2*lambda_ + 2*n_aug_));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  ///* predicted sigma points matrix
  Xsig_pred_ = MatrixXd::Zero(n_x_, n_sigma_);
}

UKF::~UKF() {}



/**
 * Generates augmented sigma points
 */
MatrixXd UKF::AugmentedSigmaPoints() {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd::Zero(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd::Zero(n_aug_, n_sigma_);

  //create augmented mean state
  x_aug << x_, 0., 0.;
  //create augmented covariance matrix
  MatrixXd Q(2,2); 
  Q <<  std_a_ * std_a_,  0, 
        0,                std_yawdd_ * std_yawdd_;
  

  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(n_aug_-n_x_, n_aug_-n_x_) = Q;


  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL(); 
  MatrixXd A_norm = A * sqrt(lambda_ + n_aug_);

  //create augmented sigma points
  Xsig_aug.colwise() += x_aug;
  Xsig_aug.block<7,7>(0,1) += A_norm;
  Xsig_aug.block<7,7>(0,n_aug_+1) -= A_norm;

  return Xsig_aug;
}

/**
 * Maps a sigma point to a predicted sigma point
 */
VectorXd UKF::SigmaPointPredict(VectorXd sigmaPoint, double delta_t)
{
    double px = sigmaPoint(0);
    double py = sigmaPoint(1);
    double v = sigmaPoint(2);
    double phi = sigmaPoint(3);
    double phidot = sigmaPoint(4);
    double nu_a = sigmaPoint(5);
    double nu_phidot = sigmaPoint(6);
    
    
    VectorXd first = VectorXd(5);
    if(fabs(phidot) > 0.001)
        first <<    (v / phidot) * (sin(phi + phidot * delta_t) - sin(phi)), 
                    (v / phidot) * (-cos(phi + phidot * delta_t) + cos(phi)),
                    0., 
                    phidot * delta_t, 
                    0.;
    else 
        first <<    v * cos(phi) * delta_t,
                    v * sin(phi) * delta_t,
                    0., 
                    phidot * delta_t, 
                    0.;
    
    VectorXd second = VectorXd(5);
    second <<   0.5 * delta_t * delta_t * cos(phi) * nu_a,
                0.5 * delta_t * delta_t * sin(phi) * nu_a,
                delta_t * nu_a,
                0.5 * delta_t * delta_t * nu_phidot,
                delta_t * nu_phidot;
    
    VectorXd output = sigmaPoint.head(5) + first + second;
    return output;
}

/**
 * Updates predicted sigma points
 */
void UKF::UpdatePredictedSigmaPoints(double delta_t) {
  // Generate sigma points
  MatrixXd Xsig_aug = UKF::AugmentedSigmaPoints();
  
  // Create matrix with predicted sigma points as columns  
  for(int i = 0; i < Xsig_aug.cols(); i++)
    Xsig_pred_.col(i) = UKF::SigmaPointPredict(Xsig_aug.col(i), delta_t);
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  bool radar_sensor = meas_package.sensor_type_ == MeasurementPackage::RADAR;
  bool laser_sensor = meas_package.sensor_type_ == MeasurementPackage::LASER;


  if(!is_initialized_) {
    time_us_ = meas_package.timestamp_;

    if(radar_sensor and use_radar_) { // We have only two types of sensors
      double rho = meas_package.raw_measurements_(0); 
      double phi = meas_package.raw_measurements_(1);
      
      x_ << rho * cos(phi), rho * sin(phi), 0.0, 0.0, 0.0;
      is_initialized_ = true;
    }
    else if(laser_sensor and use_laser_) {
      double px = meas_package.raw_measurements_(0);
      double py = meas_package.raw_measurements_(1);

      x_ << px, py, 0.0, 0.0, 0.0;
      is_initialized_ = true;
    }
          
    P_ = MatrixXd::Identity(n_x_, n_x_);
    // Done initializing, no need to predict or update
    return;
  }

  // 1. Updating time and sigma point predictions
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  UKF::UpdatePredictedSigmaPoints(dt); // This mutates Xsig_pred_

  // 2. Predicting new x and P
  UKF::Prediction(); // This mutates x_ and P_

  // 3. Updating (depends on sensor)
  if(radar_sensor and use_radar_) UKF::UpdateRadar(meas_package);
  if(laser_sensor and use_laser_) UKF::UpdateLidar(meas_package);

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction() {  
  double pi = acos(-1);
  // Predict state mean

  x_ = Xsig_pred_ * weights_;

  // Predict state covariance matrix
  P_ = MatrixXd::Zero(n_x_, n_x_);

  for(int i = 0; i < Xsig_pred_.cols(); i++) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    while (x_diff(3) > pi)  x_diff(3) -= 2.*pi;
    while (x_diff(3) < -pi) x_diff(3) += 2.*pi;

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:


  You'll also need to calculate the lidar NIS.
  */
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 2;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, n_sigma_);
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);

  //matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

  //transform sigma points into measurement space
  for(int i = 0; i < Zsig.cols(); i++)
    Zsig.col(i) << Xsig_pred_.col(i)(0), Xsig_pred_.col(i)(1);

  //calculate mean predicted measurement
  z_pred = Zsig * weights_;

  //calculate measurement covariance matrix S
  S <<  pow(std_laspx_,2),   0,                    
        0,                  pow(std_laspy_, 2);  
  
  for(int i = 0; i < Zsig.cols(); i++) {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  
  //calculate cross correlation matrix
  for(int i = 0; i < Zsig.cols(); i++)
    Tc += (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose() * weights_(i);

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  VectorXd z = meas_package.raw_measurements_;
  x_ = x_ + K * (z - z_pred);
  P_ = P_ - K * S * K.transpose();
  
}

/**
 * Maps a state to a radar measurement
 * @param VectorXd state
 */
VectorXd UKF::StateToRadarMeasurement(VectorXd state) {
  double px = state(0);
  double py = state(1);
  double v = state(2);
  double phi = state(3);
  
  VectorXd output = VectorXd(3);
  double rho = sqrt(px*px + py*py);
  output << rho, 
            atan2(py, px), 
            fabs(rho) < 0.00001 ? 0.00001 : (px*cos(phi)*v + py*sin(phi)*v)/rho;
  return output;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:



  You'll also need to calculate the radar NIS.
  */

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd::Zero(n_z, n_sigma_);
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd::Zero(n_z,n_z);

  //matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd::Zero(n_x_, n_z);

  //transform sigma points into measurement space
  for(int i = 0; i < Zsig.cols(); i++)
    Zsig.col(i) = UKF::StateToRadarMeasurement(Xsig_pred_.col(i));

  //calculate mean predicted measurement
  z_pred = Zsig * weights_;

  //calculate measurement covariance matrix S
  S <<  pow(std_radr_,2),   0,                    0,
        0,                  pow(std_radphi_, 2),  0,
        0,                  0,                    pow(std_radrd_, 2);
  
  for(int i = 0; i < Zsig.cols(); i++) {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    
    //angle normalization
    double pi = acos(-1);
    while (z_diff(1) > pi)  z_diff(1) -= 2.* pi;
    while (z_diff(1) < -pi) z_diff(1) += 2.* pi;

    S += weights_(i) * z_diff * z_diff.transpose();
  }
  
  //calculate cross correlation matrix
  for(int i = 0; i < Zsig.cols(); i++)
    Tc += (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose() * weights_(i);

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  VectorXd z = meas_package.raw_measurements_;
  x_ = x_ + K * (z - z_pred);
  P_ = P_ - K * S * K.transpose();
  }
