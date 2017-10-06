#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

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
  
  // set initialization state
  is_initialized_ = false;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
				
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;

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

	// R matrices
  R_radar_ = MatrixXd(3,3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
						  0, std_radphi_*std_radphi_, 0,
							0, 0, std_radrd_*std_radrd_;
  
	R_laser_ = MatrixXd(2,2);
	R_laser_ << std_laspx_*std_laspx_, 0,
							0, std_laspy_*std_laspy_;
							
	// set the state dimension
	n_x_ = 3
	
	// set the augmented state dimension
  n_aug_ = 7;
	
	// set the sigma point spreading parameter
	lambda_ = 3 - n_x_;
	
	// set weights
  weights_ = VectorXd(2*n_aug_ + 1);

  double weight_0 = lambda_ / (lambda_+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i < 2*n_aug_ + 1; i++) {
    double weight = 0.5 / (n_aug_+lambda_);
    weights_(i) = weight;
  }
	
	// set sigma point matrix
  Xsig_ = MatrixXd(n_x_, 2 * n_x_ + 1);

  // set augmented sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, n_sig_);

  // set predicted sigma points
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
	// first measurement
	cout << "UKF: " << endl;

	/** 
	==============
	INITIALIZATION
	==============
	*/
	if(!is_initialized_) {
		double px;
		double py;
		double vx;
		double vy;
		
		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
			std::cout<<"Initialize with RADAR"<<std::endl;

			double rho = measurement_pack.raw_measurements_[0];
			double phi = measurement_pack.raw_measurements_[1];
			double rhodot = measurement_pack.raw_measurements_[2];
			
			px = rho*cos(phi);
			py = rho*sin(phi);
			vx = rhodot*cos(phi);
			vy = rhodot*sin(phi);
		}
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			std::cout<<"Initialize with LASER"<<std::endl;

			px = measurement_pack.raw_measurements_[0];
			py = measurement_pack.raw_measurements_[1];
			vx = 0;
			vy = 0;
		}
		
		x_ << px, py, sqrt(pow(vx, 2) + pow(vy, 2)), 0, 0;		
		time_us_ = measurement_pack.timestamp_;
		
		is_initialized_ = true;
		return;
	}
		
	/** 
	==========
	PREDICTION
	==========
	*/
	// calculate delta_t since last measurement
	double delta_t = (measurement_pack.timestamp_ - time_us_) / 1000000.0;
	time_us_ = measurement_pack.timestamp_;

	Prediction(delta_t);
	
	/** 
	======
	UPDATE
	======
	*/
	if (meas_package.sensor_type_ == measurement_package::RADAR) {
		UpdateRadar(meas_package);
	}
	else {
		UpdateLidar(meas_package);
	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	AugmentedSigmaPoints();
	SigmaPointPrediction(delta_t);
	PredictMeanAndCovariance();
}

/**
 * Creates the augmented sigma points
 */
void UKF::AugmentedSigmaPoints() {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create augmented mean state
  x_aug.head(5) = x;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug; i++)
  {
    Xsig_aug_.col(i+1)       = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
}

/**
 * Predicts the sigma points
 * @param delta_t
 */
void UKF::SigmaPointPrediction(double delta_t) {
  //predict sigma points
  for (int i = 0; i< 2*n_aug_ + 1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd_) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

/**
 * Predicts the mean state and covariance matrix
 */
void UKF::PredictMeanAndCovariance() {
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);
	
  //predicted state mean
  x.fill(0.0);
  for (int i=0; i < 2*n_aug_ + 1; i++) {  //iterate over sigma points
    x = x+ weights(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    x_diff(3) = tools.WrapAngle(x_diff(3));

    P = P + weights(i) * x_diff * x_diff.transpose() ;
  }

  //write result
  x_ = x;
  P_ = P;
}

/**
 * Predict output state and measurement covariance for radar
 */
void UKF::PredictRadarMeasurement() {
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //matrix for sigma points in measurement space
  Zsig_ = MatrixXd(n_z, n_sig_);

  //transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred(0,i);
    double p_y = Xsig_pred(1,i);
    double v  = Xsig_pred(2,i);
    double yaw = Xsig_pred(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug+1; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }
	
  // predicted measurement covariance
  S_ = MatrixXd(n_z, n_z);
  S_.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    z_diff(1) = tools.wrapAngle(z_diff(1));

    S_ = S_ + weights(i) * z_diff * z_diff.transpose();
  }
  S_ = S_ + R_radar_;
}

/**
 * Predict output state and measurement covariance for laser
 */
void UKF::PredictLaserMeasurement() {
  //set measurement dimension
  int n_z = 2;

  // mean predicted measurement
  z_pred_ = VectorXd(n_z);

  //matrix for sigma points in measurement space
  Zsig_ = MatrixXd(n_z, n_sig_);

  //transform sigma points into measurement space
  for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    /*if (fabs(p_x) < p_x_min_ || fabs(p_y)< p_y_min_) {
      p_x = p_x_min_;
      p_y = p_y_min_;
    }*/

    // measurement model
    Zsig_(0,i) = p_x;
    Zsig_(1,i) = p_y;
  }

  //mean predicted measurement
  z_pred_.fill(0.0);
  for (int i=0; i < n_sig_; i++) {
      z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  // predicted measurement covariance
  S_ = MatrixXd(n_z, n_z);
  S_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S_ = S_ + R_laser_;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLaser(MeasurementPackage meas_package) {
  z_ = VectorXd(2);
  z_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
  PredictLidarMeasurement();
  UpdateState(2);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  z_ = VectorXd(3);
  z_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];
  PredictRadarMeasurement();
  UpdateState(3);
}

/**
 * Common code for update of laser and radar
 * @param n_z Dimension for measurement state (laser=2, radar=3)
 */
void UKF::UpdateState(int n_z) {

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
		
    //angle normalization
    z_diff(1) = tools.WrapAngle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    x_diff(3) = tools.WrapAngle(x_diff(3));
		
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S_.inverse();

  //residual
  VectorXd z_diff = z_ - z_pred_;

  //angle normalization
  z_diff(1) = tools.WrapAngle(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S_*K.transpose();

  // Calculate NIS
  double NIS = z_diff.transpose() * S_.inverse() * z_diff;
  if (is_radar == true) {
    NIS_radar_ = NIS;
  }
  else {
    NIS_laser_ = NIS;
  }
}
