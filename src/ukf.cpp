#include "ukf.h"
#include "Eigen/Dense"
#include "json.hpp"
#include <iostream>

using json = nlohmann::json;
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
	// Load parameters from json file
	std::ostringstream param_buf; 
	std::ifstream param_file("../data/config.json"); 
	param_buf << param_file.rdbuf(); 
	auto param = json::parse(param_buf.str());
	param_file.close();
	
	// set initialization state
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  laser_enabled_ = (bool) param["UKF"]["laser_enabled"].get<int>();

  // if this is false, radar measurements will be ignored (except during init)
  radar_enabled_ = (bool) param["UKF"]["radar_enabled"].get<int>();

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
  std_a_ = param["UKF"]["std_a"].get<float>();

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = param["UKF"]["std_yawdd"].get<float>();

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = param["UKF"]["std_laspx"].get<float>();

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = param["UKF"]["std_laspy"].get<float>();

  // Radar measurement noise standard deviation radius in m
  std_radr_ = param["UKF"]["std_radr"].get<float>();

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = param["UKF"]["std_radphi"].get<float>();

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = param["UKF"]["std_radrd"].get<float>();

	// R matrices
  R_radar_ = MatrixXd(3,3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
						  0, std_radphi_*std_radphi_, 0,
							0, 0, std_radrd_*std_radrd_;
  
	R_laser_ = MatrixXd(2,2);
	R_laser_ << std_laspx_*std_laspx_, 0,
							0, std_laspy_*std_laspy_;
							
	// set the state dimension
	n_x_ = param["UKF"]["n_x"].get<int>();
	
	// set the augmented state dimension
  n_aug_ = param["UKF"]["n_aug"].get<int>();
	
	// Number of sigma points
  n_sig_ = 2 * n_aug_ + 1;
	
	// set the sigma point spreading parameter
	lambda_ = 3 - n_x_;
	
	// set weights
  weights_ = VectorXd(2*n_aug_ + 1);

  double weight_0 = lambda_ / (lambda_+n_aug_);
  weights_(0) = weight_0;
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
	std::cout << std::endl << "UKF::ProcessMeasurment " << std::endl;

	/** 
	==============
	INITIALIZATION
	==============
	*/
	if(!is_initialized_) {
		std::cout << "  Intialization" << std::endl;
		double px;
		double py;
		double vx;
		double vy;
		
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR && radar_enabled_) {
			std::cout<<"Initialize with RADAR"<<std::endl;

			double rho = meas_package.raw_measurements_[0];
			double phi = meas_package.raw_measurements_[1];
			double rhodot = meas_package.raw_measurements_[2];
			
			px = rho*cos(phi);
			py = rho*sin(phi);
			vx = rhodot*cos(phi);
			vy = rhodot*sin(phi);
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER && laser_enabled_) {
			std::cout<<"Initialize with LASER"<<std::endl;

			px = meas_package.raw_measurements_[0];
			py = meas_package.raw_measurements_[1];
			vx = 0;
			vy = 0;
		}
		else {
			return;
		}
		
		x_ << px, py, sqrt(pow(vx, 2) + pow(vy, 2)), 0, 0;		
		time_us_ = meas_package.timestamp_;
		
		is_initialized_ = true;
		return;
	}
		
	// Skip prediction and update if sensor is not enabled
	if ((meas_package.sensor_type_ == MeasurementPackage::RADAR && !radar_enabled_) ||
			(meas_package.sensor_type_ == MeasurementPackage::LASER && !laser_enabled_)) {
		return;
  }
			
	/** 
	==========
	PREDICTION
	==========
	*/
	std::cout << "  Prediction" << std::endl;
	// calculate delta_t since last measurement
	double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
	time_us_ = meas_package.timestamp_;

	Prediction(delta_t);
	
	/** 
	======
	UPDATE
	======
	*/
  std::cout << "  Update" << std::endl;
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		UpdateRadar(meas_package);
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
		UpdateLaser(meas_package);
	}
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
	std::cout << "UKF::Prediction" << std::endl;
	AugmentedSigmaPoints();
	SigmaPointPrediction(delta_t);
	PredictMeanAndCovariance();
}

/**
 * Creates the augmented sigma points
 */
void UKF::AugmentedSigmaPoints() {
	std::cout << "UKF::AugmentedSigmaPoints" << std::endl;
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create augmented mean state
  x_aug.head(5) = x_;
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
  Xsig_aug_.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug_.col(i+1)        = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
}

/**
 * Predicts the sigma points
 * @param delta_t
 */
void UKF::SigmaPointPrediction(double delta_t) {
	std::cout << "UKF::SigmaPointPrediction" << std::endl;
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
    if (fabs(yawd) > 0.001) {
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
	std::cout << "UKF::PredictMeanAndCovariance" << std::endl;
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);
	
  //predicted state mean
  x.fill(0.0);
  for (int i=0; i < 2*n_aug_ + 1; i++) {  //iterate over sigma points
    x = x+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    //angle normalization
    x_diff(3) = tools_.WrapAngle(x_diff(3));

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }

  //write result
  x_ = x;
  P_ = P;
}

/**
 * Predict output state and measurement covariance for radar
 */
void UKF::PredictRadarMeasurement() {
	std::cout << "UKF::PredictRadarMeasurement" << std::endl;
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //matrix for sigma points in measurement space
  Zsig_ = MatrixXd(n_z, n_sig_);

  //transform sigma points into measurement space
  for (int i = 0; i < 2*n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig_(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  z_pred_ = VectorXd(n_z);
  z_pred_.fill(0.0);
  for (int i=0; i < 2*n_aug_ + 1; i++) {
      z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }
	
  // predicted measurement covariance
  S_ = MatrixXd(n_z, n_z);
  S_.fill(0.0);
  for (int i = 0; i < 2*n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    //angle normalization
    z_diff(1) = tools_.WrapAngle(z_diff(1));

    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }
  S_ = S_ + R_radar_;
}

/**
 * Predict output state and measurement covariance for laser
 */
void UKF::PredictLaserMeasurement() {
	std::cout << "UKF::PredictLaserMeasurment" << std::endl;
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
	std::cout << "UKF::UpdateLaser" << std::endl;
  z_ = VectorXd(2);
  z_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
  PredictLaserMeasurement();
  UpdateState(2);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
	std::cout << "UKF::UpdateRadar" << std::endl;
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
	std::cout << "UKF::UpdateState" << std::endl;
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
		
    //angle normalization
    z_diff(1) = tools_.WrapAngle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    x_diff(3) = tools_.WrapAngle(x_diff(3));
		
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S_.inverse();

  //residual
  VectorXd z_diff = z_ - z_pred_;

  //angle normalization
  z_diff(1) = tools_.WrapAngle(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S_*K.transpose();

  // Calculate NIS
  double NIS = z_diff.transpose() * S_.inverse() * z_diff;
	// n_z = 3 -> Radar, n_z = 2 -> Laser
  if (n_z == 3) {
    NIS_radar_ = NIS;
  }
  else {
    NIS_laser_ = NIS;
  }
}
