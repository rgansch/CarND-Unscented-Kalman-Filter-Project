#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include "tools.h"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool laser_enabled_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool radar_enabled_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;
	
	///* laser measurement noise matrix
  MatrixXd R_laser_;
	
	///* radar measurement noise matrix
  MatrixXd R_radar_;
	
	///* measurement covariance matrix
  MatrixXd S_;	

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

	///* sigma points
  MatrixXd Xsig_;
	
	///* augmented sigma points
  MatrixXd Xsig_aug_;
	
	///* sigma points in measurement space
	MatrixXd Zsig_;
	
	///* measurement vector
	VectorXd z_;
	
	///* prediticed measurement vector
	VectorXd z_pred_;
	
  ///* time when the state is true, in us
  long long time_us_;

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
	
	///* Number of sigma points
	int n_sig_;

  ///* Sigma point spreading parameter
  double lambda_;

	///* NIS for laser
	double NIS_laser_;

	///* NIS for radar
	double NIS_radar_;	
	
	///* helper functions
	Tools tools_;

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
   * AugmentedSigmaPoints
   */
	void AugmentedSigmaPoints();
	 
	 /**
   * SigmaPointPrediction
   */
	void SigmaPointPrediction(double delta_t);
	
	/**
   * PredictMeanAndCovariance
   */
	void PredictMeanAndCovariance();
	
  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

	/**
   * PredictRadarMeasurement
   */
	void PredictRadarMeasurement();
	
	/**
   * PredictLaserMeasurement
   */	
	void PredictLaserMeasurement();
	
  /**
   * UpdateLaser
   * @param meas_package The measurement at k+1
   */
  void UpdateLaser(MeasurementPackage meas_package);

  /**
   * UpdateRadar
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);
	
	/**
	* UpdateState
	* @param n_z Dimension for measurement state (laser=2, radar=3)
	*/
	void UpdateState(int n_z);
};

#endif /* UKF_H */
