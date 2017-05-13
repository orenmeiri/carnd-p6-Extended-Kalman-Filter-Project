#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;


  MatrixXd P;	// object covariance matrix
  P = MatrixXd(4, 4);
  P << 1000, 0, 0, 0,
          0, 1000, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;

  ekf_ = KalmanFilter();
  ekf_.P_ = P;

  MatrixXd F; // state transition matrix
  F = MatrixXd(4, 4);
  F << 1, 0, 1, 0,
          0, 1, 0, 1,
          0, 0, 1, 0,
          0, 0, 0, 1;
  ekf_.F_ = F;

  MatrixXd Q;	// process covariance matrix
  Q = MatrixXd(4, 4);
//    Q << 9, 0, 0, 0,
//    0, 9, 0, 0,
//    0, 0, 0, 0,
//    0, 0, 0, 0;
  ekf_.Q_ = Q;

  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
          0, 1, 0, 0;

  ekf_.x_ = VectorXd(4);
  ekf_.x_ << 1, 1, 0, 0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      ekf_.x_ << cos(measurement_pack.raw_measurements_[1]) * measurement_pack.raw_measurements_[0],
      sin(measurement_pack.raw_measurements_[1]) * measurement_pack.raw_measurements_[0], 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        cout << "initial measurement: " << endl << measurement_pack.raw_measurements_ << endl;
        ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
        cout << "initial x: " << endl << ekf_.x_ << endl;
        cout << "initial P: " << endl << ekf_.P_ << endl;
      /**
      Initialize state.
      */
    }
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;
    int noise_ax = 9;
    int noise_ay = 9;
    ekf_.Q_ <<  dt_4/(4*noise_ax), 0, dt_3/(2*noise_ax), 0,
            0, dt_4/(4*noise_ay), 0, dt_3/(2*noise_ay),
            dt_3/(2*noise_ax), 0, dt_2*noise_ax, 0,
            0, dt_3/(2*noise_ay), 0, dt_2*noise_ay;


    ekf_.Predict();

    /*****************************************************************************
     *  Update
     ****************************************************************************/
  
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
      Hj_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.H_ = Hj_;
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
      ekf_.H_ = H_laser_;
      ekf_.R_ = R_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
//  cout << "x_ = " << ekf_.x_ << endl;
//  cout << "P_ = " << ekf_.P_ << endl;
}
