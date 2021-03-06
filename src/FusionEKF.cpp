#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0, 0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0, 0, 0.0009, 0, 0, 0, 0.09;

  // H matrix laser - current state into the measurement space of the sensor
  H_laser_ << 1, 0, 0, 0, 0, 1, 0, 0;

  // Hj_ is not initialized since it is calculated explicitly before first time
  // usage

  // initial transition matrix F_ (assuming a dt of 1s) (also actually not
  // necessary)
  ekf_.F_ = MatrixXd::Identity(4, 4);
  float dt = 1;
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // initial state covariance matrix P
  ekf_.P_ = MatrixXd::Identity(4, 4);
  ekf_.P_(2, 2) = 1000;
  ekf_.P_(3, 3) = 1000;
}

FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // converting radar measurements (polar coordinates) to carthesian
      // coordinates
      float rho = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float rho_dot = measurement_pack.raw_measurements_(2);
      ekf_.x_(0) = rho * cos(phi);
      ekf_.x_(1) = rho * sin(phi);
      ekf_.x_(2) = rho_dot * cos(phi);
      ekf_.x_(3) = rho_dot * sin(phi);
    }

    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_(0) = measurement_pack.raw_measurements_(0);
      ekf_.x_(1) = measurement_pack.raw_measurements_(1);
      ekf_.x_(2) = 0.f;
      ekf_.x_(3) = 0.f;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // Time keeping: Time from ms to s -> division by 1e6
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  // Update the state transition matrix F according to the new elapsed time
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // Update the process noise covariance matrix.
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  // set the process covariance matrix Q
  // using noise_ax = 9 and noise_ay = 9 for your Q matrix.
  float noise_ax = 9.f;
  float noise_ay = 9.f;

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0, 0,
      dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay, dt_3 / 2 * noise_ax, 0,
      dt_2 * noise_ax, 0, 0, dt_3 / 2 * noise_ay, 0, dt_2 * noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates

    // H depends on current state -> CalculateJacobian
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    // using radar specific measurement noise matrix R_radar_
    ekf_.R_ = R_radar_;

    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates

    // H is always linear
    ekf_.H_ = H_laser_;

    // using laser specific measurement noise matrix R_laser_
    ekf_.R_ = R_laser_;

    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "State vector x_\n";
  cout << "px = " << ekf_.x_(0) << endl;
  cout << "py = " << ekf_.x_(1) << endl;
  cout << "vx = " << ekf_.x_(2) << endl;
  cout << "vy = " << ekf_.x_(3) << endl << endl;
  cout << "State covariance matrix P_\n";
  cout << ekf_.P_ << endl << endl << endl;
}
