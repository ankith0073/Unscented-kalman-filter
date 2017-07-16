#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "defines.h"

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

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 10 * (M_PI/180);
  //  std_yawdd_ = 0.3;

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
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {


#ifdef prints
        cout << "Inside process measurement function" << endl;
#endif
        if (is_initialized_ == false) {

            unsigned int i = 0;
            if (meas_package.sensor_type_ == 0 and use_laser_) {
                x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0, 0, 0;
                P_ = MatrixXd::Identity(5, 5);
                P_(1, 1) = 0.15 * 0.15;
                P_(0, 0) = 0.15 * 0.15;
                is_initialized_ = true;
            } else if(meas_package.sensor_type_ == 1 and use_radar_){
                float p_x = meas_package.raw_measurements_(0) * cos(meas_package.raw_measurements_(1));
                float p_y = meas_package.raw_measurements_(0) * sin(meas_package.raw_measurements_(1));

                float vx = meas_package.raw_measurements_(2) * cos(meas_package.raw_measurements_(1));
                float vy = meas_package.raw_measurements_(2) * sin(meas_package.raw_measurements_(1));
                x_ << p_x, p_y , sqrt(vx*vx+vy*vy),0,0;
                P_ = MatrixXd::Identity(5, 5);
                is_initialized_ = true;
            }
            time_us_ = meas_package.timestamp_;

            //covariance of the process noise initialization
            Q_ = MatrixXd(number_indp_noise_processes_, number_indp_noise_processes_) * 0;

            Q_ << std_a_ * std_a_, 0,
                    0, std_yawdd_ * std_yawdd_;

            /* State dimensions and the dimensions of the augmented states remain constant throughout UKF processing*/
            n_x_ = x_.size();
            n_aug_ = n_x_ + number_indp_noise_processes_;
            lambda_ = 3 - n_aug_;
            //Set the lambda parameter of UKF
            weights_ = VectorXd(2 * n_aug_ + 1);
            for (unsigned int i = 0; i < n_aug_ + 1; i++) {
                if (i == 0) {
                    weights_(0) = lambda_ / (lambda_ + n_aug_);
                } else {
                    weights_(i) = 1 / (2 * (lambda_ + n_aug_));
                    weights_(i + n_aug_) = 1 / (2 * (lambda_ + n_aug_));
                }

            }
#ifdef prints
            cout << weights_ << endl;
#endif

#ifdef print_iterations_
            iterations_ = 1;
            cout << iterations_ << "\n\n" << endl;
#endif

#ifdef NIS_test
            NIS_lidar_text_.open("LIDAR_NIS.txt", ios::out);
            NIS_radar_text_.open("RADAR_NIS.txt", ios::out);
#endif

            //update the weights which remain constant through the UKF procedure

        }


  else
  {
    float dt = (meas_package.timestamp_ - time_us_)/1000000.0;
    time_us_ = meas_package.timestamp_;

    /*Generate sigma points of X_(k|k) by appending the noise processes into state vector*/
    /*sigma point matrix of size (n_aug * 2*n_aug+1)*/
    VectorXd X_aug = VectorXd(n_aug_);
    X_aug << x_, 0, 0;

    Xsig_pred_ = MatrixXd(n_aug_, 2*n_aug_+ 1);

    /*append the P matrix to find the sigma points*/
    MatrixXd P_aug = MatrixXd(n_aug_ , n_aug_);
    MatrixXd zero1 = MatrixXd(n_x_, number_indp_noise_processes_) * 0;
    MatrixXd zero2 = MatrixXd(number_indp_noise_processes_, n_x_) * 0;

    P_aug <<  P_ , zero1,
            zero2,   Q_;

#ifdef prints
    cout << P_aug << endl;
#endif

    //find the square root matrix of P_aug
    MatrixXd A = P_aug.llt().matrixL();

    Xsig_aug_ = MatrixXd(n_aug_, 2*n_aug_+1) * 0;
    Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1) * 0;
    //Compute the sigma points
    Xsig_aug_.col(0) = X_aug;
    for(unsigned int i = 1; i < n_aug_+1; i++)
    {
      Xsig_aug_.col(i) = X_aug + sqrt(lambda_ + n_aug_) * A.col(i-1);
      Xsig_aug_.col(i+n_aug_) = X_aug - sqrt(lambda_ + n_aug_) * A.col(i-1);
    }
#ifdef prints
    cout << Xsig_aug_ << endl;
#endif

#ifdef print_iterations_
    iterations_++;
    cout << iterations_ << endl;
#endif

    //Predict the mean and covariance using the unscented transform
    Prediction(dt);

    //update
    if(meas_package.sensor_type_ == 0 && use_laser_)
    {
      UpdateLidar(meas_package);
    }
    else if (meas_package.sensor_type_ == 1 && use_radar_)
    {
      UpdateRadar(meas_package);
    }

  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  //predict state mean
  VectorXd x_pred = VectorXd(n_x_) * 0;
  MatrixXd P_pred = MatrixXd(n_x_, n_x_)*0;

  //predict sigma points

  for(unsigned int i = 0; i <  2 * n_aug_ + 1; i++)
  {
    float p_x = Xsig_aug_.col(i)(0);
    float p_y = Xsig_aug_.col(i)(1);
    float v   = Xsig_aug_.col(i)(2);
    float psi = Xsig_aug_.col(i)(3);
    float psi_dot = Xsig_aug_.col(i)(4);
    float a = Xsig_aug_.col(i)(5);
    float psi_ddot = Xsig_aug_.col(i)(6);

    VectorXd temp = VectorXd(5);
    temp(2)  = Xsig_aug_.col(i)(2) +  delta_t * a;
    temp(3)  = Xsig_aug_.col(i)(3) + psi_dot * delta_t + (delta_t * delta_t * psi_ddot / 2);
    temp(4)  = Xsig_aug_.col(i)(4) + psi_ddot * delta_t;

    if( abs(psi_dot) > 0.001)
    {
      temp(0)  = Xsig_aug_.col(i)(0) + (v/psi_dot) * (sin(psi + psi_dot * delta_t) - sin(psi)) + (delta_t * delta_t * cos(psi) * a/2);
      temp(1)  = Xsig_aug_.col(i)(1) + (v/psi_dot) * (-cos(psi + psi_dot * delta_t) + cos(psi)) + (delta_t * delta_t * sin(psi) * a/2);
    }
    else
    {
      temp(0)  = Xsig_aug_.col(i)(0) + (v*cos(psi)*delta_t) + (delta_t * delta_t * cos(psi) * a/2);
      temp(1)  = Xsig_aug_.col(i)(1) + (v*sin(psi)*delta_t) + (delta_t * delta_t * sin(psi) * a/2);
    }
    Xsig_pred_.col(i) = temp;
    x_pred = x_pred + (weights_(i) * Xsig_pred_.col(i));
  }


  for(unsigned int i = 0; i < n_aug_ +1 ; i++)
  {
      VectorXd temp =  (Xsig_pred_.col(i) - x_pred) ;
      while(temp(3) > M_PI) temp(3) = temp(3) - (2*M_PI);
      while(temp(3) < -M_PI) temp(3) = temp(3) + (2*M_PI);
    P_pred = P_pred + weights_(i) * temp * temp.transpose();
    if(i!=0)
    {
        VectorXd temp =  (Xsig_pred_.col(i+n_aug_) - x_pred) ;
        while(temp(3) > M_PI) temp(3) = temp(3) - (2*M_PI);
        while(temp(3) < -M_PI) temp(3) = temp(3) + (2*M_PI);
      P_pred = P_pred + weights_(i+n_aug_) * temp * temp.transpose();
    }

  }

  //write result
  x_ = x_pred;
  P_ = P_pred;



}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  unsigned n_z = 2;
  VectorXd Z = VectorXd(n_z) *0;
  Z << meas_package.raw_measurements_(0),meas_package.raw_measurements_(1);
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1)*0;

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z) * 0;

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  //transform sigma points into measurement space
  for(unsigned int i = 0; i < 2*n_aug_ +1; i++ )
  {
    //Predict measurements
    Zsig.col(i)(0) = Xsig_pred_.col(i)(0);
    Zsig.col(i)(1) = Xsig_pred_.col(i)(1);

    //Accumulate the mean
    z_pred = z_pred + (weights_(i) * Zsig.col(i));

  }
  //calculate mean predicted measurement
  //calculate measurement covariance matrix S
  for (unsigned int i = 0; i < n_aug_ + 1 ; i++)
  {
    //signma point matrix is antisymmetric around n_aug + 1
    if(i == 0)
    {
      S = weights_(i) * ((Zsig.col(i) - z_pred) * (Zsig.col(i) - z_pred).transpose());
    }
    else
    {
      S = S + weights_(i) * ((Zsig.col(i) - z_pred) * (Zsig.col(i) - z_pred).transpose());
      S = S + weights_(i + n_aug_ ) * ((Zsig.col(i + n_aug_) - z_pred) * (Zsig.col(i + n_aug_) - z_pred).transpose());
    }
  }
  MatrixXd R = MatrixXd(n_z,n_z) * 0;
  R << std_laspx_ * std_laspx_, 0,
          0,      std_laspy_ * std_laspy_;


  S = S + R;

  //create matrix for cross correlation Tc
  MatrixXd  Tc = MatrixXd(n_x_, n_z) * 0;

  //calculate cross correlation matrix
  Tc =  weights_(0) * (Xsig_pred_.col(0) - x_) * (Zsig.col(0) - z_pred).transpose();
  for(unsigned int i = 1; i < n_aug_ + 1; i++)
  {
    Tc = Tc + weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();
    Tc = Tc + weights_(i + n_aug_) * (Xsig_pred_.col(i+ n_aug_) - x_) * (Zsig.col(i+ n_aug_) - z_pred).transpose();
  }
  //calculate Kalman gain K;
  MatrixXd k_g = MatrixXd(n_x_, n_z);
  k_g = Tc * S.inverse();
  VectorXd innovation = (Z - z_pred);
  //update state mean and covariance matrix
  x_ =x_ + k_g*innovation;
  P_ = P_ - k_g * S * k_g.transpose();


#ifdef print_measurement_predictions
    cout << "Z      :" << Z(0) << " " << Z(1)  << endl;
    cout << "Z_pred :" << z_pred(0) << " "<< z_pred(1)  << "\n\n" << endl;
#endif

#ifdef NIS_test
  NIS_lidar_ = innovation.transpose() * S.inverse() * innovation;
  NIS_lidar_text_ << NIS_lidar_  <<endl;
#endif
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  unsigned n_z = 3;
  VectorXd Z = VectorXd(n_z) *0;
  Z << meas_package.raw_measurements_(0),meas_package.raw_measurements_(1),meas_package.raw_measurements_(2);
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1)*0;

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z) * 0;

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  //transform sigma points into measurement space
  for(unsigned int i = 0; i < 2*n_aug_ +1; i++ )
  {
    //Predict measurements
    Zsig.col(i)(0) = sqrt(Xsig_pred_.col(i)(0) * Xsig_pred_.col(i)(0) + Xsig_pred_.col(i)(1) * Xsig_pred_.col(i)(1));
    Zsig.col(i)(1) = atan2(Xsig_pred_.col(i)(1) , Xsig_pred_.col(i)(0));
    Zsig.col(i)(2) = (Xsig_pred_.col(i)(0) * Xsig_pred_.col(i)(2) * cos(Xsig_pred_.col(i)(3)) + Xsig_pred_.col(i)(1) * Xsig_pred_.col(i)(2) * sin(Xsig_pred_.col(i)(3)))/ Zsig.col(i)(0);

    //Accumulate the mean
    z_pred = z_pred + (weights_(i) * Zsig.col(i));

  }
  //calculate mean predicted measurement
  //calculate measurement covariance matrix S
  for (unsigned int i = 0; i < n_aug_ + 1 ; i++)
  {
    //signma point matrix is antisymmetric around n_aug + 1
    if(i == 0)
    {
     VectorXd innovation = (Zsig.col(i) - z_pred);
        while(innovation(1) > M_PI) innovation(1) = innovation(1) - (2*M_PI);
        while(innovation(1) < -M_PI) innovation(1) = innovation(1) + (2*M_PI);
      S = weights_(i) * (innovation * innovation.transpose());

    }
    else
    {
        VectorXd innovation = (Zsig.col(i) - z_pred);
        while(innovation(1) > M_PI) innovation(1) = innovation(1) - (2*M_PI);
        while(innovation(1) < -M_PI) innovation(1) = innovation(1) + (2*M_PI);
      S = S + weights_(i) * (innovation * innovation.transpose());

        innovation = (Zsig.col(i + n_aug_) - z_pred);
        while(innovation(1) > M_PI) innovation(1) = innovation(1) - (2*M_PI);
        while(innovation(1) < -M_PI) innovation(1) = innovation(1) + (2*M_PI);
      S = S + weights_(i + n_aug_ ) * (innovation * innovation.transpose());
    }
  }
  MatrixXd R = MatrixXd(n_z,n_z) * 0;
  R << std_radr_ * std_radr_, 0, 0 ,
          0,      std_radphi_ * std_radphi_, 0,
          0 ,       0, std_radrd_ * std_radrd_;

  S = S + R;

  //create matrix for cross correlation Tc
  MatrixXd  Tc = MatrixXd(n_x_, n_z) * 0;

  //calculate cross correlation matrix
    VectorXd innovation = (Zsig.col(0) - z_pred);
    while(innovation(1) > M_PI) innovation(1) = innovation(1) - (2*M_PI);
    while(innovation(1) < -M_PI) innovation(1) = innovation(1) + (2*M_PI);

    VectorXd temp = (Xsig_pred_.col(0) - x_);
    while(temp(3) > M_PI) temp(3) = temp(3) - (2*M_PI);
    while(temp(3) < -M_PI) temp(3) = temp(3) + (2*M_PI);
  Tc =  weights_(0) * temp * innovation.transpose();
  for(unsigned int i = 1; i < n_aug_ + 1; i++)
  {
      VectorXd innovation = (Zsig.col(i) - z_pred);
      while(innovation(1) > M_PI) innovation(1) = innovation(1) - (2*M_PI);
      while(innovation(1) < -M_PI) innovation(1) = innovation(1) + (2*M_PI);

      VectorXd temp = (Xsig_pred_.col(i) - x_);
      while(temp(3) > M_PI) temp(3) = temp(3) - (2*M_PI);
      while(temp(3) < -M_PI) temp(3) = temp(3) + (2*M_PI);
    Tc = Tc + weights_(i) * temp * innovation.transpose();

      innovation = (Zsig.col(i + n_aug_) - z_pred);
      while(innovation(1) > M_PI) innovation(1) = innovation(1) - (2*M_PI);
      while(innovation(1) < -M_PI) innovation(1) = innovation(1) + (2*M_PI);

      temp = (Xsig_pred_.col(i+ n_aug_) - x_);
      while(temp(3) > M_PI) temp(3) = temp(3) - (2*M_PI);
      while(temp(3) < -M_PI) temp(3) = temp(3) + (2*M_PI);
    Tc = Tc + weights_(i + n_aug_) * temp * innovation.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd k_g = MatrixXd(n_x_, n_z);
  k_g = Tc * S.inverse();
  innovation = (Z - z_pred);
    while(innovation(1) > M_PI)
    {
        innovation(1) = innovation(1) - (2*M_PI);
    }
    while(innovation(1) < -M_PI)
    {
        innovation(1) = innovation(1) + (2*M_PI);
    }
  //update state mean and covariance matrix
  x_ =x_ + k_g*innovation;
  P_ = P_ - k_g * S * k_g.transpose();


#ifdef print_measurement_predictions
  cout << "Z      :" << Z(0) <<" " << Z(1) <<" " << Z(2) << endl;
  cout << "Z_pred :" << z_pred(0) <<" " << z_pred(1)<< " " << z_pred(2) << "\n\n" << endl;
#endif

#ifdef NIS_test
  NIS_radar_ = innovation.transpose() * S.inverse() * innovation;
  NIS_radar_text_ << NIS_radar_  <<endl;
#endif



}
