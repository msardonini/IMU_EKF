# IMU_EKF
This code repliactes the work presented in the paper "Double Stage Kalman Filter for IMU.pdf" 
or "A Double-Stage Kalman Filter for Orientation Tracking With an Integrated Processor in 9-D IMU" 
via a Matlab code. 

The matlab code was then converted to C code using Matlab C Coder 
and then renamed to .cpp just to have it in cpp format.

It receives 9 parameters which must be initialized as follows (description given)

float P[16];                //Covariance Matrix 4x4

float quat[4];              //Quaternion Vector 4x1

float Cov_info[3];          //Noise Covariance Info 3x1 (q,r_acc,r_mag) 

                            //(q=noise variance of quaternion (equal for all)
                            
                            //(r_acc=noise variance of accelerometer measurement)
                            
                            //(r_mag=noise variance of magnetometer measurement)
                            
float omega[3];             //Gyroscpe Vector 3x1 (In radians!)

float accel[3];             //Accelerometer Vector 3x1 (It normalizes it, so doesn't matter)

float mag[3];               //Magnetometer Vector 3x1 (It normalizes it, so doesn't matter)

float dt=0.02;              //Sampling time (To be measured outside. Can be just constant though...)

int ini=1;                  //Reset/Initialize Covariance/Quaternion. 

                            //You must select 1 this at startup!!! Then 0..
                            
int use_mag=1               //Decide whether to use mag or not. (1=Yes)

Once that is done you can call it as, 

IMU_EKF(Pcov,quat,Q_info,R_info,gyro_rad,accel,mag,dt,ini,use_mag);
