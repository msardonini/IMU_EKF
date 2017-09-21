%% Clear all
clear all
close all
clc
%% Variable Declaration
P=single(eye(4));                       %Covariance Matrix
q=single(ones(4,1));                    %Quaternion Vector
Cov_info=single([0.00001;0.1;1]);       %Noise Variances information (q,r_acc,r_mag) < - - Assumes same q for all quat.
accel=single(rand(3,1));                %Accelerometer Vector (x,y,z)
omega=single(rand(3,1));                %Gyro Vector (x,y,z)
mag=single(rand(3,1));                  %Magnetometer Vector (x,y,z)
dt=single(0.02);                        %Sampling Time (must be measured outside)
ini=int8(1);                            %Initialize/Reset Covariance Matrix and Quaternion Vector
use_mag=int8(1);                        %Use magnetometer or not
%% Run Function
[P,q]=IMU_EKF(P,q,Cov_info,omega,accel,mag,dt,ini,use_mag)
Cov_info=single([1;0.1;1]);       %Noise Variances information (q,r_acc,r_mag) < - - Assumes same q for all quat.
[P,q]=IMU_EKF2(P,q,Cov_info,omega,accel,mag,dt,ini,use_mag)