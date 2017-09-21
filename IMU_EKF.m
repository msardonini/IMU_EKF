function [P,q]=IMU_EKF(P,q,Cov_info,omega,accel,mag,dt,ini,use_mag)
%P: 4x4 Covariance Matrix

%q: 4x1 Quaternion Vector

%Cov_info: 3x1 Vector
%Cov_info(1)=Covariance of unmodeled errors in the quaternion
%Cov_info(2)=Covariance of measurements in the accelerometer
%Cov_info(3)=Covariance of measurements in the magnetometer

%omega: 3x1 Gyro Rate Vector

%accel: 3x1 Acceleration Vector

%mag: 3x1 Magnetometer Vector

%dt: sampling time (not measured inside)

%ini: System Reset/Initialization

%use: Inidicator of which corrections to do
%use_mag=Use Magnetometer

%Initialize/Reset System
if ini==1
    P=single(10*eye(4));
    q=single(zeros(4,1));
    q(1)=1;
end

%Propagate Quaternion
o=omega;
Ak=single(dt/2*[0,-o(1),-o(2),-o(3);o(1),0,o(3),-o(2);o(2),-o(3),0,o(1);o(3),o(2),-o(1),0]+eye(4));
q=Ak*q;

%Propagate Covariance
P=Ak*P*Ak';
for i=1:4
    P(i,i)=P(i,i)+Cov_info(1);
end

%Correction Stage 1 - Use Accelerometer for roll/pitch correction
accel=accel/norm(accel);
Hk=[-2*q(3),2*q(4),-2*q(1),2*q(2);
    2*q(2),2*q(1),2*q(4),2*q(3);
    2*q(1),-2*q(2),-2*q(3),2*q(4)];
S=Hk*P*Hk';
for i=1:3
    S(i,i)=S(i,i)+Cov_info(2);
end
Kk=P*Hk'/S;
zk=accel;
h=[2*q(2)*q(4)-2*q(1)*q(3);
    2*q(1)*q(2)+2*q(3)*q(4);
    q(1)^2-q(2)^2-q(3)^2+q(4)^2];
qe=Kk*(zk-h);
qe(4)=0;                %Don't correct heading with accelerometer
q=q+qe;
P=P-Kk*Hk*P;

%Correction stage 2 - Use Magnetometer for heading correction
if use_mag==1
mag=mag/norm(mag);
Hk=[2*q(4),2*q(3),2*q(2),2*q(1);
    2*q(1),-2*q(2),-2*q(3),-2*q(4);
    -2*q(2),-2*q(1),2*q(4),2*q(3)];
S=Hk*P*Hk';
for i=1:3
    S(i,i)=S(i,i)+Cov_info(3);
end
Kk=P*Hk'/S;
zk=mag;
h=[2*q(2)*q(3)+2*q(1)*q(4);
    q(1)^2-q(2)^2-q(3)^2-q(4)^2;
    2*q(3)*q(4)-2*q(1)*q(2)];
qe=Kk*(zk-h);
qe(2)=0;            %Don't correct roll/pitch
qe(3)=0;            %Don't correct roll/pitch
q=q+qe;
P=P-Kk*Hk*P;
end

%Normalize Quaternion
q=q/norm(q);
end