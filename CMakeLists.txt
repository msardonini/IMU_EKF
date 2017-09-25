
cmake_minimum_required(VERSION 2.8)

project(ROBOTICSCAPE)

set(IMU_EKF_DIR lib/IMU_EKF/)

include_directories(
	${IMU_EKF_DIR}
)

set(IMU_EKF_SRCS
	${IMU_EKF_DIR}IMU_EKF.cpp
	${IMU_EKF_DIR}IMU_EKF2.cpp
)


add_library(IMU_EKF SHARED ${IMU_EKF_SRCS})
