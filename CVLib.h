#ifndef LIBQPEP_CVLIB_H
#define LIBQPEP_CVLIB_H

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Geometry> // For Quaternion
#include <Eigen/Sparse>

#ifdef USE_OPENCV
#include <opencv2/core/types.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/highgui/highgui.hpp>
#if WIN32
#include <windows.h>
#else
#if !defined(OSX_10_9) && !defined(OSX_BIG_SUR)
#include <X11/Xlib.h>
#endif
#endif

void getScreenResolution(int &width, int &height);

void plotQuatCov(cv::Mat& img,
                 const Eigen::Matrix4d& cov1,
                 const Eigen::Matrix4d& cov2,
                 const std::vector<Eigen::Vector4d>& qs,
                 const Eigen::Vector4d& mean_q,
                 const double& fontsize);

void plotTransCov(cv::Mat& img,
                 const Eigen::Matrix3d& cov1,
                 const Eigen::Matrix3d& cov2,
                 const std::vector<Eigen::Vector3d>& qs,
                 const Eigen::Vector3d& mean_q,
                 const double& fontsize);

#endif

#endif //LIBQPEP_CVLIB_H
