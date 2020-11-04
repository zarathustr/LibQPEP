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



void plotCov4d(cv::Mat& img,
               double& x_min,
               double& x_max,
               double& y_min,
               double& y_max,
               const Eigen::Matrix4d& cov,
               const std::vector<Eigen::Vector4d>& data,
               const Eigen::Vector4d& mean,
               const int& a,
               const int& b,
               const int& ellipse_size,
               const cv::Scalar& ellipse_color,
               const int& point_size,
               const cv::Scalar& point_color,
               const double& size,
               const int& linestyle,
               const Eigen::Vector2d& bias);

void plotQuatCov(cv::Mat& img,
                 const Eigen::Matrix4d& cov1,
                 const Eigen::Matrix4d& cov2,
                 const std::vector<Eigen::Vector4d>& qs,
                 const Eigen::Vector4d& mean_q);

#endif

#endif //LIBQPEP_CVLIB_H
