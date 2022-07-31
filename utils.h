// LibQPEP: A Library for Globally Optimal Solving Quadratic Pose Estimation Problems,
//          It also gives highly accurate uncertainty description of the solutions.
//
// Authors: Jin Wu and Ming Liu
// Affiliation: Hong Kong University of Science and Technology (HKUST)
// Emails: jin_wu_uestc@hotmail.com; eelium@ust.hk
// Reference: Wu, J., et al. (2022) Quadratic Pose Estimation Problems: 
//                                  Globally Optimal Solutions, 
//                                  Solvability/Observability Analysis,
//                                  and Uncertainty Description.
//                                  IEEE Transactions on Robotics.
//                                  https://doi.org/10.1109/TRO.2022.3155880

#ifndef LIBQPEP_UTILS_H
#define LIBQPEP_UTILS_H

#include <time.h>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>
#include <fstream>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Dense>
#include <Eigen/Geometry> // For Quaternion
#include <Eigen/Sparse>
#include "QPEP.h"

#ifdef USE_SUPERLU
#include <Eigen/SuperLUSupport>
#endif

#ifdef USE_OPENCL
#define VIENNACL_WITH_EIGEN 1
#include <viennacl/matrix.hpp>
#include <viennacl/vector.hpp>
#include <viennacl/compressed_matrix.hpp>
#include <viennacl/linalg/prod.hpp>
#include <viennacl/linalg/direct_solve.hpp>
#include <viennacl/linalg/lu.hpp>

template<typename T>
struct Eigen_dense_matrix
{
  typedef typename T::ERROR_NO_EIGEN_TYPE_AVAILABLE   error_type;
};

template<>
struct Eigen_dense_matrix<float>
{
  typedef Eigen::MatrixXf  type;
};
template<>
struct Eigen_dense_matrix<double>
{
  typedef Eigen::MatrixXd  type;
};

#endif


inline Eigen::VectorXd vec(const Eigen::MatrixXd X)
{
    Eigen::VectorXd res;
    int m = X.cols();
    int n = X.rows();
    res.resize(m * n);
    int counter = 0;
    for(int i = 0; i < m; ++i)
    {
        for(int j = 0; j < n; ++j) {
            res(counter) = X(j, i);
            counter = counter + 1;
        }
    }
    return res;
}

inline Eigen::Matrix3d q2R(const Eigen::Quaterniond& q)
{
    double q0 = q.x();
    double q1 = q.y();
    double q2 = q.z();
    double q3 = q.w();
    double q02 = q0 * q0;
    double q12 = q1 * q1;
    double q22 = q2 * q2;
    double q32 = q3 * q3;

    Eigen::Matrix3d R;
    R <<     q02 + q12 - q22 - q32,     2.0*q0*q3 + 2.0*q1*q2,     2.0*q1*q3 - 2.0*q0*q2,
             2.0*q1*q2 - 2.0*q0*q3,     q02 - q12 + q22 - q32,     2.0*q0*q1 + 2.0*q2*q3,
             2.0*q0*q2 + 2.0*q1*q3,     2.0*q2*q3 - 2.0*q0*q1,     q02 - q12 - q22 + q32;
    return R;
}

inline Eigen::Vector4d R2q(const Eigen::Matrix3d& R)
{
    //double det = R.determinant();
    //double orthogonality = (R.transpose() * R).trace();
    //assert(std::fabs(det - 1.0) < 1e-12 && std::fabs(orthogonality - 3.0) < 1e-12);

    double G11 = R(0, 0) + R(1, 1) + R(2, 2) - 3.0, G12 = R(1, 2) - R(2, 1), G13 = R(2, 0) - R(0, 2), G14 = R(0, 1) - R(1, 0);
    double G22 = R(0, 0) - R(1, 1) - R(2, 2) - 3.0, G23 = R(0, 1) + R(1, 0), G24 = R(0, 2) + R(2, 0);
    double G33 = R(1, 1) - R(0, 0) - R(2, 2) - 3.0, G34 = R(1, 2) + R(2, 1);

    Eigen::Vector4d qRes = Eigen::Vector4d (
            G14 * G23 * G23 - G13 * G23 * G24 - G14 * G22 * G33 + G12 * G24 * G33 + G13 * G22 * G34 - G12 * G23 * G34,
            G13 * G13 * G24 + G12 * G14 * G33 - G11 * G24 * G33 + G11 * G23 * G34 - G13 * G14 * G23 - G13 * G12 * G34,
            G13 * G14 * G22 - G12 * G14 * G23 - G12 * G13 * G24 + G11 * G23 * G24 + G12 * G12 * G34 - G11 * G22 * G34,
            - ( G13 * G13 * G22 - 2 * G12 * G13 * G23 + G11 * G23 * G23 + G12 * G12 * G33 - G11 * G22 * G33 ));
    qRes.normalize();
    if(qRes(0) < 0)
        qRes = - qRes;
    return qRes;
}

#define MAX2(A, B) (A > B ? A : B)
#define MAX4(A, B, C, D) (MAX2(MAX2(A, B), MAX2(C, D)))


typedef struct node
{
    double value;
    int index;
} node;

inline bool cmp(struct node a, struct node b)
{
    if (a.value < b.value)
    {
        return true;
    }
    return false;
}

template <typename T>
T sort_indices(std::vector<size_t> &idx, const std::vector<T> &v)
{
    node* a = new node[v.size()];
    for (int i = 0; i < v.size(); i++)
    {
        a[i].value = v[i];
        a[i].index = i;
    }

    std::sort(a, a + v.size(), cmp);
    for (int i = 0; i < v.size(); i++)
    {
        idx.push_back(a[i].index);
    }
    delete[] a;

    return 0;
}

inline double powers(double x, double order)
{
    if(order == 2.0){
        return x * x;
    }
    else if(order == 3.0){
        return x * x * x;
    }
    else if(order == 4.0) {
        double x2 = x * x;
        return x2 * x2;
    }
    else{
        return std::pow(x, order);
    }
}

typedef void (* data_func_handle)(Eigen::SparseMatrix<double>& C1, Eigen::MatrixXd& C2, const Eigen::VectorXd& data);

void readpTopdata(const std::string& filename,
                  Eigen::Matrix3d& R0,
                  Eigen::Vector3d& t0,
                  std::vector<Eigen::Vector3d>& r0,
                  std::vector<Eigen::Vector3d>& b0,
                  std::vector<Eigen::Vector3d>& nv);

void readHandEyedata(const std::string& filename,
                  Eigen::Matrix3d& R0,
                  Eigen::Vector3d& t0,
                  std::vector<Eigen::Matrix4d>& As,
                  std::vector<Eigen::Matrix4d>& Bs);

void readPnPdata(const std::string& filename,
                 Eigen::Matrix3d& R0,
                 Eigen::Vector3d& t0,
                 Eigen::Matrix3d& K,
                 std::vector<Eigen::Vector3d>& world_pt0,
                 std::vector<Eigen::Vector2d>& image_pt0);

QPEP_runtime GaussJordanElimination(
        Eigen::MatrixXd& C1_,
        const Eigen::VectorXd& data,
        const data_func_handle data_func,
        const int& size_GJ,
        const int& size_AM,
        const struct QPEP_options& opt,
        const struct QPEP_runtime& stat_);

template <typename T>
T mean(std::vector<T> data)
{
    int num = data.size();
    assert(num > 0);

    Eigen::MatrixXd s = data[0];
    double factor = 1.0 / ((double) num);
    if(num > 1)
        std::for_each(data.begin() + 1, data.end(), [&s](Eigen::MatrixXd x){s += x;});
    return s * factor;
}

Eigen::Matrix3d orthonormalize(const Eigen::Matrix3d& R);
Eigen::MatrixXd randomMatrix(const int& dim1, const int& dim2, const int& resolution);

#ifdef USE_OPENCV
#include <opencv2/core/types.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/highgui/highgui.hpp>
std::vector<cv::Point3d> Vector3dToPoint3d(std::vector<Eigen::Vector3d> pt);
std::vector<cv::Point2d> Vector2dToPoint2d(std::vector<Eigen::Vector2d> pt);
#endif

#endif //LIBQPEP_UTILS_H
