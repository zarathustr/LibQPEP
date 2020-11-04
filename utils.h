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
    double det = R.determinant();
    double orthogonality = (R.transpose() * R).trace();
    assert(abs(det - 1.0) < 1e-12 && abs(orthogonality - 3.0) < 1e-12);

    double G11 = R(0, 0) + R(1, 1) + R(2, 2) - 3.0, G12 = R(1, 2) - R(2, 1), G13 = R(2, 0) - R(0, 2), G14 = R(0, 1) - R(1, 0);
    double G22 = R(0, 0) - R(1, 1) - R(2, 2) - 3.0, G23 = R(0, 1) + R(1, 0), G24 = R(0, 2) + R(2, 0);
    double G33 = R(1, 1) - R(0, 0) - R(2, 2) - 3.0, G34 = R(1, 2) + R(2, 1);

    Eigen::Vector4d qRes = Eigen::Vector4d (
            G14 * G23 * G23 - G13 * G23 * G24 - G14 * G22 * G33 + G12 * G24 * G33 + G13 * G22 * G34 - G12 * G23 * G34,
            G13 * G13 * G24 + G12 * G14 * G33 - G11 * G24 * G33 + G11 * G23 * G34 - G13 * G14 * G23 - G13 * G12 * G34,
            G13 * G14 * G22 - G12 * G14 * G23 - G12 * G13 * G24 + G11 * G23 * G24 + G12 * G12 * G34 - G11 * G22 * G34,
            - ( G13 * G13 * G22 - 2 * G12 * G13 * G23 + G11 * G23 * G23 + G12 * G12 * G33 - G11 * G22 * G33 ));
    qRes.normalize();
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

inline QPEP_runtime GaussJordanElimination(
        Eigen::MatrixXd& C1_,
        const Eigen::VectorXd& data,
        const data_func_handle data_func,
        const int& size_GJ,
        const int& size_AM,
        const struct QPEP_options& opt,
        const struct QPEP_runtime& stat_)
{
    assert(opt.DecompositionMethod == "SparseLU" ||
           opt.DecompositionMethod == "SparseQR" ||
           opt.DecompositionMethod == "HouseholderQR" ||
           opt.DecompositionMethod == "PartialPivLU" ||
           opt.DecompositionMethod == "SVD" ||
           opt.DecompositionMethod == "BDCSVD" ||
           opt.DecompositionMethod == "Inv" ||
           opt.DecompositionMethod == "Cholesky");

    struct QPEP_runtime stat = stat_;
    C1_.resize(size_GJ, size_AM);
    C1_.setZero();
    Eigen::MatrixXd C2;
    C2.resize(size_GJ, size_AM);
    C2.setZero();

    if (opt.DecompositionMethod == "SparseLU")
    {
        clock_t time1 = clock();
        Eigen::SparseMatrix<double> C1(size_GJ, size_GJ);
        C1.setZero();
        data_func(C1, C2, data);
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
        solver.compute(C1);
        if (solver.info() != Eigen::Success) {
            std::cout << "Sparse LU Decomposition Failed!" << std::endl;
            stat.statusDecomposition = -3;
            return stat;
        }

        for (int i = 0; i < C2.cols(); ++i) {
            Eigen::VectorXd c1 = solver.solve(C2.col(i));
            if (solver.info() != Eigen::Success) {
                std::cout << "Least Squares Failed!" << std::endl;
                stat.statusDecomposition = -6;
                return stat;
            }
            C1_.col(i) = c1;
        }
        clock_t time2 = clock();
        stat.timeDecomposition = (time2 - time1) / double(CLOCKS_PER_SEC);
    }
    else if (opt.DecompositionMethod == "SparseQR")
    {
        clock_t time1 = clock();
        Eigen::SparseMatrix<double> C1(size_GJ, size_GJ);
        C1.setZero();
        data_func(C1, C2, data);
        C1.makeCompressed();
        Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
        solver.compute(C1);
        if (solver.info() != Eigen::Success) {
            std::cout << "Sparse QR Decomposition Failed!" << std::endl;
            stat.statusDecomposition = -3;
            return stat;
        }

        for (int i = 0; i < C2.cols(); ++i) {
            Eigen::VectorXd c1 = solver.solve(C2.col(i));
            if (solver.info() != Eigen::Success) {
                std::cout << "Least Squares Failed!" << std::endl;
                stat.statusDecomposition = -6;
                return stat;
            }
            C1_.col(i) = c1;
        }
        clock_t time2 = clock();
        stat.timeDecomposition = (time2 - time1) / double(CLOCKS_PER_SEC);
    }
    else if(opt.DecompositionMethod == "HouseholderQR")
    {
        clock_t time1 = clock();
        Eigen::SparseMatrix<double> C1(size_GJ, size_GJ);
        C1.setZero();
        data_func(C1, C2, data);
        Eigen::MatrixXd CC1(C1.toDense());
        Eigen::HouseholderQR<Eigen::MatrixXd > solver;
        solver.compute(CC1);
        for(int i = 0; i < C2.cols(); ++i)
        {
            Eigen::VectorXd c1 = solver.solve(C2.col(i));
            C1_.col(i) = c1;
        }
        clock_t time2 = clock();
        stat.timeDecomposition = (time2 - time1) / double(CLOCKS_PER_SEC);
    }
    else if(opt.DecompositionMethod == "PartialPivLU")
    {
        clock_t time1 = clock();
        Eigen::SparseMatrix<double> C1(size_GJ, size_GJ);
        C1.setZero();
        data_func(C1, C2, data);
        Eigen::MatrixXd CC1(C1.toDense());
        Eigen::PartialPivLU<Eigen::MatrixXd > solver;
        solver.compute(CC1);
        for(int i = 0; i < C2.cols(); ++i)
        {
            Eigen::VectorXd c1 = solver.solve(C2.col(i));
            C1_.col(i) = c1;
        }
        clock_t time2 = clock();
        stat.timeDecomposition = (time2 - time1) / double(CLOCKS_PER_SEC);
    }
    else if(opt.DecompositionMethod == "SVD") {
        clock_t time1 = clock();
        Eigen::SparseMatrix<double> C1(size_GJ, size_GJ);
        C1.setZero();
        data_func(C1, C2, data);
        Eigen::MatrixXd CC1(C1.toDense());
        Eigen::JacobiSVD<Eigen::MatrixXd> solver(CC1, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::MatrixXd singularValues = solver.singularValues();
        Eigen::MatrixXd singularValuesInv;
        singularValuesInv.resize(size_GJ, size_GJ);
        singularValuesInv.setZero();
        double pinvtoler = 1.e-9; // choose your tolerance wisely
        for (int i = 0; i < size_GJ; ++i) {
            if (singularValues(i) > pinvtoler)
                singularValuesInv(i, i) = 1.0 / singularValues(i);
            else
                singularValuesInv(i, i) = 0.0;
        }
        Eigen::MatrixXd pinv = solver.matrixV() * singularValuesInv * solver.matrixU().transpose();

        for (int i = 0; i < C2.cols(); ++i) {
            Eigen::VectorXd c1 = pinv * C2.col(i);
            C1_.col(i) = c1;
        }
        clock_t time2 = clock();
        stat.timeDecomposition = (time2 - time1) / double(CLOCKS_PER_SEC);
    }
    else if(opt.DecompositionMethod == "BDCSVD") {
        clock_t time1 = clock();
        Eigen::SparseMatrix<double> C1(size_GJ, size_GJ);
        C1.setZero();
        data_func(C1, C2, data);
        Eigen::MatrixXd CC1(C1.toDense());
        for (int i = 0; i < C2.cols(); ++i) {
            Eigen::VectorXd c1 = CC1.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(C2.col(i));
            C1_.col(i) = c1;
        }
        clock_t time2 = clock();
        stat.timeDecomposition = (time2 - time1) / double(CLOCKS_PER_SEC);
    }
    else if(opt.DecompositionMethod == "Inv") {
        clock_t time1 = clock();
        Eigen::SparseMatrix<double> C1(size_GJ, size_GJ);
        C1.setZero();
        data_func(C1, C2, data);
        Eigen::MatrixXd CC1(C1.toDense());
        Eigen::MatrixXd inv = CC1.inverse();
        for (int i = 0; i < C2.cols(); ++i) {
            Eigen::VectorXd c1 = inv * C2.col(i);
            C1_.col(i) = c1;
        }
        clock_t time2 = clock();
        stat.timeDecomposition = (time2 - time1) / double(CLOCKS_PER_SEC);
    }
    else if(opt.DecompositionMethod == "Cholesky") {
        clock_t time1 = clock();
        Eigen::SparseMatrix<double> C1(size_GJ, size_GJ);
        C1.setZero();
        data_func(C1, C2, data);
        Eigen::MatrixXd CC1(C1.toDense());
        Eigen::LLT<Eigen::MatrixXd> solver(CC1);
        for (int i = 0; i < C2.cols(); ++i) {
            Eigen::VectorXd c1 = solver.solve(C2.col(i));
            C1_.col(i) = c1;
        }
        clock_t time2 = clock();
        stat.timeDecomposition = (time2 - time1) / double(CLOCKS_PER_SEC);
    }
    return stat;
}

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

inline void readPnPdata(std::string filename,
                 Eigen::Matrix3d& R0,
                 Eigen::Vector3d& t0,
                 Eigen::Matrix3d& K,
                 std::vector<Eigen::Vector3d>& world_pt0,
                 std::vector<Eigen::Vector2d>& image_pt0)
{
    std::ifstream input(filename);
    double fx, fy, cx, cy;
    input >> R0(0, 0) >> R0(0, 1) >> R0(0, 2) >>
          R0(1, 0) >> R0(1, 1) >> R0(1, 2) >>
          R0(2, 0) >> R0(2, 1) >> R0(2, 2);
    input >> t0(0) >> t0(1) >> t0(2);
    input >> fx >> fy >> cx >> cy;
    K(0, 0) = fx;
    K(1, 1) = fy;
    K(0, 2) = cx;
    K(1, 2) = cy;
    K(2, 2) = 1.0;
    int num = 0;
    input >> num;
    world_pt0.resize(num);
    image_pt0.resize(num);
    for(int i = 0; i < num; ++i)
    {
        input >> world_pt0[i](0) >> world_pt0[i](1) >> world_pt0[i](2);
    }
    for(int i = 0; i < num; ++i)
    {
        input >> image_pt0[i](0) >> image_pt0[i](1);
    }
    input.close();
}


inline void readpTopdata(std::string filename,
                  Eigen::Matrix3d& R0,
                  Eigen::Vector3d& t0,
                  std::vector<Eigen::Vector3d>& r0,
                  std::vector<Eigen::Vector3d>& b0,
                  std::vector<Eigen::Vector3d>& nv)
{
    std::ifstream input(filename);
    input >> R0(0, 0) >> R0(0, 1) >> R0(0, 2) >>
          R0(1, 0) >> R0(1, 1) >> R0(1, 2) >>
          R0(2, 0) >> R0(2, 1) >> R0(2, 2);
    input >> t0(0) >> t0(1) >> t0(2);
    int num = 0;
    input >> num;
    r0.resize(num);
    b0.resize(num);
    nv.resize(num);
    for(int i = 0; i < num; ++i)
    {
        input >> r0[i](0) >> r0[i](1) >> r0[i](2);
    }
    for(int i = 0; i < num; ++i)
    {
        input >> b0[i](0) >> b0[i](1) >> b0[i](2);
    }
    for(int i = 0; i < num; ++i)
    {
        input >> nv[i](0) >> nv[i](1) >> nv[i](2);
    }
    input.close();
}

#ifdef USE_OPENCV
#include <opencv2/core/types.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/highgui/highgui.hpp>

inline std::vector<cv::Point3d> Vector3dToPoint3d(std::vector<Eigen::Vector3d> pt)
{
    std::vector<cv::Point3d> tmp;
    for(int i = 0; i < pt.size(); ++i)
    {
        Eigen::Vector3d point = pt[i];
        cv::Point3d vec;
        vec.x = point(0);
        vec.y = point(1);
        vec.z = point(2);
        tmp.push_back(vec);
    }
    return tmp;
}

inline std::vector<cv::Point2d> Vector2dToPoint2d(std::vector<Eigen::Vector2d> pt)
{
    std::vector<cv::Point2d> tmp;
    for(int i = 0; i < pt.size(); ++i)
    {
        Eigen::Vector2d point = pt[i];
        cv::Point2d vec;
        vec.x = point(0);
        vec.y = point(1);
        tmp.push_back(vec);
    }
    return tmp;
}
#endif

#endif //LIBQPEP_UTILS_H
