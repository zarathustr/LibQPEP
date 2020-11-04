#include <iostream>
#include "generateProjectedPoints.h"
#include "solver_WQ_1_2_3_4_5_9_13_17_33_49_approx.h"
#include "solver_WQ_approx.h"
#include "solver_WQ.h"
#include "utils.h"
#include "QPEP_grobner.h"
#include "pnp_WQD.h"
#include "pTop_WQD.h"
#include "misc_pnp_funcs.h"
#include "misc_pTop_funcs.h"
#include "QPEP_lm_single.h"
#include "QPEP.h"
#include <fstream>
#include "omp.h"

#ifdef USE_OPENCV
#include <opencv2/core/types.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/highgui/highgui.hpp>
#endif




std::vector<cv::Point3d> Vector3dToPoint3d(std::vector<Eigen::Vector3d> pt)
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

std::vector<cv::Point2d> Vector2dToPoint2d(std::vector<Eigen::Vector2d> pt)
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

void readPnPdata(std::string filename,
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


void test_generateProjectedPoints()
{
    Eigen::Matrix3d R0;
    Eigen::Vector3d t0;
    Eigen::Matrix3d K;
    K.setZero();
    std::vector<Eigen::Vector3d> world_pt0;
    std::vector<Eigen::Vector2d> image_pt0;
    const std::string filename = "/Users/zarathustra/git/LibQPEP/data/pnp_data1.txt";
    readPnPdata(filename, R0, t0, K, world_pt0, image_pt0);

    std::vector<Eigen::Vector2d> image_pt;
    std::vector<double> s;
    generateProjectedPoints(image_pt, s, world_pt0, K, R0, t0);

    for(int i = 0; i < image_pt0.size(); ++i) {
        std::cout.precision(17);
        std::cout << "image_pt0: " << image_pt0[i].transpose() << std::endl;
        std::cout << "image_pt: " << image_pt[i].transpose() << std::endl << std::endl << std::endl;
    }
}

Eigen::Matrix3d R0;
Eigen::Vector3d t0;
Eigen::Matrix3d K;
std::vector<Eigen::Vector3d> world_pt0;
std::vector<Eigen::Vector2d> image_pt0;
const std::string filename = "/Users/zarathustra/git/LibQPEP/data/pnp_data-500pt-1.txt";

void test_pnp_WQD()
{
    if(image_pt0.size() < 3)
    {
        K.setZero();
        readPnPdata(filename, R0, t0, K, world_pt0, image_pt0);
    }

    Eigen::Matrix4d XX;
    XX << R0, t0, Eigen::Vector3d::Zero(3).transpose(), 1.0;
    double problem_scale = 1e-8;
    Eigen::Matrix<double, 4, 64> W;
    Eigen::Matrix<double, 4, 4> Q;
    Eigen::Matrix<double, 3, 37> D;
    Eigen::Matrix<double, 4, 24> coef_f_q_sym;
    Eigen::Matrix<double, 1, 70> coef_J_pure;
    Eigen::Matrix<double, 3, 10> coefs_tq;
    Eigen::Matrix<double, 3, 3> pinvG;
    pnp_WQD(W, Q, D, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG, image_pt0, world_pt0, K, problem_scale);

    Eigen::Matrix<double, 3, 64> W_ = W.topRows(3);
    Eigen::Matrix<double, 3, 4> Q_ = Q.topRows(3);
    W_.row(0) = W.row(0) + W.row(1) + W.row(2);
    W_.row(1) = W.row(1) + W.row(2) + W.row(3);
    W_.row(2) = W.row(2) + W.row(3) + W.row(0);

    Q_.row(0) = Q.row(0) + Q.row(1) + Q.row(2);
    Q_.row(1) = Q.row(1) + Q.row(2) + Q.row(3);
    Q_.row(2) = Q.row(2) + Q.row(3) + Q.row(0);


    Eigen::Matrix3d R;
    Eigen::Vector3d t;
    Eigen::Matrix4d X;
    double min[27];
    struct QPEP_options opt;
    opt.ModuleName = "solver_WQ_1_2_3_4_5_9_13_17_33_49_approx";
    opt.DecompositionMethod = "SparseLU";

    struct QPEP_runtime stat =
            QPEP_WQ_grobner(R, t, X, min, W_, Q_,
                    reinterpret_cast<solver_func_handle>(solver_WQ_1_2_3_4_5_9_13_17_33_49_approx),
                    reinterpret_cast<mon_J_pure_func_handle>(mon_J_pure_pnp_func),
                    reinterpret_cast<t_func_handle>(t_pnp_func),
                    coef_J_pure, coefs_tq, pinvG, nullptr, opt);

    Eigen::Vector4d q0 = R2q(R);
    stat = QPEP_lm_single(R, t, X, q0, 100, 5e-2,
                   reinterpret_cast<eq_Jacob_func_handle>(eq_Jacob_pnp_func),
                   reinterpret_cast<t_func_handle>(t_pnp_func),
                   coef_f_q_sym, coefs_tq, pinvG, stat);

    std::cout << "True X: " << std::endl << XX << std::endl;
    std::cout << "QPEP X: " << std::endl << X << std::endl << std::endl;



//    std::vector<cv::Point2d> image_pt0_cv = Vector2dToPoint2d(image_pt0);
//    std::vector<cv::Point3d> world_pt0_cv = Vector3dToPoint3d(world_pt0);
//    cv::Mat intrinsics, distCoeffs;
//    distCoeffs = cv::Mat::zeros(4, 1, CV_64FC1);
//    for (int i = 0; i < 3; i++)
//        distCoeffs.at<double>(i, 0) = 0.0;
//    intrinsics.create(3, 3, CV_64FC1);
//    intrinsics.at<double>(0, 0) = K(0, 0);
//    intrinsics.at<double>(1, 1) = K(1, 1);
//    intrinsics.at<double>(0, 2) = K(0, 2);
//    intrinsics.at<double>(1, 2) = K(1, 2);
//    intrinsics.at<double>(2, 2) = 1;

//    cv::Mat rvec, tvec;
//    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_EPNP);
//    cv::Mat rotation_matrix = cv::Mat(3, 3, CV_64FC1, cv::Scalar::all(0));
//    cv::Rodrigues(rvec, rotation_matrix);
//    Eigen::Matrix3d Rot;
//    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
//            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
//            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
//    Eigen::Vector3d Trans;
//    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
//    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
//    std::cout << "Opencv EPnP X: " << std::endl << XX << std::endl;
//
//
//    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_P3P);
//    cv::Rodrigues(rvec, rotation_matrix);
//    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
//            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
//            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
//    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
//    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
//    std::cout << "Opencv P3P X: " << std::endl << XX << std::endl;
//
//
//    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_DLS);
//    cv::Rodrigues(rvec, rotation_matrix);
//    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
//            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
//            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
//    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
//    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
//    std::cout << "Opencv DLS X: " << std::endl << XX << std::endl;
//
//
//    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_UPNP);
//    cv::Rodrigues(rvec, rotation_matrix);
//    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
//            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
//            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
//    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
//    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
//    std::cout << "Opencv UPnP X: " << std::endl << XX << std::endl;
//
//
//    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_AP3P);
//    cv::Rodrigues(rvec, rotation_matrix);
//    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
//            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
//            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
//    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
//    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
//    std::cout << "Opencv AP3P X: " << std::endl << XX << std::endl;
//
//
//    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_IPPE);
//    cv::Rodrigues(rvec, rotation_matrix);
//    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
//            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
//            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
//    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
//    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
//    std::cout << "Opencv IPPE X: " << std::endl << XX << std::endl;
//
//
//    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_IPPE_SQUARE);
//    cv::Rodrigues(rvec, rotation_matrix);
//    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
//            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
//            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
//    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
//    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
//    std::cout << "Opencv IPPE SQUARE X: " << std::endl << XX << std::endl;
}



void readpTopdata(std::string filename,
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



std::vector<Eigen::Vector3d> rr0;
std::vector<Eigen::Vector3d> bb0;
std::vector<Eigen::Vector3d> nv0;
const std::string filename2 = "/Users/zarathustra/git/LibQPEP/data/pTop_data-100pt-2.txt";

void test_pTop_WQD() {
    if(rr0.size() < 3) {
        readpTopdata(filename2, R0, t0, rr0, bb0, nv0);
    }
    Eigen::Matrix4d XX;
    XX << R0, t0, Eigen::Vector3d::Zero(3).transpose(), 1.0;

    Eigen::Matrix<double, 4, 64> W;
    Eigen::Matrix<double, 4, 4> Q;
    Eigen::Matrix<double, 3, 28> D;
    Eigen::Matrix<double, 4, 24> coef_f_q_sym;
    Eigen::Matrix<double, 1, 85> coef_J_pure;
    Eigen::Matrix<double, 3, 11> coefs_tq;
    Eigen::Matrix<double, 3, 3> pinvG;
    pTop_WQD(W, Q, D, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG, rr0, bb0, nv0);

    Eigen::Matrix<double, 3, 64> W_ = W.topRows(3);
    Eigen::Matrix<double, 3, 4> Q_ = Q.topRows(3);
    W_.row(0) = W.row(0) + W.row(1) + W.row(2);
    W_.row(1) = W.row(1) + W.row(2) + W.row(3);
    W_.row(2) = W.row(2) + W.row(3) + W.row(0);

    Q_.row(0) = Q.row(0) + Q.row(1) + Q.row(2);
    Q_.row(1) = Q.row(1) + Q.row(2) + Q.row(3);
    Q_.row(2) = Q.row(2) + Q.row(3) + Q.row(0);

    Eigen::Matrix3d R;
    Eigen::Vector3d t;
    Eigen::Matrix4d X;
    double min[27];
    struct QPEP_options opt;
    opt.ModuleName = "solver_WQ_approx";
    opt.DecompositionMethod = "SparseLU";

    struct QPEP_runtime stat =
            QPEP_WQ_grobner(R, t, X, min, W_, Q_,
                    reinterpret_cast<solver_func_handle>(solver_WQ_approx),
                    reinterpret_cast<mon_J_pure_func_handle>(mon_J_pure_pTop_func),
                    reinterpret_cast<t_func_handle>(t_pTop_func),
                    coef_J_pure, coefs_tq, pinvG, nullptr, opt);

//    Eigen::Vector4d q0 = R2q(R);
//
//    stat = QPEP_lm_single(R, t, X, q0, 100, 5e-2,
//                   reinterpret_cast<eq_Jacob_func_handle>(eq_Jacob_pTop_func),
//                   reinterpret_cast<t_func_handle>(t_pTop_func),
//                   coef_f_q_sym, coefs_tq, pinvG, stat);
//
//    std::cout << "True X: " << std::endl << XX << std::endl;
//    std::cout << "QPEP X: " << std::endl << X << std::endl << std::endl;
}


#ifdef USE_OPENCV

void plotQuatCov(cv::Mat& img,
                 const Eigen::Matrix4d& cov1,
                 const Eigen::Matrix4d& cov2,
                 const std::vector<Eigen::Vector4d>& qs,
                 const Eigen::Vector4d& mean_q);

void test_pTop_noise(const std::string& name,
                     cv::Mat& img,
                     const int& num,
                     const double& noise,
                     const bool& verbose) {
    if(rr0.size() < 3) {
        readpTopdata(name, R0, t0, rr0, bb0, nv0);
    }
    Eigen::Matrix4d XX;
    XX << R0, t0, Eigen::Vector3d::Zero(3).transpose(), 1.0;

    std::vector<Eigen::Vector4d> qs(num);
    for(int j = 0; j < num; ++j) {
        std::vector<Eigen::Vector3d> rr(rr0);
        std::vector<Eigen::Vector3d> bb(bb0);
        std::vector<Eigen::Vector3d> nv(nv0);
        for(int i = 0; i < rr.size(); ++i) {
            rr[i] += noise * Eigen::Vector3d::Random();
            bb[i] += noise * Eigen::Vector3d::Random();
            nv[i] += noise * Eigen::Vector3d::Random();
        }
        Eigen::Matrix<double, 4, 64> W;
        Eigen::Matrix<double, 4, 4> Q;
        Eigen::Matrix<double, 3, 28> D;
        Eigen::Matrix<double, 4, 24> coef_f_q_sym;
        Eigen::Matrix<double, 1, 85> coef_J_pure;
        Eigen::Matrix<double, 3, 11> coefs_tq;
        Eigen::Matrix<double, 3, 3> pinvG;
        pTop_WQD(W, Q, D, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG, rr, bb, nv);

        Eigen::Matrix<double, 3, 64> W_ = W.topRows(3);
        Eigen::Matrix<double, 3, 4> Q_ = Q.topRows(3);
        W_.row(0) = W.row(0) + W.row(1) + W.row(2);
        W_.row(1) = W.row(1) + W.row(2) + W.row(3);
        W_.row(2) = W.row(2) + W.row(3) + W.row(0);

        Q_.row(0) = Q.row(0) + Q.row(1) + Q.row(2);
        Q_.row(1) = Q.row(1) + Q.row(2) + Q.row(3);
        Q_.row(2) = Q.row(2) + Q.row(3) + Q.row(0);

        Eigen::Matrix3d R;
        Eigen::Vector3d t;
        Eigen::Matrix4d X;
        double min[27];
        struct QPEP_options opt;
        opt.ModuleName = "solver_WQ_approx";
        opt.DecompositionMethod = "SparseLU";

        struct QPEP_runtime stat =
                QPEP_WQ_grobner(R, t, X, min, W_, Q_,
                                reinterpret_cast<solver_func_handle>(solver_WQ_approx),
                                reinterpret_cast<mon_J_pure_func_handle>(mon_J_pure_pTop_func),
                                reinterpret_cast<t_func_handle>(t_pTop_func),
                                coef_J_pure, coefs_tq, pinvG, nullptr, opt);

        Eigen::Vector4d q0 = R2q(R);

        stat = QPEP_lm_single(R, t, X, q0, 100, 5e-2,
                              reinterpret_cast<eq_Jacob_func_handle>(eq_Jacob_pTop_func),
                              reinterpret_cast<t_func_handle>(t_pTop_func),
                              coef_f_q_sym, coefs_tq, pinvG, stat);

        q0 = R2q(R);
        qs[j] = q0;

        if(verbose)
        {
            std::cout << "True X: " << std::endl << XX << std::endl;
            std::cout << "QPEP X: " << std::endl << X << std::endl << std::endl;
        }
    }

    Eigen::Vector4d mean_q = mean(qs);

    std::vector<Eigen::Vector4d> qs_(qs);
    std::for_each(qs_.begin(), qs_.end(), [=](Eigen::Vector4d& x){x = x - mean_q;});
    Eigen::Matrix4d res(Eigen::Matrix4d::Zero());
    std::for_each(qs_.begin(), qs_.end(), [&res](const Eigen::Vector4d& x){res += x * x.transpose();});
    Eigen::Matrix4d Sigma_data_stat = res * 1.0 / ((double) qs_.size());

    plotQuatCov(img, Sigma_data_stat, Sigma_data_stat, qs, mean_q);
}



cv::RotatedRect getErrorEllipse(double chisquare_val, cv::Point2f mean, cv::Mat covmat){

    //Get the eigenvalues and eigenvectors
    cv::Mat eigenvalues, eigenvectors;
    cv::eigen(covmat, eigenvalues, eigenvectors);

    //Calculate the angle between the largest eigenvector and the x-axis
    double angle = atan2(eigenvectors.at<double>(0,1), eigenvectors.at<double>(0,0));

    //Shift the angle to the [0, 2pi] interval instead of [-pi, pi]
    if(angle < 0)
        angle += 6.28318530718;

    //Conver to degrees instead of radians
    angle = 180*angle/3.14159265359;

    //Calculate the size of the minor and major axes
    double halfmajoraxissize= chisquare_val*sqrt(eigenvalues.at<double>(0));
    double halfminoraxissize= chisquare_val*sqrt(eigenvalues.at<double>(1));

    //Return the oriented ellipse
    //The -angle is used because OpenCV defines the angle clockwise instead of anti-clockwise
    return cv::RotatedRect(mean, cv::Size2f(halfmajoraxissize, halfminoraxissize), angle);

}

void plotDataPoint(cv::Mat& img,
                   const int& x,
                   const int& y,
                   const int& r,
                   const cv::Scalar& color)
{
    cv::circle(img, cv::Point(x, y), r, color, -1);
    return;
}

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
               const Eigen::Vector2d& bias)
{
    cv::Mat cov_cv = (cv::Mat_<double>(2, 2) << cov(a - 1, a - 1), cov(a - 1, b - 1), cov(b - 1, a - 1), cov(b - 1, b - 1));
    cv::Point2d mean_cv(mean(a - 1), mean(b - 1));
    double scale = 5e6;


    std::vector<double> demeaned_x(data.size());
    std::vector<double> demeaned_y(data.size());
    for(int i = 0; i < data.size(); ++i)
    {
        demeaned_x[i] = data[i](a - 1) - mean_cv.x;
        demeaned_y[i] = data[i](b - 1) - mean_cv.y;
    }
    std::sort(demeaned_x.begin(), demeaned_x.end());
    std::sort(demeaned_y.begin(), demeaned_y.end());
    x_min = (demeaned_x[0]);
    x_max = (demeaned_x[data.size() - 1]);
    y_min = (demeaned_y[0]);
    y_max = (demeaned_y[data.size() - 1]);
    scale = 1.0 / MAX4(fabs(x_min), fabs(x_max), fabs(y_min), fabs(y_max));
    scale = scale * size * 0.9;

    double chiqaured = 6.636 * scale;
    cv::RotatedRect ellipse = getErrorEllipse(chiqaured, cv::Point2f(bias.x(), bias.y()), cov_cv);
    cv::ellipse(img, ellipse, ellipse_color, 2, linestyle);

    for(int i = 0; i < data.size(); ++i)
    {
        plotDataPoint(img,
                      (int)((data[i](a - 1) - mean_cv.x) * scale + bias.x()),
                      (int)((data[i](b - 1) - mean_cv.y) * scale + bias.y()),
                      point_size,
                      point_color);
    }
}

void plotQuatCov(cv::Mat& img,
                 const Eigen::Matrix4d& cov1,
                 const Eigen::Matrix4d& cov2,
                 const std::vector<Eigen::Vector4d>& qs,
                 const Eigen::Vector4d& mean_q)
{
    int num_data = qs.size();
    int margin = img.rows / 50;
    double size = img.rows / 6.0 - margin * 2;
    double spacing = 5.0;
    double textlen = 200.0;

    Eigen::Vector2d bias;
    double x_min, x_max, y_min, y_max;
    std::string str_min, str_max;
    char str[512];

    for(int i = 0; i < 4; ++i)
    {
        for(int j = i + 1; j < 4; ++j)
        {
            bias(0) = (2.0 * j - 1.0) * img.rows / 6.0;
            bias(1) = (2.0 * i + 1.0) * img.cols / 6.0;
            plotCov4d(img, x_min, x_max, y_min, y_max,
                      cov1, qs, mean_q, i + 1, j + 1,
                      12, cv::Scalar(0, 0, 0) ,
                      2, cv::Scalar(1, 1, 1) * 120, size, cv::LINE_8, bias);

            double rect_xx, rect_xy;
            double rect_yx, rect_yy;
            rect_xx = (j - 1) * img.cols / 3.0;
            rect_yx = j * img.cols / 3.0;
            rect_xy = i * img.cols / 3.0;
            rect_yy = (i + 1) * img.cols / 3.0;

            cv::rectangle( img,
                           cv::Point( rect_xx + margin, rect_xy + margin ),
                           cv::Point( rect_yx - margin, rect_yy - margin),
                           cv::Scalar( 0, 0, 0 ), 2);

            std::sprintf(str, "q%d", i);
            str_min = std::string(str);
            putText(img, str_min, cv::Point(rect_xx, (rect_xy + rect_yy) / 2.0),
                    cv::FONT_HERSHEY_COMPLEX, 1, cv::Scalar(0, 0, 0), 2);
            std::sprintf(str, "q%d", j);
            str_min = std::string(str);
            putText(img, str_min, cv::Point((rect_yx + rect_xx) / 2.0, rect_yy - margin + 5 * spacing),
                    cv::FONT_HERSHEY_COMPLEX, 1, cv::Scalar(0, 0, 0), 2);


            std::sprintf(str, "%.3e", y_min);
            str_min = std::string(str);
            std::sprintf(str, "%.3e", y_max);
            str_max = std::string(str);
            putText(img, str_min, cv::Point(rect_xx + margin + 1 * spacing, rect_xy + margin + 5 * spacing),
                    cv::FONT_HERSHEY_COMPLEX, 1, cv::Scalar(0, 0, 0), 2);

            putText(img, str_max, cv::Point(rect_xx + margin + 1 * spacing, rect_yy - margin - 2 * spacing),
                    cv::FONT_HERSHEY_COMPLEX, 1, cv::Scalar(0, 0, 0), 2);

            std::sprintf(str, "%.3e", x_min);
            str_min = std::string(str);
            std::sprintf(str, "%.3e", x_max);
            str_max = std::string(str);
            putText(img, str_min, cv::Point(rect_xx + spacing, rect_yy - margin + 5 * spacing),
                    cv::FONT_HERSHEY_COMPLEX, 1, cv::Scalar(0, 0, 0), 2);

            putText(img, str_max, cv::Point(rect_yx - textlen, rect_yy - margin + 5 * spacing),
                    cv::FONT_HERSHEY_COMPLEX, 1, cv::Scalar(0, 0, 0), 2);
        }
    }

    margin = img.rows / 200;
    plotDataPoint(img,
                  0 * img.cols / 3 + 5 * margin, 2 * img.rows / 3 + 3 * margin,
                  8, cv::Scalar(1, 1, 1) * 120);

    putText(img, "Data Points", cv::Point(0 * img.cols / 3 + 10 * margin, 2 * img.rows / 3 + 5 * margin),
            cv::FONT_HERSHEY_COMPLEX, 1, cv::Scalar(0, 0, 0), 2);

    cv::line( img,
          cv::Point(0 * img.cols / 3 + 1 * margin, 2 * img.rows / 3 + 12 * margin),
          cv::Point(0 * img.cols / 3 + 8 * margin, 2 * img.rows / 3 + 12 * margin),
          cv::Scalar( 0, 0, 0 ), 2);

    putText(img, "Stat Quaternion Covariance", cv::Point(0 * img.cols / 3 + 10 * margin, 2 * img.rows / 3 + 14 * margin),
            cv::FONT_HERSHEY_COMPLEX, 1, cv::Scalar(0, 0, 0), 2);
}
#endif

int main(int argc,char ** argv) {

#ifdef USE_OPENCV
    const int row = 1920;
    const int col = 1920;
    cv::Mat imageDraw = cv::Mat::zeros(row, col, CV_8UC3);
    cv::Mat ColorMask(row, col, CV_8UC3, cv::Scalar(1, 1, 1) * 255);
    cv::addWeighted(imageDraw, 0.0, ColorMask, 1.0, 0, imageDraw);

    const bool verbose = false;
    test_pTop_noise("/Users/zarathustra/git/LibQPEP/data/pTop_data-100pt-3.txt",
                    imageDraw, 1500, 1e-4, verbose);

    imshow("imageDraw", imageDraw);
    cv::waitKey(0);
#endif


    std::cout.precision(16);
//    test_generateProjectedPoints();

    double time = 0.0;
    clock_t time1 = clock();
    double loops = 10.0;

    for(int i = 0; i < (int) loops; ++i)
    {
//        test_pnp_WQD();
        test_pTop_WQD();
    }
    clock_t time2 = clock();
    time = time2 - time1;
    std::cout << "Time: " << time / loops / double(CLOCKS_PER_SEC) << std::endl;


    return 0;
}
