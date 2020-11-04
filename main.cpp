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
#include "omp.h"

#ifdef USE_OPENCV
#include <opencv2/core/types.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "CVLib.h"
#endif


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

static Eigen::Matrix3d R0;
static Eigen::Vector3d t0;
static Eigen::Matrix3d K;
static std::vector<Eigen::Vector3d> world_pt0;
static std::vector<Eigen::Vector2d> image_pt0;
std::vector<Eigen::Vector3d> rr0;
std::vector<Eigen::Vector3d> bb0;
std::vector<Eigen::Vector3d> nv0;

void test_pnp_WQD(const std::string& filename,
                  const bool& verbose)
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

    if(verbose) {
        std::cout << "True X: " << std::endl << XX << std::endl;
        std::cout << "QPEP X: " << std::endl << X << std::endl << std::endl;
    }


#ifdef USE_OPENCV
    std::vector<cv::Point2d> image_pt0_cv = Vector2dToPoint2d(image_pt0);
    std::vector<cv::Point3d> world_pt0_cv = Vector3dToPoint3d(world_pt0);
    cv::Mat intrinsics, distCoeffs;
    distCoeffs = cv::Mat::zeros(4, 1, CV_64FC1);
    for (int i = 0; i < 3; i++)
        distCoeffs.at<double>(i, 0) = 0.0;
    intrinsics.create(3, 3, CV_64FC1);
    intrinsics.at<double>(0, 0) = K(0, 0);
    intrinsics.at<double>(1, 1) = K(1, 1);
    intrinsics.at<double>(0, 2) = K(0, 2);
    intrinsics.at<double>(1, 2) = K(1, 2);
    intrinsics.at<double>(2, 2) = 1;

    cv::Mat rvec, tvec;
    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_EPNP);
    cv::Mat rotation_matrix = cv::Mat(3, 3, CV_64FC1, cv::Scalar::all(0));
    cv::Rodrigues(rvec, rotation_matrix);
    Eigen::Matrix3d Rot;
    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
    Eigen::Vector3d Trans;
    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
    if(verbose) {
        std::cout << "Opencv EPnP X: " << std::endl << XX << std::endl;
    }


    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_P3P);
    cv::Rodrigues(rvec, rotation_matrix);
    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
    if(verbose) {
        std::cout << "Opencv P3P X: " << std::endl << XX << std::endl;
    }


    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_DLS);
    cv::Rodrigues(rvec, rotation_matrix);
    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
    if(verbose) {
        std::cout << "Opencv DLS X: " << std::endl << XX << std::endl;
    }


    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_UPNP);
    cv::Rodrigues(rvec, rotation_matrix);
    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
    if(verbose) {
        std::cout << "Opencv UPnP X: " << std::endl << XX << std::endl;
    }


    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_AP3P);
    cv::Rodrigues(rvec, rotation_matrix);
    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
    if(verbose) {
        std::cout << "Opencv AP3P X: " << std::endl << XX << std::endl;
    }


    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_IPPE);
    cv::Rodrigues(rvec, rotation_matrix);
    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
    if(verbose) {
        std::cout << "Opencv IPPE X: " << std::endl << XX << std::endl;
    }


    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_IPPE_SQUARE);
    cv::Rodrigues(rvec, rotation_matrix);
    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
    if(verbose) {
        std::cout << "Opencv IPPE Square X: " << std::endl << XX << std::endl;
    }
#endif
}



void test_pTop_WQD(const std::string& filename,
                   const bool& verbose) {
    if(rr0.size() < 3) {
        readpTopdata(filename, R0, t0, rr0, bb0, nv0);
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

    Eigen::Vector4d q0 = R2q(R);

    stat = QPEP_lm_single(R, t, X, q0, 100, 5e-2,
                   reinterpret_cast<eq_Jacob_func_handle>(eq_Jacob_pTop_func),
                   reinterpret_cast<t_func_handle>(t_pTop_func),
                   coef_f_q_sym, coefs_tq, pinvG, stat);

    if(verbose)
    {
        std::cout << "True X: " << std::endl << XX << std::endl;
        std::cout << "QPEP X: " << std::endl << X << std::endl << std::endl;
    }
}


#ifdef USE_OPENCV

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
//        test_pnp_WQD("/Users/zarathustra/git/LibQPEP/data/pTop_data-5pt-1.txt", false);
        test_pTop_WQD("/Users/zarathustra/git/LibQPEP/data/pTop_data-100pt-3.txt", false);
    }
    clock_t time2 = clock();
    time = time2 - time1;
    std::cout << "Time: " << time / loops / double(CLOCKS_PER_SEC) << std::endl;


    return 0;
}
