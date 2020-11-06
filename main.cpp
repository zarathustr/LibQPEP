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

#include "CovEstimation.h"
#include <Eigen/../unsupported/Eigen/KroneckerProduct>



void test_generateProjectedPoints()
{
    Eigen::Matrix3d R0;
    Eigen::Vector3d t0;
    Eigen::Matrix3d K;
    K.setZero();
    std::vector<Eigen::Vector3d> world_pt0;
    std::vector<Eigen::Vector2d> image_pt0;
    const std::string filename = "../data/pnp_data1.txt";
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

#ifdef USE_OPENCV4
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
    Eigen::Matrix<double, 3, 9> G;
    Eigen::Vector3d c;
    Eigen::Matrix<double, 4, 24> coef_f_q_sym;
    Eigen::Matrix<double, 1, 85> coef_J_pure;
    Eigen::Matrix<double, 3, 11> coefs_tq;
    Eigen::Matrix<double, 3, 3> pinvG;
    pTop_WQD(W, Q, D, G, c, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG, rr0, bb0, nv0);

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


Eigen::MatrixXd covx(const std::vector<Eigen::MatrixXd>& x,
                     const std::vector<Eigen::MatrixXd>& y)
{
    int len = x.size();
    std::vector<Eigen::MatrixXd> res(len);
    Eigen::MatrixXd mean_x = mean(x);
    Eigen::MatrixXd mean_y = mean(y);
    for(int i = 0; i < len; ++i)
    {
        res[i].resize(x[0].rows(), y[0].rows());
        res[i] = (x[i] - mean_x) * (y[i] - mean_y).transpose();
    }
    return mean(res);
}

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

    std::vector<Eigen::MatrixXd> qs(num);
    std::vector<Eigen::MatrixXd> ys(num);
    std::vector<Eigen::MatrixXd > vs(num);
    std::vector<Eigen::MatrixXd > ds(num);
    std::vector<Eigen::MatrixXd > gs(num);
    std::vector<Eigen::MatrixXd > cs(num);
    std::vector<Eigen::MatrixXd > ts(num);

//#pragma omp parallel
    {
//#pragma omp for
        for(int j = 0; j < num; ++j) {
            std::vector<Eigen::Vector3d> rr(rr0);
            std::vector<Eigen::Vector3d> bb(bb0);
            std::vector<Eigen::Vector3d> nv(nv0);
            for (int i = 0; i < rr.size(); ++i) {
                rr[i] += noise * Eigen::Vector3d::Random();
                bb[i] += noise * Eigen::Vector3d::Random();
                nv[i] += noise * Eigen::Vector3d::Random();
            }

            Eigen::Matrix<double, 4, 64> W;
            Eigen::Matrix<double, 4, 4> Q;
            Eigen::Matrix<double, 3, 28> D;
            Eigen::Matrix<double, 3, 9> G;
            Eigen::Vector3d c;
            Eigen::Matrix<double, 4, 24> coef_f_q_sym;
            Eigen::Matrix<double, 1, 85> coef_J_pure;
            Eigen::Matrix<double, 3, 11> coefs_tq;
            Eigen::Matrix<double, 3, 3> pinvG;
            pTop_WQD(W, Q, D, G, c, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG, rr, bb, nv);

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
            qs[j].resize(4, 1);
            qs[j] = q0;
            ts[j].resize(3, 1);
            ts[j] = t;
            Eigen::Matrix<double, 28, 1> v;
            Eigen::Matrix<double, 9, 1> y;
            v.setZero();
            y.setZero();
            v_func_pTop(v, q0);
            y_func_pTop(y, q0);
            vs[j].resize(28, 1);
            vs[j] = v;
            ys[j].resize(9, 1);
            ys[j] = y;
            ds[j].resize(3 * 28, 1);
            ds[j] = vec(D);
            gs[j].resize(3 * 9, 1);
            gs[j] = vec(G);
            cs[j].resize(3, 1);
            cs[j] = c;

            if (verbose) {
                std::cout << "True X: " << std::endl << XX << std::endl;
                std::cout << "QPEP X: " << std::endl << X << std::endl << std::endl;
            }
        }
    };

    Eigen::Vector4d mean_q = mean(qs);

//    std::vector<Eigen::MatrixXd> qs_(qs);
//    std::for_each(qs_.begin(), qs_.end(), [=](Eigen::Vector4d& x){x = x - mean_q;});
//    Eigen::Matrix4d res(Eigen::Matrix4d::Zero());
//    std::for_each(qs_.begin(), qs_.end(), [&res](const Eigen::Vector4d& x){res += x * x.transpose();});
//    Eigen::Matrix4d Sigma_data_stat_ = res * 1.0 / ((double) qs_.size());
    Eigen::Matrix4d Sigma_q_stat = covx(qs, qs);
    Eigen::MatrixXd Sigma_t_stat = covx(ts, ts);
    Eigen::MatrixXd Sigma_v_stat = covx(vs, vs);
    Eigen::MatrixXd Sigma_y_stat = covx(ys, ys);
    Eigen::MatrixXd Sigma_d_stat = covx(ds, ds);
    Eigen::MatrixXd Sigma_g_stat = covx(gs, gs);
    Eigen::MatrixXd Sigma_c_stat = covx(cs, cs);
    Eigen::MatrixXd Sigma_gc_stat = covx(gs, cs);
    Eigen::MatrixXd Sigma_dg_stat = covx(ds, gs);
    Eigen::MatrixXd Sigma_dc_stat = covx(ds, cs);
    Eigen::MatrixXd Sigma_cg_stat = Sigma_gc_stat.transpose();
    Eigen::MatrixXd Sigma_gd_stat = Sigma_dg_stat.transpose();
    Eigen::MatrixXd Sigma_cd_stat = Sigma_dc_stat.transpose();




    std::vector<Eigen::Vector3d> rr(rr0);
    std::vector<Eigen::Vector3d> bb(bb0);
    std::vector<Eigen::Vector3d> nv(nv0);
    for (int i = 0; i < rr.size(); ++i) {
        rr[i] += noise * Eigen::Vector3d::Random();
        bb[i] += noise * Eigen::Vector3d::Random();
        nv[i] += noise * Eigen::Vector3d::Random();
    }

    Eigen::Matrix<double, 4, 64> W;
    Eigen::Matrix<double, 4, 4> Q;
    Eigen::Matrix<double, 3, 28> D;
    Eigen::Matrix<double, 3, 9> G;
    Eigen::Vector3d c;
    Eigen::Matrix<double, 4, 24> coef_f_q_sym;
    Eigen::Matrix<double, 1, 85> coef_J_pure;
    Eigen::Matrix<double, 3, 11> coefs_tq;
    Eigen::Matrix<double, 3, 3> pinvG;
    pTop_WQD(W, Q, D, G, c, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG, rr, bb, nv);

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
    Eigen::Vector4d q = R2q(R);
    Eigen::Matrix<double, 28, 1> v;
    Eigen::Matrix<double, 9, 1> y;
    v.setZero();
    y.setZero();
    v_func_pTop(v, q);
    y_func_pTop(y, q);
    Eigen::MatrixXd V_cal = Eigen::kroneckerProduct(v.transpose(), Eigen::Matrix3d::Identity());
    Eigen::MatrixXd Y_cal = Eigen::kroneckerProduct(y.transpose(), Eigen::Matrix3d::Identity());
    Eigen::MatrixXd partial_v_q_sym_val, partial_y_q_sym_val;
    partial_v_q_sym_val.resize(28, 4);
    partial_y_q_sym_val.resize(9, 4);
    partial_v_q_sym_val.setZero();
    partial_y_q_sym_val.setZero();
    jv_func_pTop(partial_v_q_sym_val, q);
    jy_func_pTop(partial_y_q_sym_val, q);
    Eigen::MatrixXd F = D * partial_v_q_sym_val + G * partial_y_q_sym_val;
    Eigen::MatrixXd cov_left = V_cal * Sigma_d_stat * V_cal.transpose() + Y_cal * Sigma_g_stat * Y_cal.transpose() + Sigma_c_stat +
                               V_cal * Sigma_dg_stat * Y_cal.transpose() + V_cal * Sigma_dc_stat + Y_cal * Sigma_gc_stat +
                               Sigma_cd_stat * V_cal.transpose() + Sigma_cg_stat * Y_cal.transpose() + Y_cal * Sigma_gd_stat * V_cal.transpose();

    Eigen::Matrix4d cov;

    Eigen::MatrixXd cov_left_ = cov_left;
    double scaling = 1.0;
    double cov_tr = cov_left_.trace();
//    if(cov_tr <= 1e-5)
//    {
        scaling = 3e-5 / (cov_tr);
        cov_left_ *= scaling;
//    }
    csdp_cov(cov, F, cov_left_, q);
    cov = cov / scaling;
    std::cout << "Estimated Covariance:" << std::endl << cov << std::endl;

    std::vector<Eigen::Vector4d> qs__(num);
    for(int i = 0; i < num; ++i)
    {
        qs__[i] = qs[i];
    }
    double scale = fabs(cov_left.trace() / (F * cov * F.transpose()).trace());
    plotQuatCov(img, Sigma_q_stat, scale * cov, qs__, mean_q);

    std::cout << "Stat Covariance:" << std::endl << Sigma_q_stat << std::endl;
}
#endif

int main(int argc,char ** argv) {

    double time = 0.0;
    clock_t time1 = clock(), time2;
    double loops = 100.0;

//    Eigen::Vector4d q;
//    q << 0.319153034662463, 0.0137947360154964, 0.906900328581383, 0.274741405221315;
//    Eigen::Matrix3d cov_left;
//    cov_left << 4.32090178480567e-09, 1.9881316485166e-09, 7.57378543793111e-10,
//                1.98813164851657e-09, 5.7051245100625e-08, 2.83263963656342e-08,
//                7.57378543793107e-10, 2.83263963656342e-08, 1.43515132711022e-08;
//    Eigen::Matrix<double, 3, 4> F;
//    F << -0.0133788091155435,       -0.0319461589550735,       -0.0223229328790909,      -0.00211035440564425,
//         -0.150273289955035,      0.00919882917604315,        -0.657044040069996,        -0.239866632558646,
//         -0.0175483625464778,       0.00864823916809603,        -0.160369725637978,       -0.0703673441123501;
//
//    Eigen::Matrix4d cov;
//    time1 = clock();
//    for(int i = 0; i < loops; ++i)
//        csdp_cov(cov, F, cov_left, q);
//    time2 = clock();
//    std::cout << "Covariance: " << cov << std::endl;
//    time = time2 - time1;
//    std::cout << "CSDP Time: " << time / loops / double(CLOCKS_PER_SEC) << std::endl;


#ifdef USE_OPENCV
    const int row = 1920;
    const int col = 1920;
    cv::Mat imageDraw = cv::Mat::zeros(row, col, CV_8UC3);
    cv::Mat ColorMask(row, col, CV_8UC3, cv::Scalar(1, 1, 1) * 255);
    cv::addWeighted(imageDraw, 0.0, ColorMask, 1.0, 0, imageDraw);

    const bool verbose = false;
    test_pTop_noise("../data/pTop_data-100pt-1.txt",
                    imageDraw, 1500, 1e-5, verbose);

    imshow("imageDraw", imageDraw);
    cv::waitKey(0);
#endif


    std::cout.precision(16);
//    test_generateProjectedPoints();

    time1 = clock();
    loops = 10.0;

    for(int i = 0; i < (int) loops; ++i)
    {
//        test_pnp_WQD("../data/pTop_data-5pt-1.txt", false);
        test_pTop_WQD("../data/pTop_data-100pt-1.txt", false);
    }
    time2 = clock();
    time = time2 - time1;
    std::cout << "Time: " << time / loops / double(CLOCKS_PER_SEC) << std::endl;


    return 0;
}
