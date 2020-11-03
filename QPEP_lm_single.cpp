#include "QPEP_lm_single.h"
#include "utils.h"

Eigen::VectorXd absvec(const Eigen::VectorXd vec)
{
    Eigen::VectorXd t(vec);
    for(int i = 0; i < vec.size(); ++i)
    {
        t(i) = abs(t(i));
    }
    return t;
}

Eigen::Matrix4d inv4(Eigen::Matrix4d A)
{
    double A11 = A(0, 0), A12 = A(0, 1), A13 = A(0, 2), A14 = A(0, 3);
    double A22 = A(1, 1), A23 = A(1, 2), A24 = A(1, 3);
    double A33 = A(2, 2), A34 = A(2, 3);
    double A44 = A(3, 3);
    double t0 = (A13*A13)*(A24*A24)+(A14*A14)*(A23*A23)+(A12*A12)*(A34*A34)-A11*A22*(A34*A34)-A11*(A24*A24)*A33-(A14*A14)*A22*A33-A11*(A23*A23)*A44-(A13*A13)*A22*A44-(A12*A12)*A33*A44-A13*A14*A23*A24*2.0-A12*A13*A24*A34*2.0-A12*A14*A23*A34*2.0+A12*A14*A24*A33*2.0+A13*A14*A22*A34*2.0+A11*A23*A24*A34*2.0+A12*A13*A23*A44*2.0+A11*A22*A33*A44;
    Eigen::Matrix4d T;
    T(0, 0) = -A22*(A34*A34)-(A24*A24)*A33-(A23*A23)*A44+A23*A24*A34*2.0+A22*A33*A44;
    T(0, 1) = A12*(A34*A34)-A13*A24*A34-A14*A23*A34+A14*A24*A33+A13*A23*A44-A12*A33*A44;
    T(0, 2) = A13*(A24*A24)-A14*A23*A24-A12*A24*A34+A14*A22*A34+A12*A23*A44-A13*A22*A44;
    T(0, 3) = A14*(A23*A23)-A13*A23*A24-A12*A23*A34+A12*A24*A33+A13*A22*A34-A14*A22*A33;
    T(1, 1) = -A11*(A34*A34)-(A14*A14)*A33-(A13*A13)*A44+A13*A14*A34*2.0+A11*A33*A44;
    T(1, 2) = (A14*A14)*A23-A13*A14*A24-A12*A14*A34+A11*A24*A34+A12*A13*A44-A11*A23*A44;
    T(1, 3) = (A13*A13)*A24-A13*A14*A23-A12*A13*A34+A12*A14*A33+A11*A23*A34-A11*A24*A33;
    T(2, 2) = -A11*(A24*A24)-(A14*A14)*A22-(A12*A12)*A44+A12*A14*A24*2.0+A11*A22*A44;
    T(2, 3) = (A12*A12)*A34-A12*A13*A24-A12*A14*A23+A13*A14*A22+A11*A23*A24-A11*A22*A34;
    T(3, 3) = -A11*(A23*A23)-(A13*A13)*A22-(A12*A12)*A33+A12*A13*A23*2.0+A11*A22*A33;
    for(int i = 0; i < 4; ++i)
        for(int j = i; j < 4; ++j) {
            T(i, j) = T(i, j) / t0;
            T(j, i) = T(i, j);
        }

    return T;
}

struct QPEP_runtime QPEP_lm_single(Eigen::Matrix3d& R,
                    Eigen::Vector3d& t,
                    Eigen::Matrix4d& X,
                    const Eigen::Vector4d& q0,
                    const int& max_iter,
                    const double& mu,
                    const eq_Jacob_func_handle& eq_Jacob_func,
                    const t_func_handle& t_func,
                    const Eigen::MatrixXd& coef_f_q_sym,
                    const Eigen::MatrixXd& coefs_tq,
                    const Eigen::MatrixXd& pinvG,
                    const struct QPEP_runtime& stat_)
{
    clock_t time1 = clock();
    struct QPEP_runtime stat = stat_;
    double mu_ = mu;
    Eigen::Vector4d qq0(q0);
    Eigen::Vector4d last_q;
    last_q.setZero();
    double last_err = 1e15;
    Eigen::Vector4d residual;
    double err;
    Eigen::Matrix<double, 4, 1> res, grad;
    Eigen::Matrix<double, 4, 4> Jacob, Hess;
    for(int j = 0; j < max_iter; ++j)
    {
        residual = absvec(last_q) - absvec(qq0);
        err = residual.norm();
        if(err < 1e-20)
        {
            break;
        }
        else if(err >= 1e-20 && err < last_err / 2.0)
        {
            mu_ = mu_ / 2.0;
        }
        last_err = err;
        last_q = qq0;
        eq_Jacob_func(res, Jacob, coef_f_q_sym, qq0);
        grad = Jacob.transpose() * res;
        Hess = (Jacob.transpose() * Jacob + mu_ * Eigen::Matrix4d::Identity());
        qq0 = qq0 - inv4(Hess) * grad;
        qq0.normalize();
    }
    Eigen::Quaterniond qq(qq0);
    t = t_func(pinvG, coefs_tq, qq);
    R = q2R(qq);
    X << R, t, Eigen::Vector3d::Zero(3).transpose(), 1.0;
    clock_t time2 = clock();
    stat.timeLM = (time2 - time1) / double(CLOCKS_PER_SEC);
    return stat;
}
