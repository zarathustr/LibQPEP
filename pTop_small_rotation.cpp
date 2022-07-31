#include "pTop_small_rotation.h"
#include "utils.h"

Eigen::Matrix<double, 6, 6> AA_pTop_func(const Eigen::Vector3d &rr, const Eigen::Vector3d &bb, const Eigen::Vector3d &nv)
{
    double r1 = rr(0);
    double r2 = rr(1);
    double r3 = rr(2);
    double b1 = bb(0);
    double b2 = bb(1);
    double b3 = bb(2);
    double n1 = nv(0);
    double n2 = nv(1);
    double n3 = nv(2);

    double symobj[6][6];

    std::memset(symobj, 0, 6 * 6 * sizeof(double));
    symobj[0][0] = (n2*n2)*(r3*r3)*2.0+(n3*n3)*(r2*r2)*2.0-n2*n3*r2*r3*4.0;
    symobj[0][1] = n1*n2*(r3*r3)*-2.0-(n3*n3)*r1*r2*2.0+n1*n3*r2*r3*2.0+n2*n3*r1*r3*2.0;
    symobj[0][2] = n1*n3*(r2*r2)*-2.0-(n2*n2)*r1*r3*2.0+n1*n2*r2*r3*2.0+n2*n3*r1*r2*2.0;
    symobj[0][3] = n1*n2*r3*-2.0+n1*n3*r2*2.0;
    symobj[0][4] = (n2*n2)*r3*-2.0+n2*n3*r2*2.0;
    symobj[0][5] = (n3*n3)*r2*2.0-n2*n3*r3*2.0;
    symobj[1][0] = n1*n2*(r3*r3)*-2.0-(n3*n3)*r1*r2*2.0+n1*n3*r2*r3*2.0+n2*n3*r1*r3*2.0;
    symobj[1][1] = (n1*n1)*(r3*r3)*2.0+(n3*n3)*(r1*r1)*2.0-n1*n3*r1*r3*4.0;
    symobj[1][2] = n2*n3*(r1*r1)*-2.0-(n1*n1)*r2*r3*2.0+n1*n2*r1*r3*2.0+n1*n3*r1*r2*2.0;
    symobj[1][3] = (n1*n1)*r3*2.0-n1*n3*r1*2.0;
    symobj[1][4] = n1*n2*r3*2.0-n2*n3*r1*2.0;
    symobj[1][5] = (n3*n3)*r1*-2.0+n1*n3*r3*2.0;
    symobj[2][0] = n1*n3*(r2*r2)*-2.0-(n2*n2)*r1*r3*2.0+n1*n2*r2*r3*2.0+n2*n3*r1*r2*2.0;
    symobj[2][1] = n2*n3*(r1*r1)*-2.0-(n1*n1)*r2*r3*2.0+n1*n2*r1*r3*2.0+n1*n3*r1*r2*2.0;
    symobj[2][2] = (n1*n1)*(r2*r2)*2.0+(n2*n2)*(r1*r1)*2.0-n1*n2*r1*r2*4.0;
    symobj[2][3] = (n1*n1)*r2*-2.0+n1*n2*r1*2.0;
    symobj[2][4] = (n2*n2)*r1*2.0-n1*n2*r2*2.0;
    symobj[2][5] = n1*n3*r2*-2.0+n2*n3*r1*2.0;
    symobj[3][0] = n1*n2*r3*-2.0+n1*n3*r2*2.0;
    symobj[3][1] = (n1*n1)*r3*2.0-n1*n3*r1*2.0;
    symobj[3][2] = (n1*n1)*r2*-2.0+n1*n2*r1*2.0;
    symobj[3][3] = (n1*n1)*2.0;
    symobj[3][4] = n1*n2*2.0;
    symobj[3][5] = n1*n3*2.0;
    symobj[4][0] = (n2*n2)*r3*-2.0+n2*n3*r2*2.0;
    symobj[4][1] = n1*n2*r3*2.0-n2*n3*r1*2.0;
    symobj[4][2] = (n2*n2)*r1*2.0-n1*n2*r2*2.0;
    symobj[4][3] = n1*n2*2.0;
    symobj[4][4] = (n2*n2)*2.0;
    symobj[4][5] = n2*n3*2.0;
    symobj[5][0] = (n3*n3)*r2*2.0-n2*n3*r3*2.0;
    symobj[5][1] = (n3*n3)*r1*-2.0+n1*n3*r3*2.0;
    symobj[5][2] = n1*n3*r2*-2.0+n2*n3*r1*2.0;
    symobj[5][3] = n1*n3*2.0;
    symobj[5][4] = n2*n3*2.0;
    symobj[5][5] = (n3*n3)*2.0;

    Eigen::Matrix<double, 6, 6> tmp;
    for(int i = 0; i < 6; ++i)
        for(int j = 0; j < 6; ++j)
            tmp(i, j) = symobj[i][j];
    return tmp;
}



Eigen::Matrix<double, 6, 1> bb_pTop_func(const Eigen::Vector3d &rr, const Eigen::Vector3d &bb, const Eigen::Vector3d &nv)
{
    double r1 = rr(0);
    double r2 = rr(1);
    double r3 = rr(2);
    double b1 = bb(0);
    double b2 = bb(1);
    double b3 = bb(2);
    double n1 = nv(0);
    double n2 = nv(1);
    double n3 = nv(2);

    double symobj[6][1];

    std::memset(symobj, 0, 6 * 1 * sizeof(double));
    symobj[0][0] = b2*(n2*n2)*r3*-2.0+b3*(n3*n3)*r2*2.0-n2*n3*(r2*r2)*2.0+n2*n3*(r3*r3)*2.0+(n2*n2)*r2*r3*2.0-(n3*n3)*r2*r3*2.0-b1*n1*n2*r3*2.0+b1*n1*n3*r2*2.0+b2*n2*n3*r2*2.0-b3*n2*n3*r3*2.0+n1*n2*r1*r3*2.0-n1*n3*r1*r2*2.0;
    symobj[1][0] = b1*(n1*n1)*r3*2.0-b3*(n3*n3)*r1*2.0+n1*n3*(r1*r1)*2.0-n1*n3*(r3*r3)*2.0-(n1*n1)*r1*r3*2.0+(n3*n3)*r1*r3*2.0-b1*n1*n3*r1*2.0+b2*n1*n2*r3*2.0-b2*n2*n3*r1*2.0+b3*n1*n3*r3*2.0-n1*n2*r2*r3*2.0+n2*n3*r1*r2*2.0;
    symobj[2][0] = b1*(n1*n1)*r2*-2.0+b2*(n2*n2)*r1*2.0-n1*n2*(r1*r1)*2.0+n1*n2*(r2*r2)*2.0+(n1*n1)*r1*r2*2.0-(n2*n2)*r1*r2*2.0+b1*n1*n2*r1*2.0-b2*n1*n2*r2*2.0-b3*n1*n3*r2*2.0+b3*n2*n3*r1*2.0+n1*n3*r2*r3*2.0-n2*n3*r1*r3*2.0;
    symobj[3][0] = b1*(n1*n1)*2.0-(n1*n1)*r1*2.0+b2*n1*n2*2.0+b3*n1*n3*2.0-n1*n2*r2*2.0-n1*n3*r3*2.0;
    symobj[4][0] = b2*(n2*n2)*2.0-(n2*n2)*r2*2.0+b1*n1*n2*2.0+b3*n2*n3*2.0-n1*n2*r1*2.0-n2*n3*r3*2.0;
    symobj[5][0] = b3*(n3*n3)*2.0-(n3*n3)*r3*2.0+b1*n1*n3*2.0+b2*n2*n3*2.0-n1*n3*r1*2.0-n2*n3*r2*2.0;

    Eigen::Matrix<double, 6, 1> tmp;
    for(int i = 0; i < 6; ++i)
        for(int j = 0; j < 1; ++j)
            tmp(i, j) = symobj[i][j];
    return tmp;
}

static Eigen::Matrix3d skew(Eigen::MatrixXd v)
{
    Eigen::Matrix3d S;
    S << 0,   -v(2),  v(1),
            v(2),  0,    -v(0),
            -v(1), v(0),   0;
    return S;
}

void pTop_small_rotation(Eigen::Matrix4d& X,
                         const std::vector<Eigen::Vector3d> &rr, const std::vector<Eigen::Vector3d> &bb, const std::vector<Eigen::Vector3d> &nv)
{
    Eigen::Matrix<double, 6, 6> AA;
    Eigen::Matrix<double, 6, 1> bbb;
    AA.setZero();
    bbb.setZero();

    for(int i = 0; i < rr.size(); ++i)
    {
        AA += AA_pTop_func(rr[i], bb[i], nv[i]);
        bbb += bb_pTop_func(rr[i], bb[i], nv[i]);
    }

    Eigen::Matrix<double, 6, 1> xi = AA.inverse() * bbb;
    Eigen::Matrix3d R = orthonormalize(Eigen::Matrix3d::Identity() + skew(xi.topRows(3)));
    Eigen::Vector3d t = xi.bottomRows(3);
    X << R, t, Eigen::Vector3d::Zero(3).transpose(), 1.0;
}


void pTop_small_rotation_refine(Eigen::Matrix4d& X,
                                const std::vector<Eigen::Vector3d> &rr, const std::vector<Eigen::Vector3d> &bb, const std::vector<Eigen::Vector3d> &nv,
                                const Eigen::Matrix4d& X0,
                                const int& num)
{
    std::vector<Eigen::Vector3d> RR;
    RR.resize(rr.size());
    Eigen::Matrix4d X0_ = X0;

    for(int j = 0; j < num; ++j)
    {
        X = X0_;
        for(int i = 0; i < RR.size(); ++i) {
            RR[i] = X.topLeftCorner(3, 3) * rr[i] + X.topRightCorner(3, 1);
        }
        pTop_small_rotation(X, RR, bb, nv);
        X0_ = X * X0_;
    }
    X = X0_;
}