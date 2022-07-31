#include "hand_eye_small_rotation.h"
#include "utils.h"



Eigen::Matrix<double, 6, 6> AA_hand_eye_func(const Eigen::Matrix4d& A, const Eigen::Matrix4d& B)
{
    double A11 = A(0, 0);
    double A12 = A(0, 1);
    double A13 = A(0, 2);
    double A14 = A(0, 3);
    double A21 = A(1, 0);
    double A22 = A(1, 1);
    double A23 = A(1, 2);
    double A24 = A(1, 3);
    double A31 = A(2, 0);
    double A32 = A(2, 1);
    double A33 = A(2, 2);
    double A34 = A(2, 3);
    double B11 = B(0, 0);
    double B12 = B(0, 1);
    double B13 = B(0, 2);
    double B14 = B(0, 3);
    double B21 = B(1, 0);
    double B22 = B(1, 1);
    double B23 = B(1, 2);
    double B24 = B(1, 3);
    double B31 = B(2, 0);
    double B32 = B(2, 1);
    double B33 = B(2, 2);
    double B34 = B(2, 3);

    double symobj[6][6];

    std::memset(symobj, 0, 6 * 6 * sizeof(double));
    symobj[0][0] = (A12*A12)*2.0+(A13*A13)*2.0+(A22*A22)*2.0+(A23*A23)*2.0+(A32*A32)*2.0+(A33*A33)*2.0+(B21*B21)*2.0+(B22*B22)*2.0+(B23*B23)*2.0+(B24*B24)*2.0+(B31*B31)*2.0+(B32*B32)*2.0+(B33*B33)*2.0+(B34*B34)*2.0-A22*B33*4.0+A23*B32*4.0+A32*B23*4.0-A33*B22*4.0;
    symobj[0][1] = A11*A12*-2.0-A21*A22*2.0-A31*A32*2.0+A12*B33*2.0-A13*B32*2.0-A32*B13*2.0+A33*B12*2.0+A21*B33*2.0-A23*B31*2.0-A31*B23*2.0+A33*B21*2.0-B11*B21*2.0-B12*B22*2.0-B13*B23*2.0-B14*B24*2.0;
    symobj[0][2] = A11*A13*-2.0-A21*A23*2.0-A31*A33*2.0-A12*B23*2.0+A13*B22*2.0+A22*B13*2.0-A23*B12*2.0-A21*B32*2.0+A22*B31*2.0+A31*B22*2.0-A32*B21*2.0-B11*B31*2.0-B12*B32*2.0-B13*B33*2.0-B14*B34*2.0;
    symobj[0][3] = A21*B34*2.0-A31*B24*2.0;
    symobj[0][4] = B34*-2.0+A22*B34*2.0-A32*B24*2.0;
    symobj[0][5] = B24*2.0+A23*B34*2.0-A33*B24*2.0;
    symobj[1][0] = A11*A12*-2.0-A21*A22*2.0-A31*A32*2.0+A12*B33*2.0-A13*B32*2.0-A32*B13*2.0+A33*B12*2.0+A21*B33*2.0-A23*B31*2.0-A31*B23*2.0+A33*B21*2.0-B11*B21*2.0-B12*B22*2.0-B13*B23*2.0-B14*B24*2.0;
    symobj[1][1] = (A11*A11)*2.0+(A13*A13)*2.0+(A21*A21)*2.0+(A23*A23)*2.0+(A31*A31)*2.0+(A33*A33)*2.0+(B11*B11)*2.0+(B12*B12)*2.0+(B13*B13)*2.0+(B14*B14)*2.0+(B31*B31)*2.0+(B32*B32)*2.0+(B33*B33)*2.0+(B34*B34)*2.0-A11*B33*4.0+A13*B31*4.0+A31*B13*4.0-A33*B11*4.0;
    symobj[1][2] = A12*A13*-2.0-A22*A23*2.0-A32*A33*2.0+A11*B23*2.0-A13*B21*2.0-A21*B13*2.0+A23*B11*2.0+A11*B32*2.0-A12*B31*2.0-A31*B12*2.0+A32*B11*2.0-B21*B31*2.0-B22*B32*2.0-B23*B33*2.0-B24*B34*2.0;
    symobj[1][3] = B34*2.0-A11*B34*2.0+A31*B14*2.0;
    symobj[1][4] = A12*B34*-2.0+A32*B14*2.0;
    symobj[1][5] = B14*-2.0-A13*B34*2.0+A33*B14*2.0;
    symobj[2][0] = A11*A13*-2.0-A21*A23*2.0-A31*A33*2.0-A12*B23*2.0+A13*B22*2.0+A22*B13*2.0-A23*B12*2.0-A21*B32*2.0+A22*B31*2.0+A31*B22*2.0-A32*B21*2.0-B11*B31*2.0-B12*B32*2.0-B13*B33*2.0-B14*B34*2.0;
    symobj[2][1] = A12*A13*-2.0-A22*A23*2.0-A32*A33*2.0+A11*B23*2.0-A13*B21*2.0-A21*B13*2.0+A23*B11*2.0+A11*B32*2.0-A12*B31*2.0-A31*B12*2.0+A32*B11*2.0-B21*B31*2.0-B22*B32*2.0-B23*B33*2.0-B24*B34*2.0;
    symobj[2][2] = (A11*A11)*2.0+(A12*A12)*2.0+(A21*A21)*2.0+(A22*A22)*2.0+(A31*A31)*2.0+(A32*A32)*2.0+(B11*B11)*2.0+(B12*B12)*2.0+(B13*B13)*2.0+(B14*B14)*2.0+(B21*B21)*2.0+(B22*B22)*2.0+(B23*B23)*2.0+(B24*B24)*2.0-A11*B22*4.0+A12*B21*4.0+A21*B12*4.0-A22*B11*4.0;
    symobj[2][3] = B24*-2.0+A11*B24*2.0-A21*B14*2.0;
    symobj[2][4] = B14*2.0+A12*B24*2.0-A22*B14*2.0;
    symobj[2][5] = A13*B24*2.0-A23*B14*2.0;
    symobj[3][0] = A21*B34*2.0-A31*B24*2.0;
    symobj[3][1] = B34*2.0-A11*B34*2.0+A31*B14*2.0;
    symobj[3][2] = B24*-2.0+A11*B24*2.0-A21*B14*2.0;
    symobj[3][3] = A11*-4.0+(A11*A11)*2.0+(A21*A21)*2.0+(A31*A31)*2.0+2.0;
    symobj[3][4] = A12*-2.0-A21*2.0+A11*A12*2.0+A21*A22*2.0+A31*A32*2.0;
    symobj[3][5] = A13*-2.0-A31*2.0+A11*A13*2.0+A21*A23*2.0+A31*A33*2.0;
    symobj[4][0] = B34*-2.0+A22*B34*2.0-A32*B24*2.0;
    symobj[4][1] = A12*B34*-2.0+A32*B14*2.0;
    symobj[4][2] = B14*2.0+A12*B24*2.0-A22*B14*2.0;
    symobj[4][3] = A12*-2.0-A21*2.0+A11*A12*2.0+A21*A22*2.0+A31*A32*2.0;
    symobj[4][4] = A22*-4.0+(A12*A12)*2.0+(A22*A22)*2.0+(A32*A32)*2.0+2.0;
    symobj[4][5] = A23*-2.0-A32*2.0+A12*A13*2.0+A22*A23*2.0+A32*A33*2.0;
    symobj[5][0] = B24*2.0+A23*B34*2.0-A33*B24*2.0;
    symobj[5][1] = B14*-2.0-A13*B34*2.0+A33*B14*2.0;
    symobj[5][2] = A13*B24*2.0-A23*B14*2.0;
    symobj[5][3] = A13*-2.0-A31*2.0+A11*A13*2.0+A21*A23*2.0+A31*A33*2.0;
    symobj[5][4] = A23*-2.0-A32*2.0+A12*A13*2.0+A22*A23*2.0+A32*A33*2.0;
    symobj[5][5] = A33*-4.0+(A13*A13)*2.0+(A23*A23)*2.0+(A33*A33)*2.0+2.0;

    Eigen::Matrix<double, 6, 6> tmp;
    for(int i = 0; i < 6; ++i)
        for(int j = 0; j < 6; ++j)
            tmp(i, j) = symobj[i][j];
    return tmp;
}



Eigen::Matrix<double, 6, 1> bb_hand_eye_func(const Eigen::Matrix4d& A, const Eigen::Matrix4d& B)
{
    double A11 = A(0, 0);
    double A12 = A(0, 1);
    double A13 = A(0, 2);
    double A14 = A(0, 3);
    double A21 = A(1, 0);
    double A22 = A(1, 1);
    double A23 = A(1, 2);
    double A24 = A(1, 3);
    double A31 = A(2, 0);
    double A32 = A(2, 1);
    double A33 = A(2, 2);
    double A34 = A(2, 3);
    double B11 = B(0, 0);
    double B12 = B(0, 1);
    double B13 = B(0, 2);
    double B14 = B(0, 3);
    double B21 = B(1, 0);
    double B22 = B(1, 1);
    double B23 = B(1, 2);
    double B24 = B(1, 3);
    double B31 = B(2, 0);
    double B32 = B(2, 1);
    double B33 = B(2, 2);
    double B34 = B(2, 3);

    double symobj[6][1];

    std::memset(symobj, 0, 6 * 1 * sizeof(double));
    symobj[0][0] = A12*B13*-2.0+A13*B12*2.0-A22*B23*2.0+A23*B22*2.0-A21*B31*2.0+A31*B21*2.0-A22*B32*2.0+A32*B22*2.0-A23*B33*2.0+A33*B23*2.0-A24*B34*2.0+A34*B24*2.0-A32*B33*2.0+A33*B32*2.0;
    symobj[1][0] = A11*B13*2.0-A13*B11*2.0+A11*B31*2.0-A31*B11*2.0+A12*B32*2.0+A21*B23*2.0-A23*B21*2.0-A32*B12*2.0+A13*B33*2.0-A33*B13*2.0+A14*B34*2.0-A34*B14*2.0+A31*B33*2.0-A33*B31*2.0;
    symobj[2][0] = A11*B12*-2.0+A12*B11*2.0-A11*B21*2.0+A21*B11*2.0-A12*B22*2.0+A22*B12*2.0-A13*B23*2.0+A23*B13*2.0-A14*B24*2.0+A24*B14*2.0-A21*B22*2.0+A22*B21*2.0-A31*B32*2.0+A32*B31*2.0;
    symobj[3][0] = A14*2.0-B14*2.0-A11*A14*2.0-A21*A24*2.0-A31*A34*2.0+A11*B14*2.0+A21*B24*2.0+A31*B34*2.0;
    symobj[4][0] = A24*2.0-B24*2.0-A12*A14*2.0-A22*A24*2.0-A32*A34*2.0+A12*B14*2.0+A22*B24*2.0+A32*B34*2.0;
    symobj[5][0] = A34*2.0-B34*2.0-A13*A14*2.0-A23*A24*2.0-A33*A34*2.0+A13*B14*2.0+A23*B24*2.0+A33*B34*2.0;

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

void hand_eye_small_rotation(Eigen::Matrix4d& X,
                             const std::vector<Eigen::Matrix4d>& A,
                             const std::vector<Eigen::Matrix4d>& B)
{
    Eigen::Matrix<double, 6, 6> AA;
    Eigen::Matrix<double, 6, 1> bb;
    AA.setZero();
    bb.setZero();

    for(int i = 0; i < A.size(); ++i)
    {
        AA += AA_hand_eye_func(A[i], B[i]);
        bb += bb_hand_eye_func(A[i], B[i]);
    }

    Eigen::Matrix<double, 6, 1> xi = AA.inverse() * bb;
    Eigen::Matrix3d R = orthonormalize(Eigen::Matrix3d::Identity() + skew(xi.topRows(3)));
    Eigen::Vector3d t = xi.bottomRows(3);
    X << R, t, Eigen::Vector3d::Zero(3).transpose(), 1.0;
}


// AX = XB
// A * X * X0 = X * X0 * B
// A * X = X * C
// C = X0 * B * inv(X0)
void hand_eye_small_rotation_refine(Eigen::Matrix4d& X,
                             const std::vector<Eigen::Matrix4d>& A,
                             const std::vector<Eigen::Matrix4d>& B,
                             const Eigen::Matrix4d& X0,
                             const int& num)
{
    std::vector<Eigen::Matrix4d> BB;
    BB.resize(B.size());
    Eigen::Matrix4d X0_ = X0;

    for(int j = 0; j < num; ++j)
    {
        X = X0_;
        for(int i = 0; i < B.size(); ++i) {
            BB[i] = X * B[i] * X.inverse();
        }
        hand_eye_small_rotation(X, A, BB);
        X0_ = X * X0_;
    }
    X = X0_;
}