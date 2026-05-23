// LibQPEP: A Library for Globally Optimal Solving Quadratic Pose Estimation Problems,
//          It also gives highly accurate uncertainty description of the solutions.
//
// Author: Jin Wu
// Affiliation: Hong Kong University of Science and Technology (HKUST)
// Emails: jin_wu_uestc@hotmail.com; jwucp@connect.ust.hk
// Reference: Wu, J., et al. (2022) Quadratic Pose Estimation Problems: 
//                                  Globally Optimal Solutions, 
//                                  Solvability/Observability Analysis,
//                                  and Uncertainty Description.
//                                  IEEE Transactions on Robotics.
//                                  https://doi.org/10.1109/TRO.2022.3155880
//
//
// CovEstimation.cpp: Covariance estimation module


#include "CovEstimation.h"
#include "declarations.h"

void csdp_cov(Eigen::Matrix4d& cov,
              const Eigen::MatrixXd& F,
              const Eigen::Matrix3d& cov_left,
              const Eigen::Vector4d& q)
{
    struct blockmatrix C;
    double *b;
    struct constraintmatrix *constraints;
    double q1 = q(0);
    double q2 = q(1);
    double q3 = q(2);
    double q4 = q(3);

    struct blockmatrix X, Z;
    double *y;
    double pobj, dobj;

    struct sparseblock *blockptr;
    int ret;

    C.nblocks = 2;
    C.blocks = (struct blockrec *) malloc(3 * sizeof(struct blockrec));

    C.blocks[1].blockcategory = MATRIX;
    C.blocks[1].blocksize = 4;
    C.blocks[1].data.mat = (double *) malloc((4 * 4) * sizeof(double));
    memset(C.blocks[1].data.mat, 0, (4 * 4) * sizeof(double));

    C.blocks[2].blockcategory = DIAG;
    C.blocks[2].blocksize = 12;
    C.blocks[2].data.vec = (double *) malloc((12 + 1) * sizeof(double));

    C.blocks[2].data.vec[1] = cov_left(0, 0);
    C.blocks[2].data.vec[2] = cov_left(0, 1);
    C.blocks[2].data.vec[3] = cov_left(0, 2);
    C.blocks[2].data.vec[4] = cov_left(1, 1);
    C.blocks[2].data.vec[5] = cov_left(1, 2);
    C.blocks[2].data.vec[6] = cov_left(2, 2);
    C.blocks[2].data.vec[7] = - cov_left(0, 0);
    C.blocks[2].data.vec[8] = - cov_left(0, 1);
    C.blocks[2].data.vec[9] = - cov_left(0, 2);
    C.blocks[2].data.vec[10] = - cov_left(1, 1);
    C.blocks[2].data.vec[11] = - cov_left(1, 2);
    C.blocks[2].data.vec[12] = - cov_left(2, 2);

    double FF1_1 = F(0, 0);
    double FF1_2 = F(0, 1);
    double FF1_3 = F(0, 2);
    double FF1_4 = F(0, 3);
    double FF2_1 = F(1, 0);
    double FF2_2 = F(1, 1);
    double FF2_3 = F(1, 2);
    double FF2_4 = F(1, 3);
    double FF3_1 = F(2, 0);
    double FF3_2 = F(2, 1);
    double FF3_3 = F(2, 2);
    double FF3_4 = F(2, 3);

    Eigen::Matrix<double, 6, 10> H;
    H(0, 0) = FF1_1*FF1_1;
    H(0, 1) = FF1_1*FF1_2*2.0;
    H(0, 2) = FF1_1*FF1_3*2.0;
    H(0, 3) = FF1_1*FF1_4*2.0;
    H(0, 4) = FF1_2*FF1_2;
    H(0, 5) = FF1_2*FF1_3*2.0;
    H(0, 6) = FF1_2*FF1_4*2.0;
    H(0, 7) = FF1_3*FF1_3;
    H(0, 8) = FF1_3*FF1_4*2.0;
    H(0, 9) = FF1_4*FF1_4;
    H(1, 0) = FF1_1*FF2_1;
    H(1, 1) = FF1_1*FF2_2+FF1_2*FF2_1;
    H(1, 2) = FF1_1*FF2_3+FF1_3*FF2_1;
    H(1, 3) = FF1_1*FF2_4+FF1_4*FF2_1;
    H(1, 4) = FF1_2*FF2_2;
    H(1, 5) = FF1_2*FF2_3+FF1_3*FF2_2;
    H(1, 6) = FF1_2*FF2_4+FF1_4*FF2_2;
    H(1, 7) = FF1_3*FF2_3;
    H(1, 8) = FF1_3*FF2_4+FF1_4*FF2_3;
    H(1, 9) = FF1_4*FF2_4;
    H(2, 0) = FF1_1*FF3_1;
    H(2, 1) = FF1_1*FF3_2+FF1_2*FF3_1;
    H(2, 2) = FF1_1*FF3_3+FF1_3*FF3_1;
    H(2, 3) = FF1_1*FF3_4+FF1_4*FF3_1;
    H(2, 4) = FF1_2*FF3_2;
    H(2, 5) = FF1_2*FF3_3+FF1_3*FF3_2;
    H(2, 6) = FF1_2*FF3_4+FF1_4*FF3_2;
    H(2, 7) = FF1_3*FF3_3;
    H(2, 8) = FF1_3*FF3_4+FF1_4*FF3_3;
    H(2, 9) = FF1_4*FF3_4;
    H(3, 0) = FF2_1*FF2_1;
    H(3, 1) = FF2_1*FF2_2*2.0;
    H(3, 2) = FF2_1*FF2_3*2.0;
    H(3, 3) = FF2_1*FF2_4*2.0;
    H(3, 4) = FF2_2*FF2_2;
    H(3, 5) = FF2_2*FF2_3*2.0;
    H(3, 6) = FF2_2*FF2_4*2.0;
    H(3, 7) = FF2_3*FF2_3;
    H(3, 8) = FF2_3*FF2_4*2.0;
    H(3, 9) = FF2_4*FF2_4;
    H(4, 0) = FF2_1*FF3_1;
    H(4, 1) = FF2_1*FF3_2+FF2_2*FF3_1;
    H(4, 2) = FF2_1*FF3_3+FF2_3*FF3_1;
    H(4, 3) = FF2_1*FF3_4+FF2_4*FF3_1;
    H(4, 4) = FF2_2*FF3_2;
    H(4, 5) = FF2_2*FF3_3+FF2_3*FF3_2;
    H(4, 6) = FF2_2*FF3_4+FF2_4*FF3_2;
    H(4, 7) = FF2_3*FF3_3;
    H(4, 8) = FF2_3*FF3_4+FF2_4*FF3_3;
    H(4, 9) = FF2_4*FF3_4;
    H(5, 0) = FF3_1*FF3_1;
    H(5, 1) = FF3_1*FF3_2*2.0;
    H(5, 2) = FF3_1*FF3_3*2.0;
    H(5, 3) = FF3_1*FF3_4*2.0;
    H(5, 4) = FF3_2*FF3_2;
    H(5, 5) = FF3_2*FF3_3*2.0;
    H(5, 6) = FF3_2*FF3_4*2.0;
    H(5, 7) = FF3_3*FF3_3;
    H(5, 8) = FF3_3*FF3_4*2.0;
    H(5, 9) = FF3_4*FF3_4;

    b = (double *) malloc((10 + 1) * sizeof(double));
    b[1] = - q1*q1;
    b[2] = - 2.0 * q1*q2;
    b[3] = - 2.0 * q1*q3;
    b[4] = - 2.0 * q1*q4;
    b[5] = - q2*q2;
    b[6] = - 2.0 * q2*q3;
    b[7] = - 2.0 * q2*q4;
    b[8] = - q3*q3;
    b[9] = - 2.0 * q3*q4;
    b[10] = - q4*q4;

    constraints = (struct constraintmatrix *) malloc((10 + 1)*sizeof(struct constraintmatrix));
    int ijlist[10][2] = {
            {1, 1},
            {1, 2},
            {1, 3},
            {1, 4},
            {2, 2},
            {2, 3},
            {2, 4},
            {3, 3},
            {3, 4},
            {4, 4}
    };

    for(int i = 1; i <= 10; ++i)
    {
        constraints[i].blocks = NULL;

        blockptr = (struct sparseblock *) malloc(sizeof(struct sparseblock));
        blockptr->blocknum = 2;
        blockptr->blocksize = 12;
        blockptr->constraintnum = i;
        blockptr->next = NULL;
        blockptr->nextbyblock = NULL;
        blockptr->entries = (double *) malloc((12 + 1) * sizeof(double));
        blockptr->iindices = (int *) malloc((12 + 1) * sizeof(int));
        blockptr->jindices = (int *) malloc((12 + 1) * sizeof(int));
        blockptr->numentries = 12;
        int counter = 1;
        for(int j = 0; j < 6; ++j)
        {
            blockptr->iindices[counter] = j + 1;
            blockptr->jindices[counter] = j + 1;
            blockptr->entries[counter] = - H(j, i - 1);
            counter = counter + 1;
        }
        for(int j = 6; j < 12; ++j)
        {
            blockptr->iindices[counter] = j + 1;
            blockptr->jindices[counter] = j + 1;
            blockptr->entries[counter] = H(j - 6, i - 1);
            counter = counter + 1;
        }

        blockptr->next = constraints[i].blocks;
        constraints[i].blocks = blockptr;

        blockptr = (struct sparseblock *) malloc(sizeof(struct sparseblock));

        blockptr->blocknum = 1;
        blockptr->blocksize = 4;
        blockptr->constraintnum = i;
        blockptr->next = NULL;
        blockptr->nextbyblock = NULL;
        blockptr->entries = (double *) malloc((1 + 1) * sizeof(double));
        blockptr->iindices = (int *) malloc((1 + 1) * sizeof(int));
        blockptr->jindices = (int *) malloc((1 + 1) * sizeof(int));
        blockptr->numentries = 1;
        blockptr->iindices[1] = ijlist[i - 1][0];
        blockptr->jindices[1] = ijlist[i - 1][1];
        blockptr->entries[1] = - 1.0;
        blockptr->next = constraints[i].blocks;
        constraints[i].blocks = blockptr;
    }

//    write_prob("prob.dat-s", 4, 10, C, b, constraints);

    initsoln(16, 10, C, b, constraints, &X, &y, &Z);
    ret = easy_sdp(16,10, C, b, constraints, 0.0, &X, &y, &Z, &pobj, &dobj);

//    if (ret == 0)
//        printf("The objective value is %.7e \n", (dobj + pobj) / 2);
//    else
//        printf("SDP failed.\n");

    cov(0, 0) = Z.blocks[1].data.mat[ijtok(1, 1, Z.blocks[1].blocksize)];
    cov(0, 1) = Z.blocks[1].data.mat[ijtok(1, 2, Z.blocks[1].blocksize)];
    cov(0, 2) = Z.blocks[1].data.mat[ijtok(1, 3, Z.blocks[1].blocksize)];
    cov(0, 3) = Z.blocks[1].data.mat[ijtok(1, 4, Z.blocks[1].blocksize)];
    cov(1, 1) = Z.blocks[1].data.mat[ijtok(2, 2, Z.blocks[1].blocksize)];
    cov(1, 2) = Z.blocks[1].data.mat[ijtok(2, 3, Z.blocks[1].blocksize)];
    cov(1, 3) = Z.blocks[1].data.mat[ijtok(2, 4, Z.blocks[1].blocksize)];
    cov(2, 2) = Z.blocks[1].data.mat[ijtok(3, 3, Z.blocks[1].blocksize)];
    cov(2, 3) = Z.blocks[1].data.mat[ijtok(3, 4, Z.blocks[1].blocksize)];
    cov(3, 3) = Z.blocks[1].data.mat[ijtok(4, 4, Z.blocks[1].blocksize)];
    for(int i = 0; i < 4; ++i)
        for(int j = i + 1; j < 4; ++j)
            cov(j, i) = cov(i, j);

//    write_sol("prob.sol", 16, 10, X, y, Z);
    free_prob(16, 10, C, b, constraints, X, y, Z);
}
