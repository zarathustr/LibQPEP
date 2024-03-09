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

#ifndef LIBQPEP_QPEP_H
#define LIBQPEP_QPEP_H

#include <string>

struct QPEP_options {
    std::string ModuleName = "";
    std::string DecompositionMethod = "";
};

struct QPEP_runtime {
    double timeDecompositionDataPrepare = 1e15;
    double timeDecomposition = 1e15;
    double timeGrobner = 1e15;
    double timeEigen = 1e15;
    double timeLM = 1e15;

    double lossGrobner = 1e15;
    double lossLM = 1e15;

    int statusDecomposition = -1;
    int statusGrobner = -1;
    int statusEigen = -1;
    int statusLM = -1;
};

#endif //LIBQPEP_QPEP_H
