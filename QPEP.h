#ifndef LIBQPEP_QPEP_H
#define LIBQPEP_QPEP_H

#include <string>

struct QPEP_options {
    std::string ModuleName = "";
    std::string DecompositionMethod = "";
};

struct QPEP_runtime {
    double timeDataPrepare = 1e15;
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
