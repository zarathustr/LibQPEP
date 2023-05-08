// LibQPEP: A Library for Globally Optimal Solving Quadratic Pose Estimation Problems,
//          It also gives highly accurate uncertainty description of the solutions.
//
// Authors: Jin Wu and Ming Liu
// Affiliation: Hong Kong University of Science and Technology (HKUST)
// Emails: jin_wu_uestc@hotmail.com; eelium@ust.hk
// Reference: Wu, J., et al. (2022) Quadratic Pose Estimation Problems: 
//                                  Globally Optimal Solutions, 
//                                  Solvability/Observability Analysis,
//                                  and Uncertainty Description.
//                                  IEEE Transactions on Robotics.
//                                  https://doi.org/10.1109/TRO.2022.3155880
//
//
// generateProjectedPoints.cpp: Generate image points for PnP problems


#include "utils.h"
#include <random>
#include <ctime>
#ifdef USE_OPENCL
#include <viennacl/linalg/cg.hpp>
#include <viennacl/linalg/gmres.hpp>
#include <viennacl/linalg/bicgstab.hpp>
#include <viennacl/linalg/ilu.hpp>
#include <viennacl/linalg/jacobi_precond.hpp>
#endif

#ifdef USE_OPENCV

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
#endif

Eigen::Matrix3d orthonormalize(const Eigen::Matrix3d& R)
{
    Eigen::JacobiSVD<Eigen::Matrix3d> solver(R, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector3d singularValues = solver.singularValues();
    Eigen::Matrix3d S;
    S.setIdentity();
    if(singularValues(2) < 0) {
        S(2, 2) = -1.0;
    }
    return solver.matrixU() * S * solver.matrixV().transpose();
}

Eigen::MatrixXd randomMatrix(const int& dim1, const int& dim2, const int& resolution)
{
    static std::default_random_engine e(time(0));
    double res = resolution;
    static std::normal_distribution<double> n(0, res) ;
    Eigen::MatrixXd m = Eigen::MatrixXd::Zero(dim1, dim2).unaryExpr([](double dummy){return n(e);});
    return m * (1.0 / res);
}

QPEP_runtime GaussJordanElimination(
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
           opt.DecompositionMethod == "SuperLU" ||
           opt.DecompositionMethod == "HouseholderQR" ||
           opt.DecompositionMethod == "PartialPivLU" ||
           opt.DecompositionMethod == "SVD" ||
           opt.DecompositionMethod == "BDCSVD" ||
           opt.DecompositionMethod == "Inv" ||
           opt.DecompositionMethod == "Cholesky" ||
           opt.DecompositionMethod == "LinSolve" ||
           opt.DecompositionMethod == "ViennaCL-ConjugateGradient" ||
           opt.DecompositionMethod == "ViennaCL-BiCGStab" ||
           opt.DecompositionMethod == "ViennaCL-GMRES" ||
           opt.DecompositionMethod == "ViennaCL-ILUT" ||
           opt.DecompositionMethod == "ViennaCL-Jacobi");

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
        bool is_symmetric = true;
        for(int i = 0; i < size_GJ; ++i)
        {
            for(int j = i + 1; j < size_GJ; ++j)
            {
                if(CC1(i, j) != CC1(j, i))
                {
                    is_symmetric = false;
                    break;
                }
            }
        }
        clock_t time2 = clock();
        //std::cout << "GJ Elimination " << opt.DecompositionMethod << " data_func time: " << (time2 - time1) / double(CLOCKS_PER_SEC) << std::endl;

        Eigen::PartialPivLU<Eigen::MatrixXd > solver;
        if(is_symmetric)
        {
            solver.compute(CC1);
            C1_ = solver.solve(C2);
        }
        else
        {
            solver.compute(CC1.transpose() * CC1);
            C1_ = solver.solve(CC1.transpose() * C2);
        }
        clock_t time3 = clock();
        stat.timeDecomposition = (time3 - time1) / double(CLOCKS_PER_SEC);
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
#if EIGEN_VERSION_AT_LEAST(3,3,0)
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
#endif
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
    else if(opt.DecompositionMethod == "LinSolve") {
        clock_t time1 = clock();
        Eigen::SparseMatrix<double> C1(size_GJ, size_GJ);
        C1.setZero();
        data_func(C1, C2, data);
        clock_t time2 = clock();
        stat.timeDecompositionDataPrepare = (time2 - time1) / double(CLOCKS_PER_SEC);
        Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
        solver.compute(C1);
        Eigen::SparseMatrix<double> I(size_GJ, size_GJ);
        I.setIdentity();
        C1_ = solver.solve(I).toDense() * C2;
        clock_t time3 = clock();
        stat.timeDecomposition = (time3 - time2) / double(CLOCKS_PER_SEC);
    }
#ifdef USE_SUPERLU
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
#endif
#ifdef USE_OPENCL
    else if(opt.DecompositionMethod == "ViennaCL-ConjugateGradient" ||
            opt.DecompositionMethod == "ViennaCL-BiCGStab" ||
            opt.DecompositionMethod == "ViennaCL-GMRES")
    {
        clock_t time1 = clock();

        Eigen::SparseMatrix<double> C1(size_GJ, size_GJ);
        data_func(C1, C2, data);

        viennacl::compressed_matrix<double> vcl_matrix(size_GJ, size_GJ);
        viennacl::copy(C1, vcl_matrix);
        viennacl::vector<double> vcl_rhs(size_GJ);
        viennacl::vector<double> vcl_result;
        for(int i = 0; i < size_AM; ++i)
        {
            Eigen::Map<Eigen::VectorXd> C2_map(C2.col(i).data(), size_GJ);
            /* Set up matrix and vectors here */
            viennacl::copy(C2_map, vcl_rhs);
            if(opt.DecompositionMethod == "ViennaCL-ConjugateGradient")
                vcl_result = viennacl::linalg::solve(vcl_matrix, vcl_rhs, viennacl::linalg::cg_tag());
            else if(opt.DecompositionMethod == "ViennaCL-BiCGStab")
                vcl_result = viennacl::linalg::solve(vcl_matrix, vcl_rhs, viennacl::linalg::bicgstab_tag());
            else if(opt.DecompositionMethod == "ViennaCL-GMRES")
                vcl_result = viennacl::linalg::solve(vcl_matrix, vcl_rhs, viennacl::linalg::gmres_tag());
            viennacl::copy(vcl_result, C2_map);
            C1_.col(i) = C2_map;
        }

        clock_t time2 = clock();
        stat.timeDecomposition = (time2 - time1) / double(CLOCKS_PER_SEC);
    }
    else if(opt.DecompositionMethod == "ViennaCL-ILUT")
    {
        clock_t time1 = clock();
        Eigen::SparseMatrix<double> C1(size_GJ, size_GJ);
        data_func(C1, C2, data);

        viennacl::compressed_matrix<double> vcl_matrix(size_GJ, size_GJ);
        viennacl::copy(C1, vcl_matrix);
        viennacl::linalg::ilut_precond<viennacl::compressed_matrix<double>> vcl_ilut(vcl_matrix, viennacl::linalg::ilut_tag());
        viennacl::vector<double> vcl_rhs(size_GJ);
        viennacl::vector<double> vcl_result;
        for(int i = 0; i < size_AM; ++i)
        {
            Eigen::Map<Eigen::VectorXd> C2_map(C2.col(i).data(), size_GJ);
            /* Set up matrix and vectors here */
            viennacl::copy(C2_map, vcl_rhs);
            //solve (e.g. using conjugate gradient solver)
            vcl_result = viennacl::linalg::solve(vcl_matrix, vcl_rhs, viennacl::linalg::bicgstab_tag(), vcl_ilut); //preconditioner here
            viennacl::copy(vcl_result, C2_map);
            C1_.col(i) = C2_map;
        }

        clock_t time2 = clock();
        stat.timeDecomposition = (time2 - time1) / double(CLOCKS_PER_SEC);
    }
    else if(opt.DecompositionMethod == "ViennaCL-Jacobi")
    {
        clock_t time1 = clock();
        Eigen::SparseMatrix<double> C1(size_GJ, size_GJ);
        data_func(C1, C2, data);

        viennacl::compressed_matrix<double> vcl_matrix(size_GJ, size_GJ);
        viennacl::copy(C1, vcl_matrix);
        viennacl::linalg::jacobi_precond<viennacl::compressed_matrix<double>> vcl_jacobi(vcl_matrix, viennacl::linalg::jacobi_tag());
        viennacl::vector<double> vcl_rhs(size_GJ);
        viennacl::vector<double> vcl_result;
        for(int i = 0; i < size_AM; ++i)
        {
            Eigen::Map<Eigen::VectorXd> C2_map(C2.col(i).data(), size_GJ);
            /* Set up matrix and vectors here */
            viennacl::copy(C2_map, vcl_rhs);
            //solve (e.g. using conjugate gradient solver)
            vcl_result = viennacl::linalg::solve(vcl_matrix, vcl_rhs, viennacl::linalg::cg_tag(), vcl_jacobi); //preconditioner here
            viennacl::copy(vcl_result, C2_map);
            C1_.col(i) = C2_map;
        }

        clock_t time2 = clock();
        stat.timeDecomposition = (time2 - time1) / double(CLOCKS_PER_SEC);
    }
#endif

    return stat;
}

void readPnPdata(const std::string& filename,
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


void readpTopdata(const std::string& filename,
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



void readHandEyedata(const std::string& filename,
                  Eigen::Matrix3d& R0,
                  Eigen::Vector3d& t0,
                  std::vector<Eigen::Matrix4d>& As,
                  std::vector<Eigen::Matrix4d>& Bs)
{
    std::ifstream input(filename);
    input >> R0(0, 0) >> R0(0, 1) >> R0(0, 2) >>
          R0(1, 0) >> R0(1, 1) >> R0(1, 2) >>
          R0(2, 0) >> R0(2, 1) >> R0(2, 2);
    input >> t0(0) >> t0(1) >> t0(2);
    int num = 0;
    input >> num;
    As.resize(num);
    Bs.resize(num);

    for(int i = 0; i < num; ++i)
    {
        input >> As[i](0, 0) >> As[i](0, 1) >> As[i](0, 2) >> As[i](0, 3)
              >> As[i](1, 0) >> As[i](1, 1) >> As[i](1, 2) >> As[i](1, 3)
              >> As[i](2, 0) >> As[i](2, 1) >> As[i](2, 2) >> As[i](2, 3)
              >> As[i](3, 0) >> As[i](3, 1) >> As[i](3, 2) >> As[i](3, 3);
    }
    for(int i = 0; i < num; ++i)
    {
        input >> Bs[i](0, 0) >> Bs[i](0, 1) >> Bs[i](0, 2) >> Bs[i](0, 3)
              >> Bs[i](1, 0) >> Bs[i](1, 1) >> Bs[i](1, 2) >> Bs[i](1, 3)
              >> Bs[i](2, 0) >> Bs[i](2, 1) >> Bs[i](2, 2) >> Bs[i](2, 3)
              >> Bs[i](3, 0) >> Bs[i](3, 1) >> Bs[i](3, 2) >> Bs[i](3, 3);
    }
    input.close();
}
