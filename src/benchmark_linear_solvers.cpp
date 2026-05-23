// LibQPEP linear solver backend benchmark.

#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "QPEP_grobner.h"
#include "hand_eye_WQD.h"
#include "misc_hand_eye_funcs.h"
#include "misc_pTop_funcs.h"
#include "misc_pnp_funcs.h"
#include "pTop_WQD.h"
#include "pnp_WQD.h"
#include "solver_WQ_1_2_3_4_5_9_13_17_33_49_approx.h"
#include "solver_WQ_approx.h"
#include "utils.h"

struct BenchmarkCase
{
    std::string problem;
    std::string module;
    Eigen::MatrixXd W;
    Eigen::MatrixXd Q;
    Eigen::MatrixXd coef_J_pure;
    Eigen::MatrixXd coefs_tq;
    Eigen::MatrixXd pinvG;
    solver_func_handle solver_func = nullptr;
    mon_J_pure_func_handle mon_J_pure_func = nullptr;
    t_func_handle t_func = nullptr;
};

struct BenchmarkStats
{
    std::string backend;
    int loops = 0;
    int successes = 0;
    int failures = 0;
    double wall_time = 0.0;
    double data_prepare_time = 0.0;
    double decomposition_time = 0.0;
    double grobner_time = 0.0;
    double eigen_time = 0.0;
    double loss = 0.0;
};

static void print_usage()
{
    std::cout << "Usage: ./LibQPEP-benchmark-linear-solvers [--problem PnP|PtoP|Hand-eye] [--data file] [--loops N] [--backends b1,b2,...]" << std::endl;
    std::cout << "Default problem: PtoP, default loops: 5, default backends: all available." << std::endl;
    std::cout << "Available linear solver backends: " << LinearSolverBackendsDescription() << std::endl;
}

static std::string format_double(const double value)
{
    std::ostringstream out;
    out << std::setprecision(6) << value;
    return out.str();
}

static std::string format_average(const double value, const bool has_value)
{
    return has_value ? format_double(value) : "n/a";
}

static void print_table(const std::string& title,
                        const std::vector<std::string>& headers,
                        const std::vector<std::vector<std::string>>& rows)
{
    std::vector<size_t> widths(headers.size(), 0);
    for(size_t i = 0; i < headers.size(); ++i) {
        widths[i] = headers[i].size();
    }
    for(const std::vector<std::string>& row : rows) {
        for(size_t i = 0; i < row.size() && i < widths.size(); ++i) {
            widths[i] = std::max(widths[i], row[i].size());
        }
    }

    const auto print_border = [&widths]() {
        std::cout << "+";
        for(size_t width : widths) {
            std::cout << std::string(width + 2, '-') << "+";
        }
        std::cout << std::endl;
    };

    const auto print_row = [&widths](const std::vector<std::string>& row) {
        std::cout << "|";
        for(size_t i = 0; i < widths.size(); ++i) {
            const std::string cell = i < row.size() ? row[i] : "";
            std::cout << " " << std::left << std::setw(widths[i]) << cell << " |";
        }
        std::cout << std::right << std::endl;
    };

    std::cout << std::endl << title << std::endl;
    print_border();
    print_row(headers);
    print_border();
    for(const std::vector<std::string>& row : rows) {
        print_row(row);
    }
    print_border();
}

static std::string resolve_data_file(const std::string& data_file)
{
    const std::string src_dir(CURRENT_SRC_DIR);
    if(data_file.empty()) {
        return data_file;
    }
    if(data_file[0] == '/') {
        if(data_file.compare(0, 6, "/data/") == 0) {
            return src_dir + data_file;
        }
        return data_file;
    }
    return src_dir + "/" + data_file;
}

static void make_reduced_WQ(Eigen::MatrixXd& W_reduced,
                            Eigen::MatrixXd& Q_reduced,
                            const Eigen::MatrixXd& W,
                            const Eigen::MatrixXd& Q)
{
    W_reduced = W.topRows(3);
    Q_reduced = Q.topRows(3);

    W_reduced.row(0) = W.row(0) + W.row(1) + W.row(2);
    W_reduced.row(1) = W.row(1) + W.row(2) + W.row(3);
    W_reduced.row(2) = W.row(2) + W.row(3) + W.row(0);

    Q_reduced.row(0) = Q.row(0) + Q.row(1) + Q.row(2);
    Q_reduced.row(1) = Q.row(1) + Q.row(2) + Q.row(3);
    Q_reduced.row(2) = Q.row(2) + Q.row(3) + Q.row(0);
}

static BenchmarkCase prepare_pnp_case(const std::string& data_file)
{
    Eigen::Matrix3d R0;
    Eigen::Vector3d t0;
    Eigen::Matrix3d K;
    std::vector<Eigen::Vector3d> world_pt;
    std::vector<Eigen::Vector2d> image_pt;
    K.setZero();
    readPnPdata(data_file, R0, t0, K, world_pt, image_pt);

    Eigen::Matrix<double, 4, 64> W;
    Eigen::Matrix<double, 4, 4> Q;
    Eigen::Matrix<double, 3, 37> D;
    Eigen::Matrix<double, 4, 24> coef_f_q_sym;
    Eigen::Matrix<double, 1, 70> coef_J_pure;
    Eigen::Matrix<double, 3, 10> coefs_tq;
    Eigen::Matrix<double, 3, 3> pinvG;
    pnp_WQD(W, Q, D, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG,
            image_pt, world_pt, K, 1e-8);

    BenchmarkCase test_case;
    test_case.problem = "PnP";
    test_case.module = "solver_WQ_1_2_3_4_5_9_13_17_33_49_approx";
    make_reduced_WQ(test_case.W, test_case.Q, W, Q);
    test_case.coef_J_pure = coef_J_pure;
    test_case.coefs_tq = coefs_tq;
    test_case.pinvG = pinvG;
    test_case.solver_func = reinterpret_cast<solver_func_handle>(solver_WQ_1_2_3_4_5_9_13_17_33_49_approx);
    test_case.mon_J_pure_func = reinterpret_cast<mon_J_pure_func_handle>(mon_J_pure_pnp_func);
    test_case.t_func = reinterpret_cast<t_func_handle>(t_pnp_func);
    return test_case;
}

static BenchmarkCase prepare_ptop_case(const std::string& data_file)
{
    Eigen::Matrix3d R0;
    Eigen::Vector3d t0;
    std::vector<Eigen::Vector3d> rr;
    std::vector<Eigen::Vector3d> bb;
    std::vector<Eigen::Vector3d> nv;
    readpTopdata(data_file, R0, t0, rr, bb, nv);

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

    BenchmarkCase test_case;
    test_case.problem = "PtoP";
    test_case.module = "solver_WQ_approx";
    make_reduced_WQ(test_case.W, test_case.Q, W, Q);
    test_case.coef_J_pure = coef_J_pure;
    test_case.coefs_tq = coefs_tq;
    test_case.pinvG = pinvG;
    test_case.solver_func = reinterpret_cast<solver_func_handle>(solver_WQ_approx);
    test_case.mon_J_pure_func = reinterpret_cast<mon_J_pure_func_handle>(mon_J_pure_pTop_func);
    test_case.t_func = reinterpret_cast<t_func_handle>(t_pTop_func);
    return test_case;
}

static BenchmarkCase prepare_hand_eye_case(const std::string& data_file)
{
    Eigen::Matrix3d R0;
    Eigen::Vector3d t0;
    std::vector<Eigen::Matrix4d> As;
    std::vector<Eigen::Matrix4d> Bs;
    readHandEyedata(data_file, R0, t0, As, Bs);

    Eigen::Matrix<double, 4, 64> W;
    Eigen::Matrix<double, 4, 4> Q;
    Eigen::Matrix<double, 3, 28> D;
    Eigen::Matrix<double, 3, 9> G;
    Eigen::Vector3d c;
    Eigen::Matrix<double, 4, 24> coef_f_q_sym;
    Eigen::Matrix<double, 1, 85> coef_J_pure;
    Eigen::Matrix<double, 3, 11> coefs_tq;
    Eigen::Matrix<double, 3, 3> pinvG;
    hand_eye_WQD(W, Q, D, G, c, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG, As, Bs);

    BenchmarkCase test_case;
    test_case.problem = "Hand-eye";
    test_case.module = "solver_WQ_approx";
    make_reduced_WQ(test_case.W, test_case.Q, W, Q);
    test_case.coef_J_pure = coef_J_pure;
    test_case.coefs_tq = coefs_tq;
    test_case.pinvG = pinvG;
    test_case.solver_func = reinterpret_cast<solver_func_handle>(solver_WQ_approx);
    test_case.mon_J_pure_func = reinterpret_cast<mon_J_pure_func_handle>(mon_J_pure_hand_eye_func);
    test_case.t_func = reinterpret_cast<t_func_handle>(t_hand_eye_func);
    return test_case;
}

static std::vector<std::string> parse_backend_list(const std::string& backends)
{
    if(backends.empty() || backends == "all") {
        return AvailableLinearSolverBackends();
    }

    std::vector<std::string> parsed;
    std::stringstream stream(backends);
    std::string item;
    while(std::getline(stream, item, ','))
    {
        std::string normalized = NormalizeLinearSolverBackend(item);
        if(normalized.empty())
        {
            std::cout << "Unknown linear solver backend: " << item << std::endl;
            std::cout << "Available linear solver backends: " << LinearSolverBackendsDescription() << std::endl;
            std::exit(1);
        }
        parsed.push_back(normalized);
    }
    return parsed;
}

static BenchmarkStats benchmark_backend(const BenchmarkCase& test_case,
                                        const std::string& backend,
                                        const int loops)
{
    BenchmarkStats result;
    result.backend = backend;
    result.loops = loops;

    for(int i = 0; i < loops; ++i)
    {
        Eigen::Matrix3d R;
        Eigen::Vector3d t;
        Eigen::Matrix4d X;
        double min_values[40];
        QPEP_options opt;
        opt.ModuleName = test_case.module;
        opt.DecompositionMethod = backend;

        const clock_t start = clock();
        QPEP_runtime stat = QPEP_WQ_grobner(R, t, X, min_values,
                                            test_case.W, test_case.Q,
                                            test_case.solver_func,
                                            test_case.mon_J_pure_func,
                                            test_case.t_func,
                                            test_case.coef_J_pure,
                                            test_case.coefs_tq,
                                            test_case.pinvG,
                                            nullptr, opt);
        const clock_t end = clock();
        result.wall_time += (end - start) / double(CLOCKS_PER_SEC);

        if(stat.statusDecomposition != 0 || stat.statusGrobner != 0 || stat.statusEigen != 0)
        {
            ++result.failures;
            continue;
        }

        ++result.successes;
        result.data_prepare_time += stat.timeDecompositionDataPrepare;
        result.decomposition_time += stat.timeDecomposition;
        result.grobner_time += stat.timeGrobner;
        result.eigen_time += stat.timeEigen;
        result.loss += stat.lossGrobner;
    }

    result.wall_time /= loops;
    if(result.successes == 0)
    {
        return result;
    }

    result.data_prepare_time /= result.successes;
    result.decomposition_time /= result.successes;
    result.grobner_time /= result.successes;
    result.eigen_time /= result.successes;
    result.loss /= result.successes;
    return result;
}

int main(int argc, char** argv)
{
    std::string problem = "PtoP";
    std::string data_file = "/data/pTop_data-4096pt-1.txt";
    std::string backends = "all";
    int loops = 5;

    for(int i = 1; i < argc; ++i)
    {
        std::string arg(argv[i]);
        if(arg == "-h" || arg == "--help")
        {
            print_usage();
            return 0;
        }
        else if(arg == "--list-linear-solvers" || arg == "--list-backends")
        {
            std::cout << LinearSolverBackendsDescription() << std::endl;
            return 0;
        }
        else if(arg == "--problem" && i + 1 < argc)
        {
            problem = argv[++i];
        }
        else if(arg.compare(0, 10, "--problem=") == 0)
        {
            problem = arg.substr(10);
        }
        else if(arg == "--data" && i + 1 < argc)
        {
            data_file = argv[++i];
        }
        else if(arg.compare(0, 7, "--data=") == 0)
        {
            data_file = arg.substr(7);
        }
        else if(arg == "--loops" && i + 1 < argc)
        {
            loops = std::atoi(argv[++i]);
        }
        else if(arg.compare(0, 8, "--loops=") == 0)
        {
            loops = std::atoi(arg.substr(8).c_str());
        }
        else if(arg == "--backends" && i + 1 < argc)
        {
            backends = argv[++i];
        }
        else if(arg.compare(0, 11, "--backends=") == 0)
        {
            backends = arg.substr(11);
        }
        else
        {
            std::cout << "Unknown argument: " << arg << std::endl;
            print_usage();
            return 1;
        }
    }

    if(loops <= 0)
    {
        std::cout << "--loops must be positive." << std::endl;
        return 1;
    }

    if(problem == "PnP" && data_file == "/data/pTop_data-4096pt-1.txt") {
        data_file = "/data/pnp_data-500pt-1.txt";
    } else if(problem == "Hand-eye" && data_file == "/data/pTop_data-4096pt-1.txt") {
        data_file = "/data/hand_eye_data-20pt-1.txt";
    }

    const std::string full_file = resolve_data_file(data_file);
    BenchmarkCase test_case;
    if(problem == "PnP")
    {
        test_case = prepare_pnp_case(full_file);
    }
    else if(problem == "PtoP")
    {
        test_case = prepare_ptop_case(full_file);
    }
    else if(problem == "Hand-eye")
    {
        test_case = prepare_hand_eye_case(full_file);
    }
    else
    {
        std::cout << "Unknown problem: " << problem << std::endl;
        print_usage();
        return 1;
    }

    print_table("Benchmark Case",
                {"Problem", "Data File", "Loops"},
                {{test_case.problem, full_file, std::to_string(loops)}});

    const std::vector<std::string> backend_list = parse_backend_list(backends);
    std::vector<BenchmarkStats> results;
    for(const std::string& backend : backend_list)
    {
        results.push_back(benchmark_backend(test_case, backend, loops));
    }

    std::vector<std::vector<std::string>> status_rows;
    std::vector<std::vector<std::string>> timing_rows;
    for(const BenchmarkStats& result : results)
    {
        const bool has_success = result.successes > 0;
        status_rows.push_back({
                result.backend,
                std::to_string(result.loops),
                std::to_string(result.successes),
                std::to_string(result.failures),
                format_average(result.loss, has_success)
        });
        timing_rows.push_back({
                result.backend,
                format_double(result.wall_time),
                format_average(result.data_prepare_time, has_success),
                format_average(result.decomposition_time, has_success),
                format_average(result.grobner_time, has_success),
                format_average(result.eigen_time, has_success)
        });
    }

    print_table("Solve Status",
                {"Backend", "Loops", "Successes", "Failures", "Avg Loss"},
                status_rows);
    print_table("Average Timings (s)",
                {"Backend", "Wall", "Data Prep", "Decomposition", "Grobner", "Eigen"},
                timing_rows);

    return 0;
}
