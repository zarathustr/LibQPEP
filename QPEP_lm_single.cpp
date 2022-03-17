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
// QPEP_lm_single.cpp: Refines QPEPs by Levenberg-Marquardt algorithm (LMA)



#include "QPEP_lm_single.h"
#include "utils.h"


#include "linearLeastSquares.h"
#include "rt_nonfinite.h"
#include "xnrm2.h"
#include "xzlarf.h"
#include "xzlarfg.h"
#include <cmath>

typedef bool boolean_T;
double rtInf = std::numeric_limits<double>::max();
double rtNaN = NAN;

static int div_nde_s32_floor(int numerator, int denominator)
{
    int b_numerator;
    if (((numerator < 0) != (denominator < 0)) &&
        (numerator % denominator != 0)) {
        b_numerator = -1;
    } else {
        b_numerator = 0;
    }
    return numerator / denominator + b_numerator;
}



            void xzlarf(int m, int n, int iv0, double tau, double C[45], int ic0,
                        double work[5])
            {
                int i;
                int ia;
                int lastc;
                int lastv;
                if (tau != 0.0) {
                    boolean_T exitg2;
                    lastv = m;
                    i = iv0 + m;
                    while ((lastv > 0) && (C[i - 2] == 0.0)) {
                        lastv--;
                        i--;
                    }
                    lastc = n - 1;
                    exitg2 = false;
                    while ((!exitg2) && (lastc + 1 > 0)) {
                        int exitg1;
                        i = ic0 + lastc * 9;
                        ia = i;
                        do {
                            exitg1 = 0;
                            if (ia <= (i + lastv) - 1) {
                                if (C[ia - 1] != 0.0) {
                                    exitg1 = 1;
                                } else {
                                    ia++;
                                }
                            } else {
                                lastc--;
                                exitg1 = 2;
                            }
                        } while (exitg1 == 0);
                        if (exitg1 == 1) {
                            exitg2 = true;
                        }
                    }
                } else {
                    lastv = 0;
                    lastc = -1;
                }
                if (lastv > 0) {
                    double c;
                    int b_i;
                    int iac;
                    if (lastc + 1 != 0) {
                        if (0 <= lastc) {
                            std::memset(&work[0], 0, (lastc + 1) * sizeof(double));
                        }
                        b_i = ic0 + 9 * lastc;
                        for (iac = ic0; iac <= b_i; iac += 9) {
                            c = 0.0;
                            i = (iac + lastv) - 1;
                            for (ia = iac; ia <= i; ia++) {
                                c += C[ia - 1] * C[((iv0 + ia) - iac) - 1];
                            }
                            i = div_nde_s32_floor(iac - ic0, 9);
                            work[i] += c;
                        }
                    }
                    if (!(-tau == 0.0)) {
                        i = ic0;
                        for (iac = 0; iac <= lastc; iac++) {
                            if (work[iac] != 0.0) {
                                c = work[iac] * -tau;
                                b_i = lastv + i;
                                for (ia = i; ia < b_i; ia++) {
                                    C[ia - 1] += C[((iv0 + ia) - i) - 1] * c;
                                }
                            }
                            i += 9;
                        }
                    }
                }
            }


            double xnrm2(int n, const double x[45], int ix0)
            {
                double scale;
                double y;
                int kend;
                y = 0.0;
                scale = 3.3121686421112381E-170;
                kend = (ix0 + n) - 1;
                for (int k{ix0}; k <= kend; k++) {
                    double absxk;
                    absxk = std::abs(x[k - 1]);
                    if (absxk > scale) {
                        double t;
                        t = scale / absxk;
                        y = y * t * t + 1.0;
                        scale = absxk;
                    } else {
                        double t;
                        t = absxk / scale;
                        y += t * t;
                    }
                }
                return scale * std::sqrt(y);
            }
static double rt_hypotd_snf(double u0, double u1)
{
    double a;
    double y;
    a = std::abs(u0);
    y = std::abs(u1);
    if (a < y) {
        a /= y;
        y *= std::sqrt(a * a + 1.0);
    } else if (a > y) {
        y /= a;
        y = a * std::sqrt(y * y + 1.0);
    } else if (!std::isnan(y)) {
        y = a * 1.4142135623730951;
    }
    return y;
}

            double xzlarfg(int n, double *alpha1, double x[45], int ix0)
            {
                double tau;
                double xnorm;
                tau = 0.0;
                xnorm = xnrm2(n - 1, x, ix0);
                if (xnorm != 0.0) {
                    double beta1;
                    beta1 = rt_hypotd_snf(*alpha1, xnorm);
                    if (*alpha1 >= 0.0) {
                        beta1 = -beta1;
                    }
                    if (std::abs(beta1) < 1.0020841800044864E-292) {
                        int i;
                        int k;
                        int knt;
                        knt = -1;
                        i = (ix0 + n) - 2;
                        do {
                            knt++;
                            for (k = ix0; k <= i; k++) {
                                x[k - 1] *= 9.9792015476736E+291;
                            }
                            beta1 *= 9.9792015476736E+291;
                            *alpha1 *= 9.9792015476736E+291;
                        } while (!(std::abs(beta1) >= 1.0020841800044864E-292));
                        beta1 = rt_hypotd_snf(*alpha1, xnrm2(n - 1, x, ix0));
                        if (*alpha1 >= 0.0) {
                            beta1 = -beta1;
                        }
                        tau = (beta1 - *alpha1) / beta1;
                        xnorm = 1.0 / (*alpha1 - beta1);
                        for (k = ix0; k <= i; k++) {
                            x[k - 1] *= xnorm;
                        }
                        for (k = 0; k <= knt; k++) {
                            beta1 *= 1.0020841800044864E-292;
                        }
                        *alpha1 = beta1;
                    } else {
                        int i;
                        tau = (beta1 - *alpha1) / beta1;
                        xnorm = 1.0 / (*alpha1 - beta1);
                        i = (ix0 + n) - 2;
                        for (int k{ix0}; k <= i; k++) {
                            x[k - 1] *= xnorm;
                        }
                        *alpha1 = beta1;
                    }
                }
                return tau;
            }


namespace coder {
    namespace optim {
        namespace coder {
            namespace levenbergMarquardt {
                void linearLeastSquares(double lhs[45], double rhs[9], double dx[5])
                {
                    double jpvt[5];
                    double tau[5];
                    double vn1[5];
                    double vn2[5];
                    double work[5];
                    double d;
                    double temp;
                    int b_i;
                    int i;
                    int ii;
                    int ix;
                    int iy;
                    int j;
                    int k;
                    int nfxd;
                    int temp_tmp;
                    for (i = 0; i < 5; i++) {
                        jpvt[i] = 0.0;
                    }
                    nfxd = -1;
                    for (j = 0; j < 5; j++) {
                        if (jpvt[j] != 0.0) {
                            nfxd++;
                            if (j + 1 != nfxd + 1) {
                                ix = j * 9;
                                iy = nfxd * 9;
                                for (k = 0; k < 9; k++) {
                                    temp_tmp = ix + k;
                                    temp = lhs[temp_tmp];
                                    b_i = iy + k;
                                    lhs[temp_tmp] = lhs[b_i];
                                    lhs[b_i] = temp;
                                }
                                jpvt[j] = jpvt[nfxd];
                                jpvt[nfxd] = static_cast<double>(j) + 1.0;
                            } else {
                                jpvt[j] = static_cast<double>(j) + 1.0;
                            }
                        } else {
                            jpvt[j] = static_cast<double>(j) + 1.0;
                        }
                        tau[j] = 0.0;
                        work[j] = 0.0;
                    }
                    for (i = 0; i <= nfxd; i++) {
                        ii = i * 9 + i;
                        temp = lhs[ii];
                        d = xzlarfg(9 - i, &temp, lhs, ii + 2);
                        tau[i] = d;
                        lhs[ii] = temp;
                        if (i + 1 < 5) {
                            temp = lhs[ii];
                            lhs[ii] = 1.0;
                            xzlarf(9 - i, 4 - i, ii + 1, d, lhs, ii + 10, work);
                            lhs[ii] = temp;
                        }
                    }
                    if (nfxd + 1 < 5) {
                        for (i = 0; i < 5; i++) {
                            work[i] = 0.0;
                            vn1[i] = 0.0;
                            vn2[i] = 0.0;
                        }
                        b_i = nfxd + 2;
                        for (j = b_i; j < 6; j++) {
                            d = xnrm2(8 - nfxd, lhs, (nfxd + (j - 1) * 9) + 2);
                            vn1[j - 1] = d;
                            vn2[j - 1] = d;
                        }
                        for (i = b_i; i < 6; i++) {
                            double s;
                            int ip1;
                            ip1 = i + 1;
                            j = (i - 1) * 9;
                            ii = (j + i) - 1;
                            iy = 6 - i;
                            nfxd = -1;
                            if (6 - i > 1) {
                                temp = std::abs(vn1[i - 1]);
                                for (k = 2; k <= iy; k++) {
                                    s = std::abs(vn1[(i + k) - 2]);
                                    if (s > temp) {
                                        nfxd = k - 2;
                                        temp = s;
                                    }
                                }
                            }
                            iy = i + nfxd;
                            if (iy + 1 != i) {
                                ix = iy * 9;
                                for (k = 0; k < 9; k++) {
                                    temp_tmp = ix + k;
                                    temp = lhs[temp_tmp];
                                    nfxd = j + k;
                                    lhs[temp_tmp] = lhs[nfxd];
                                    lhs[nfxd] = temp;
                                }
                                temp = jpvt[iy];
                                jpvt[iy] = jpvt[i - 1];
                                jpvt[i - 1] = temp;
                                vn1[iy] = vn1[i - 1];
                                vn2[iy] = vn2[i - 1];
                            }
                            temp = lhs[ii];
                            d = xzlarfg(10 - i, &temp, lhs, ii + 2);
                            tau[i - 1] = d;
                            lhs[ii] = temp;
                            if (i < 5) {
                                temp = lhs[ii];
                                lhs[ii] = 1.0;
                                xzlarf(10 - i, 5 - i, ii + 1, d, lhs, ii + 10,
                                                            work);
                                lhs[ii] = temp;
                            }
                            for (j = ip1; j < 6; j++) {
                                iy = i + (j - 1) * 9;
                                d = vn1[j - 1];
                                if (d != 0.0) {
                                    temp = std::abs(lhs[iy - 1]) / d;
                                    temp = 1.0 - temp * temp;
                                    if (temp < 0.0) {
                                        temp = 0.0;
                                    }
                                    s = d / vn2[j - 1];
                                    s = temp * (s * s);
                                    if (s <= 1.4901161193847656E-8) {
                                        d = xnrm2(9 - i, lhs, iy + 1);
                                        vn1[j - 1] = d;
                                        vn2[j - 1] = d;
                                    } else {
                                        vn1[j - 1] = d * std::sqrt(temp);
                                    }
                                }
                            }
                        }
                    }
                    for (j = 0; j < 5; j++) {
                        if (tau[j] != 0.0) {
                            temp = rhs[j];
                            b_i = j + 2;
                            for (i = b_i; i < 10; i++) {
                                temp += lhs[(i + 9 * j) - 1] * rhs[i - 1];
                            }
                            temp *= tau[j];
                            if (temp != 0.0) {
                                rhs[j] -= temp;
                                for (i = b_i; i < 10; i++) {
                                    rhs[i - 1] -= lhs[(i + 9 * j) - 1] * temp;
                                }
                            }
                        }
                    }
                    for (i = 0; i < 5; i++) {
                        dx[i] = rhs[i];
                    }
                    for (j = 4; j >= 0; j--) {
                        iy = j + j * 9;
                        dx[j] /= lhs[iy];
                        for (i = 0; i < j; i++) {
                            ix = (j - i) - 1;
                            dx[ix] -= dx[j] * lhs[(iy - i) - 1];
                        }
                    }
                    for (b_i = 0; b_i < 5; b_i++) {
                        tau[b_i] = dx[b_i];
                    }
                    for (b_i = 0; b_i < 5; b_i++) {
                        dx[static_cast<int>(jpvt[b_i]) - 1] = tau[b_i];
                    }
                }

            } // namespace levenbergMarquardt
        } // namespace coder
    } // namespace optim
} // namespace coder


// Function Declarations
static void anon(const double W[256], const double Q[16], const double x[5],
                 double varargout_1[4], double varargout_2[20]);

// Function Definitions
//
// Arguments    : const double W[256]
//                const double Q[16]
//                const double x[5]
//                double varargout_1[4]
//                double varargout_2[20]
// Return Type  : void
//
static void anon(const double W[256], const double Q[16], const double x[5],
                 double varargout_1[4], double varargout_2[20])
{
    static const signed char b[16]{1, 0, 0, 0, 0, 1, 0, 0,
                                   0, 0, 1, 0, 0, 0, 0, 1};
    double b_b[64];
    double K[16];
    double varargout_2_tmp;
    int i1;
    int kidx;
    kidx = -1;
    for (i1 = 0; i1 < 4; i1++) {
        K[kidx + 1] = x[i1] * x[0];
        K[kidx + 2] = x[i1] * x[1];
        K[kidx + 3] = x[i1] * x[2];
        K[kidx + 4] = x[i1] * x[3];
        kidx += 4;
    }
    kidx = -1;
    for (i1 = 0; i1 < 4; i1++) {
        for (int i2{0}; i2 < 16; i2++) {
            b_b[(kidx + i2) + 1] = x[i1] * K[i2];
        }
        kidx += 16;
    }
    for (kidx = 0; kidx < 16; kidx++) {
        K[kidx] = Q[kidx] - x[4] * static_cast<double>(b[kidx]);
    }
    for (kidx = 0; kidx < 4; kidx++) {
        varargout_2_tmp = 0.0;
        for (i1 = 0; i1 < 64; i1++) {
            varargout_2_tmp += W[kidx + (i1 << 2)] * b_b[i1];
        }
        varargout_1[kidx] =
                varargout_2_tmp -
                (((K[kidx] * x[0] + K[kidx + 4] * x[1]) + K[kidx + 8] * x[2]) +
                 K[kidx + 12] * x[3]);
    }
    double ab_varargout_2_tmp;
    double ac_varargout_2_tmp;
    double ad_varargout_2_tmp;
    double b_varargout_2_tmp;
    double bb_varargout_2_tmp;
    double bc_varargout_2_tmp;
    double bd_varargout_2_tmp;
    double c_varargout_2_tmp;
    double cb_varargout_2_tmp;
    double cc_varargout_2_tmp;
    double cd_varargout_2_tmp;
    double d_varargout_2_tmp;
    double db_varargout_2_tmp;
    double dc_varargout_2_tmp;
    double dd_varargout_2_tmp;
    double e_varargout_2_tmp;
    double eb_varargout_2_tmp;
    double ec_varargout_2_tmp;
    double ed_varargout_2_tmp;
    double f_varargout_2_tmp;
    double fb_varargout_2_tmp;
    double fc_varargout_2_tmp;
    double fd_varargout_2_tmp;
    double g_varargout_2_tmp;
    double gb_varargout_2_tmp;
    double gc_varargout_2_tmp;
    double gd_varargout_2_tmp;
    double h_varargout_2_tmp;
    double hb_varargout_2_tmp;
    double hc_varargout_2_tmp;
    double i_varargout_2_tmp;
    double ib_varargout_2_tmp;
    double ic_varargout_2_tmp;
    double j_varargout_2_tmp;
    double jb_varargout_2_tmp;
    double jc_varargout_2_tmp;
    double k_varargout_2_tmp;
    double kb_varargout_2_tmp;
    double kc_varargout_2_tmp;
    double l_varargout_2_tmp;
    double lb_varargout_2_tmp;
    double lc_varargout_2_tmp;
    double m_varargout_2_tmp;
    double mb_varargout_2_tmp;
    double mc_varargout_2_tmp;
    double n_varargout_2_tmp;
    double nb_varargout_2_tmp;
    double nc_varargout_2_tmp;
    double o_varargout_2_tmp;
    double ob_varargout_2_tmp;
    double oc_varargout_2_tmp;
    double p_varargout_2_tmp;
    double pb_varargout_2_tmp;
    double pc_varargout_2_tmp;
    double q_varargout_2_tmp;
    double qb_varargout_2_tmp;
    double qc_varargout_2_tmp;
    double r_varargout_2_tmp;
    double rb_varargout_2_tmp;
    double rc_varargout_2_tmp;
    double s_varargout_2_tmp;
    double sb_varargout_2_tmp;
    double sc_varargout_2_tmp;
    double t_varargout_2_tmp;
    double tb_varargout_2_tmp;
    double tc_varargout_2_tmp;
    double u_varargout_2_tmp;
    double ub_varargout_2_tmp;
    double uc_varargout_2_tmp;
    double v_varargout_2_tmp;
    double vb_varargout_2_tmp;
    double vc_varargout_2_tmp;
    double w_varargout_2_tmp;
    double wb_varargout_2_tmp;
    double wc_varargout_2_tmp;
    double x_varargout_2_tmp;
    double xb_varargout_2_tmp;
    double xc_varargout_2_tmp;
    double y_varargout_2_tmp;
    double yb_varargout_2_tmp;
    double yc_varargout_2_tmp;
    varargout_2_tmp = x[1] * x[1];
    b_varargout_2_tmp = x[2] * x[2];
    c_varargout_2_tmp = x[3] * x[3];
    d_varargout_2_tmp = x[0] * x[0];
    varargout_2[0] =
            (((((((((((((((((((((((((((((((((((((-Q[0] + x[4]) +
                                                W[0] * d_varargout_2_tmp * 3.0) +
                                               W[20] * varargout_2_tmp) +
                                              W[40] * b_varargout_2_tmp) +
                                             W[60] * c_varargout_2_tmp) +
                                            W[68] * varargout_2_tmp) +
                                           W[80] * varargout_2_tmp) +
                                          W[136] * b_varargout_2_tmp) +
                                         W[160] * b_varargout_2_tmp) +
                                        W[204] * c_varargout_2_tmp) +
                                       W[240] * c_varargout_2_tmp) +
                                      x[0] * W[4] * x[1] * 2.0) +
                                     x[0] * W[8] * x[2] * 2.0) +
                                    x[0] * W[16] * x[1] * 2.0) +
                                   x[0] * W[12] * x[3] * 2.0) +
                                  x[1] * W[24] * x[2]) +
                                 x[0] * W[32] * x[2] * 2.0) +
                                x[1] * W[28] * x[3]) +
                               x[1] * W[36] * x[2]) +
                              x[0] * W[48] * x[3] * 2.0) +
                             x[2] * W[44] * x[3]) +
                            x[1] * W[52] * x[3]) +
                           x[0] * W[64] * x[1] * 2.0) +
                          x[2] * W[56] * x[3]) +
                         x[1] * W[72] * x[2]) +
                        x[1] * W[76] * x[3]) +
                       x[1] * W[96] * x[2]) +
                      x[1] * W[112] * x[3]) +
                     x[0] * W[128] * x[2] * 2.0) +
                    x[1] * W[132] * x[2]) +
                   x[1] * W[144] * x[2]) +
                  x[2] * W[140] * x[3]) +
                 x[2] * W[176] * x[3]) +
                x[0] * W[192] * x[3] * 2.0) +
               x[1] * W[196] * x[3]) +
              x[2] * W[200] * x[3]) +
             x[1] * W[208] * x[3]) +
            x[2] * W[224] * x[3];
    varargout_2[1] =
            ((((((((((((((((((((((((((((((((((((-Q[1] +
                                                W[1] * d_varargout_2_tmp * 3.0) +
                                               W[21] * varargout_2_tmp) +
                                              W[41] * b_varargout_2_tmp) +
                                             W[61] * c_varargout_2_tmp) +
                                            W[69] * varargout_2_tmp) +
                                           W[81] * varargout_2_tmp) +
                                          W[137] * b_varargout_2_tmp) +
                                         W[161] * b_varargout_2_tmp) +
                                        W[205] * c_varargout_2_tmp) +
                                       W[241] * c_varargout_2_tmp) +
                                      x[0] * W[5] * x[1] * 2.0) +
                                     x[0] * W[9] * x[2] * 2.0) +
                                    x[0] * W[17] * x[1] * 2.0) +
                                   x[0] * W[13] * x[3] * 2.0) +
                                  x[1] * W[25] * x[2]) +
                                 x[0] * W[33] * x[2] * 2.0) +
                                x[1] * W[29] * x[3]) +
                               x[1] * W[37] * x[2]) +
                              x[0] * W[49] * x[3] * 2.0) +
                             x[2] * W[45] * x[3]) +
                            x[1] * W[53] * x[3]) +
                           x[0] * W[65] * x[1] * 2.0) +
                          x[2] * W[57] * x[3]) +
                         x[1] * W[73] * x[2]) +
                        x[1] * W[77] * x[3]) +
                       x[1] * W[97] * x[2]) +
                      x[1] * W[113] * x[3]) +
                     x[0] * W[129] * x[2] * 2.0) +
                    x[1] * W[133] * x[2]) +
                   x[1] * W[145] * x[2]) +
                  x[2] * W[141] * x[3]) +
                 x[2] * W[177] * x[3]) +
                x[0] * W[193] * x[3] * 2.0) +
               x[1] * W[197] * x[3]) +
              x[2] * W[201] * x[3]) +
             x[1] * W[209] * x[3]) +
            x[2] * W[225] * x[3];
    varargout_2[2] =
            ((((((((((((((((((((((((((((((((((((-Q[2] +
                                                W[2] * d_varargout_2_tmp * 3.0) +
                                               W[22] * varargout_2_tmp) +
                                              W[42] * b_varargout_2_tmp) +
                                             W[62] * c_varargout_2_tmp) +
                                            W[70] * varargout_2_tmp) +
                                           W[82] * varargout_2_tmp) +
                                          W[138] * b_varargout_2_tmp) +
                                         W[162] * b_varargout_2_tmp) +
                                        W[206] * c_varargout_2_tmp) +
                                       W[242] * c_varargout_2_tmp) +
                                      x[0] * W[6] * x[1] * 2.0) +
                                     x[0] * W[10] * x[2] * 2.0) +
                                    x[0] * W[18] * x[1] * 2.0) +
                                   x[0] * W[14] * x[3] * 2.0) +
                                  x[1] * W[26] * x[2]) +
                                 x[0] * W[34] * x[2] * 2.0) +
                                x[1] * W[30] * x[3]) +
                               x[1] * W[38] * x[2]) +
                              x[0] * W[50] * x[3] * 2.0) +
                             x[2] * W[46] * x[3]) +
                            x[1] * W[54] * x[3]) +
                           x[0] * W[66] * x[1] * 2.0) +
                          x[2] * W[58] * x[3]) +
                         x[1] * W[74] * x[2]) +
                        x[1] * W[78] * x[3]) +
                       x[1] * W[98] * x[2]) +
                      x[1] * W[114] * x[3]) +
                     x[0] * W[130] * x[2] * 2.0) +
                    x[1] * W[134] * x[2]) +
                   x[1] * W[146] * x[2]) +
                  x[2] * W[142] * x[3]) +
                 x[2] * W[178] * x[3]) +
                x[0] * W[194] * x[3] * 2.0) +
               x[1] * W[198] * x[3]) +
              x[2] * W[202] * x[3]) +
             x[1] * W[210] * x[3]) +
            x[2] * W[226] * x[3];
    varargout_2[3] =
            ((((((((((((((((((((((((((((((((((((-Q[3] +
                                                W[3] * d_varargout_2_tmp * 3.0) +
                                               W[23] * varargout_2_tmp) +
                                              W[43] * b_varargout_2_tmp) +
                                             W[63] * c_varargout_2_tmp) +
                                            W[71] * varargout_2_tmp) +
                                           W[83] * varargout_2_tmp) +
                                          W[139] * b_varargout_2_tmp) +
                                         W[163] * b_varargout_2_tmp) +
                                        W[207] * c_varargout_2_tmp) +
                                       W[243] * c_varargout_2_tmp) +
                                      x[0] * W[7] * x[1] * 2.0) +
                                     x[0] * W[11] * x[2] * 2.0) +
                                    x[0] * W[19] * x[1] * 2.0) +
                                   x[0] * W[15] * x[3] * 2.0) +
                                  x[1] * W[27] * x[2]) +
                                 x[0] * W[35] * x[2] * 2.0) +
                                x[1] * W[31] * x[3]) +
                               x[1] * W[39] * x[2]) +
                              x[0] * W[51] * x[3] * 2.0) +
                             x[2] * W[47] * x[3]) +
                            x[1] * W[55] * x[3]) +
                           x[0] * W[67] * x[1] * 2.0) +
                          x[2] * W[59] * x[3]) +
                         x[1] * W[75] * x[2]) +
                        x[1] * W[79] * x[3]) +
                       x[1] * W[99] * x[2]) +
                      x[1] * W[115] * x[3]) +
                     x[0] * W[131] * x[2] * 2.0) +
                    x[1] * W[135] * x[2]) +
                   x[1] * W[147] * x[2]) +
                  x[2] * W[143] * x[3]) +
                 x[2] * W[179] * x[3]) +
                x[0] * W[195] * x[3] * 2.0) +
               x[1] * W[199] * x[3]) +
              x[2] * W[203] * x[3]) +
             x[1] * W[211] * x[3]) +
            x[2] * W[227] * x[3];
    e_varargout_2_tmp = x[0] * W[24];
    f_varargout_2_tmp = x[0] * W[36];
    g_varargout_2_tmp = x[0] * W[72];
    h_varargout_2_tmp = x[0] * W[96];
    i_varargout_2_tmp = x[0] * W[132];
    j_varargout_2_tmp = x[0] * W[144];
    k_varargout_2_tmp = x[0] * W[28];
    l_varargout_2_tmp = x[0] * W[52];
    m_varargout_2_tmp = x[0] * W[76];
    n_varargout_2_tmp = x[0] * W[112];
    o_varargout_2_tmp = x[0] * W[196];
    p_varargout_2_tmp = x[0] * W[208];
    varargout_2[4] =
            ((((((((((((((((((((((((((((((((((((-Q[4] + W[4] * d_varargout_2_tmp) +
                                               W[16] * d_varargout_2_tmp) +
                                              W[64] * d_varargout_2_tmp) +
                                             W[84] * varargout_2_tmp * 3.0) +
                                            W[104] * b_varargout_2_tmp) +
                                           W[124] * c_varargout_2_tmp) +
                                          W[152] * b_varargout_2_tmp) +
                                         W[164] * b_varargout_2_tmp) +
                                        W[220] * c_varargout_2_tmp) +
                                       W[244] * c_varargout_2_tmp) +
                                      x[0] * W[20] * x[1] * 2.0) +
                                     e_varargout_2_tmp * x[2]) +
                                    k_varargout_2_tmp * x[3]) +
                                   f_varargout_2_tmp * x[2]) +
                                  l_varargout_2_tmp * x[3]) +
                                 x[0] * W[68] * x[1] * 2.0) +
                                g_varargout_2_tmp * x[2]) +
                               x[0] * W[80] * x[1] * 2.0) +
                              m_varargout_2_tmp * x[3]) +
                             x[1] * W[88] * x[2] * 2.0) +
                            h_varargout_2_tmp * x[2]) +
                           x[1] * W[92] * x[3] * 2.0) +
                          x[1] * W[100] * x[2] * 2.0) +
                         n_varargout_2_tmp * x[3]) +
                        x[2] * W[108] * x[3]) +
                       x[1] * W[116] * x[3] * 2.0) +
                      x[2] * W[120] * x[3]) +
                     i_varargout_2_tmp * x[2]) +
                    j_varargout_2_tmp * x[2]) +
                   x[1] * W[148] * x[2] * 2.0) +
                  x[2] * W[156] * x[3]) +
                 x[2] * W[180] * x[3]) +
                o_varargout_2_tmp * x[3]) +
               p_varargout_2_tmp * x[3]) +
              x[1] * W[212] * x[3] * 2.0) +
             x[2] * W[216] * x[3]) +
            x[2] * W[228] * x[3];
    q_varargout_2_tmp = x[0] * W[25];
    r_varargout_2_tmp = x[0] * W[37];
    s_varargout_2_tmp = x[0] * W[73];
    t_varargout_2_tmp = x[0] * W[97];
    u_varargout_2_tmp = x[0] * W[133];
    v_varargout_2_tmp = x[0] * W[145];
    w_varargout_2_tmp = x[0] * W[29];
    x_varargout_2_tmp = x[0] * W[53];
    y_varargout_2_tmp = x[0] * W[77];
    ab_varargout_2_tmp = x[0] * W[113];
    bb_varargout_2_tmp = x[0] * W[197];
    cb_varargout_2_tmp = x[0] * W[209];
    varargout_2[5] =
            (((((((((((((((((((((((((((((((((((((-Q[5] + x[4]) +
                                                W[5] * d_varargout_2_tmp) +
                                               W[17] * d_varargout_2_tmp) +
                                              W[65] * d_varargout_2_tmp) +
                                             W[85] * varargout_2_tmp * 3.0) +
                                            W[105] * b_varargout_2_tmp) +
                                           W[125] * c_varargout_2_tmp) +
                                          W[153] * b_varargout_2_tmp) +
                                         W[165] * b_varargout_2_tmp) +
                                        W[221] * c_varargout_2_tmp) +
                                       W[245] * c_varargout_2_tmp) +
                                      x[0] * W[21] * x[1] * 2.0) +
                                     q_varargout_2_tmp * x[2]) +
                                    w_varargout_2_tmp * x[3]) +
                                   r_varargout_2_tmp * x[2]) +
                                  x_varargout_2_tmp * x[3]) +
                                 x[0] * W[69] * x[1] * 2.0) +
                                s_varargout_2_tmp * x[2]) +
                               x[0] * W[81] * x[1] * 2.0) +
                              y_varargout_2_tmp * x[3]) +
                             x[1] * W[89] * x[2] * 2.0) +
                            t_varargout_2_tmp * x[2]) +
                           x[1] * W[93] * x[3] * 2.0) +
                          x[1] * W[101] * x[2] * 2.0) +
                         ab_varargout_2_tmp * x[3]) +
                        x[2] * W[109] * x[3]) +
                       x[1] * W[117] * x[3] * 2.0) +
                      x[2] * W[121] * x[3]) +
                     u_varargout_2_tmp * x[2]) +
                    v_varargout_2_tmp * x[2]) +
                   x[1] * W[149] * x[2] * 2.0) +
                  x[2] * W[157] * x[3]) +
                 x[2] * W[181] * x[3]) +
                bb_varargout_2_tmp * x[3]) +
               cb_varargout_2_tmp * x[3]) +
              x[1] * W[213] * x[3] * 2.0) +
             x[2] * W[217] * x[3]) +
            x[2] * W[229] * x[3];
    db_varargout_2_tmp = x[0] * W[26];
    eb_varargout_2_tmp = x[0] * W[38];
    fb_varargout_2_tmp = x[0] * W[74];
    gb_varargout_2_tmp = x[0] * W[98];
    hb_varargout_2_tmp = x[0] * W[134];
    ib_varargout_2_tmp = x[0] * W[146];
    jb_varargout_2_tmp = x[0] * W[30];
    kb_varargout_2_tmp = x[0] * W[54];
    lb_varargout_2_tmp = x[0] * W[78];
    mb_varargout_2_tmp = x[0] * W[114];
    nb_varargout_2_tmp = x[0] * W[198];
    ob_varargout_2_tmp = x[0] * W[210];
    varargout_2[6] =
            ((((((((((((((((((((((((((((((((((((-Q[6] + W[6] * d_varargout_2_tmp) +
                                               W[18] * d_varargout_2_tmp) +
                                              W[66] * d_varargout_2_tmp) +
                                             W[86] * varargout_2_tmp * 3.0) +
                                            W[106] * b_varargout_2_tmp) +
                                           W[126] * c_varargout_2_tmp) +
                                          W[154] * b_varargout_2_tmp) +
                                         W[166] * b_varargout_2_tmp) +
                                        W[222] * c_varargout_2_tmp) +
                                       W[246] * c_varargout_2_tmp) +
                                      x[0] * W[22] * x[1] * 2.0) +
                                     db_varargout_2_tmp * x[2]) +
                                    jb_varargout_2_tmp * x[3]) +
                                   eb_varargout_2_tmp * x[2]) +
                                  kb_varargout_2_tmp * x[3]) +
                                 x[0] * W[70] * x[1] * 2.0) +
                                fb_varargout_2_tmp * x[2]) +
                               x[0] * W[82] * x[1] * 2.0) +
                              lb_varargout_2_tmp * x[3]) +
                             x[1] * W[90] * x[2] * 2.0) +
                            gb_varargout_2_tmp * x[2]) +
                           x[1] * W[94] * x[3] * 2.0) +
                          x[1] * W[102] * x[2] * 2.0) +
                         mb_varargout_2_tmp * x[3]) +
                        x[2] * W[110] * x[3]) +
                       x[1] * W[118] * x[3] * 2.0) +
                      x[2] * W[122] * x[3]) +
                     hb_varargout_2_tmp * x[2]) +
                    ib_varargout_2_tmp * x[2]) +
                   x[1] * W[150] * x[2] * 2.0) +
                  x[2] * W[158] * x[3]) +
                 x[2] * W[182] * x[3]) +
                nb_varargout_2_tmp * x[3]) +
               ob_varargout_2_tmp * x[3]) +
              x[1] * W[214] * x[3] * 2.0) +
             x[2] * W[218] * x[3]) +
            x[2] * W[230] * x[3];
    pb_varargout_2_tmp = x[0] * W[27];
    qb_varargout_2_tmp = x[0] * W[39];
    rb_varargout_2_tmp = x[0] * W[75];
    sb_varargout_2_tmp = x[0] * W[99];
    tb_varargout_2_tmp = x[0] * W[135];
    ub_varargout_2_tmp = x[0] * W[147];
    vb_varargout_2_tmp = x[0] * W[31];
    wb_varargout_2_tmp = x[0] * W[55];
    xb_varargout_2_tmp = x[0] * W[79];
    yb_varargout_2_tmp = x[0] * W[115];
    ac_varargout_2_tmp = x[0] * W[199];
    bc_varargout_2_tmp = x[0] * W[211];
    varargout_2[7] =
            ((((((((((((((((((((((((((((((((((((-Q[7] + W[7] * d_varargout_2_tmp) +
                                               W[19] * d_varargout_2_tmp) +
                                              W[67] * d_varargout_2_tmp) +
                                             W[87] * varargout_2_tmp * 3.0) +
                                            W[107] * b_varargout_2_tmp) +
                                           W[127] * c_varargout_2_tmp) +
                                          W[155] * b_varargout_2_tmp) +
                                         W[167] * b_varargout_2_tmp) +
                                        W[223] * c_varargout_2_tmp) +
                                       W[247] * c_varargout_2_tmp) +
                                      x[0] * W[23] * x[1] * 2.0) +
                                     pb_varargout_2_tmp * x[2]) +
                                    vb_varargout_2_tmp * x[3]) +
                                   qb_varargout_2_tmp * x[2]) +
                                  wb_varargout_2_tmp * x[3]) +
                                 x[0] * W[71] * x[1] * 2.0) +
                                rb_varargout_2_tmp * x[2]) +
                               x[0] * W[83] * x[1] * 2.0) +
                              xb_varargout_2_tmp * x[3]) +
                             x[1] * W[91] * x[2] * 2.0) +
                            sb_varargout_2_tmp * x[2]) +
                           x[1] * W[95] * x[3] * 2.0) +
                          x[1] * W[103] * x[2] * 2.0) +
                         yb_varargout_2_tmp * x[3]) +
                        x[2] * W[111] * x[3]) +
                       x[1] * W[119] * x[3] * 2.0) +
                      x[2] * W[123] * x[3]) +
                     tb_varargout_2_tmp * x[2]) +
                    ub_varargout_2_tmp * x[2]) +
                   x[1] * W[151] * x[2] * 2.0) +
                  x[2] * W[159] * x[3]) +
                 x[2] * W[183] * x[3]) +
                ac_varargout_2_tmp * x[3]) +
               bc_varargout_2_tmp * x[3]) +
              x[1] * W[215] * x[3] * 2.0) +
             x[2] * W[219] * x[3]) +
            x[2] * W[231] * x[3];
    cc_varargout_2_tmp = x[0] * W[44];
    dc_varargout_2_tmp = x[0] * W[56];
    ec_varargout_2_tmp = x[1] * W[108];
    fc_varargout_2_tmp = x[1] * W[120];
    gc_varargout_2_tmp = x[0] * W[140];
    hc_varargout_2_tmp = x[1] * W[156];
    ic_varargout_2_tmp = x[0] * W[176];
    jc_varargout_2_tmp = x[1] * W[180];
    kc_varargout_2_tmp = x[0] * W[200];
    lc_varargout_2_tmp = x[1] * W[216];
    mc_varargout_2_tmp = x[0] * W[224];
    nc_varargout_2_tmp = x[1] * W[228];
    varargout_2[8] =
            ((((((((((((((((((((((((((((((((((((-Q[8] + W[8] * d_varargout_2_tmp) +
                                               W[32] * d_varargout_2_tmp) +
                                              W[88] * varargout_2_tmp) +
                                             W[100] * varargout_2_tmp) +
                                            W[128] * d_varargout_2_tmp) +
                                           W[148] * varargout_2_tmp) +
                                          W[168] * b_varargout_2_tmp * 3.0) +
                                         W[188] * c_varargout_2_tmp) +
                                        W[236] * c_varargout_2_tmp) +
                                       W[248] * c_varargout_2_tmp) +
                                      e_varargout_2_tmp * x[1]) +
                                     f_varargout_2_tmp * x[1]) +
                                    x[0] * W[40] * x[2] * 2.0) +
                                   cc_varargout_2_tmp * x[3]) +
                                  dc_varargout_2_tmp * x[3]) +
                                 g_varargout_2_tmp * x[1]) +
                                h_varargout_2_tmp * x[1]) +
                               x[1] * W[104] * x[2] * 2.0) +
                              ec_varargout_2_tmp * x[3]) +
                             fc_varargout_2_tmp * x[3]) +
                            i_varargout_2_tmp * x[1]) +
                           x[0] * W[136] * x[2] * 2.0) +
                          j_varargout_2_tmp * x[1]) +
                         gc_varargout_2_tmp * x[3]) +
                        x[1] * W[152] * x[2] * 2.0) +
                       x[0] * W[160] * x[2] * 2.0) +
                      hc_varargout_2_tmp * x[3]) +
                     x[1] * W[164] * x[2] * 2.0) +
                    ic_varargout_2_tmp * x[3]) +
                   x[2] * W[172] * x[3] * 2.0) +
                  jc_varargout_2_tmp * x[3]) +
                 x[2] * W[184] * x[3] * 2.0) +
                kc_varargout_2_tmp * x[3]) +
               lc_varargout_2_tmp * x[3]) +
              mc_varargout_2_tmp * x[3]) +
             nc_varargout_2_tmp * x[3]) +
            x[2] * W[232] * x[3] * 2.0;
    e_varargout_2_tmp = x[0] * W[45];
    f_varargout_2_tmp = x[0] * W[57];
    g_varargout_2_tmp = x[1] * W[109];
    h_varargout_2_tmp = x[1] * W[121];
    i_varargout_2_tmp = x[0] * W[141];
    j_varargout_2_tmp = x[1] * W[157];
    oc_varargout_2_tmp = x[0] * W[177];
    pc_varargout_2_tmp = x[1] * W[181];
    qc_varargout_2_tmp = x[0] * W[201];
    rc_varargout_2_tmp = x[1] * W[217];
    sc_varargout_2_tmp = x[0] * W[225];
    tc_varargout_2_tmp = x[1] * W[229];
    varargout_2[9] =
            ((((((((((((((((((((((((((((((((((((-Q[9] + W[9] * d_varargout_2_tmp) +
                                               W[33] * d_varargout_2_tmp) +
                                              W[89] * varargout_2_tmp) +
                                             W[101] * varargout_2_tmp) +
                                            W[129] * d_varargout_2_tmp) +
                                           W[149] * varargout_2_tmp) +
                                          W[169] * b_varargout_2_tmp * 3.0) +
                                         W[189] * c_varargout_2_tmp) +
                                        W[237] * c_varargout_2_tmp) +
                                       W[249] * c_varargout_2_tmp) +
                                      q_varargout_2_tmp * x[1]) +
                                     r_varargout_2_tmp * x[1]) +
                                    x[0] * W[41] * x[2] * 2.0) +
                                   e_varargout_2_tmp * x[3]) +
                                  f_varargout_2_tmp * x[3]) +
                                 s_varargout_2_tmp * x[1]) +
                                t_varargout_2_tmp * x[1]) +
                               x[1] * W[105] * x[2] * 2.0) +
                              g_varargout_2_tmp * x[3]) +
                             h_varargout_2_tmp * x[3]) +
                            u_varargout_2_tmp * x[1]) +
                           x[0] * W[137] * x[2] * 2.0) +
                          v_varargout_2_tmp * x[1]) +
                         i_varargout_2_tmp * x[3]) +
                        x[1] * W[153] * x[2] * 2.0) +
                       x[0] * W[161] * x[2] * 2.0) +
                      j_varargout_2_tmp * x[3]) +
                     x[1] * W[165] * x[2] * 2.0) +
                    oc_varargout_2_tmp * x[3]) +
                   x[2] * W[173] * x[3] * 2.0) +
                  pc_varargout_2_tmp * x[3]) +
                 x[2] * W[185] * x[3] * 2.0) +
                qc_varargout_2_tmp * x[3]) +
               rc_varargout_2_tmp * x[3]) +
              sc_varargout_2_tmp * x[3]) +
             tc_varargout_2_tmp * x[3]) +
            x[2] * W[233] * x[3] * 2.0;
    q_varargout_2_tmp = x[0] * W[46];
    r_varargout_2_tmp = x[0] * W[58];
    s_varargout_2_tmp = x[1] * W[110];
    t_varargout_2_tmp = x[1] * W[122];
    u_varargout_2_tmp = x[0] * W[142];
    v_varargout_2_tmp = x[1] * W[158];
    uc_varargout_2_tmp = x[0] * W[178];
    vc_varargout_2_tmp = x[1] * W[182];
    wc_varargout_2_tmp = x[0] * W[202];
    xc_varargout_2_tmp = x[1] * W[218];
    yc_varargout_2_tmp = x[0] * W[226];
    ad_varargout_2_tmp = x[1] * W[230];
    varargout_2[10] =
            (((((((((((((((((((((((((((((((((((((-Q[10] + x[4]) +
                                                W[10] * d_varargout_2_tmp) +
                                               W[34] * d_varargout_2_tmp) +
                                              W[90] * varargout_2_tmp) +
                                             W[102] * varargout_2_tmp) +
                                            W[130] * d_varargout_2_tmp) +
                                           W[150] * varargout_2_tmp) +
                                          W[170] * b_varargout_2_tmp * 3.0) +
                                         W[190] * c_varargout_2_tmp) +
                                        W[238] * c_varargout_2_tmp) +
                                       W[250] * c_varargout_2_tmp) +
                                      db_varargout_2_tmp * x[1]) +
                                     eb_varargout_2_tmp * x[1]) +
                                    x[0] * W[42] * x[2] * 2.0) +
                                   q_varargout_2_tmp * x[3]) +
                                  r_varargout_2_tmp * x[3]) +
                                 fb_varargout_2_tmp * x[1]) +
                                gb_varargout_2_tmp * x[1]) +
                               x[1] * W[106] * x[2] * 2.0) +
                              s_varargout_2_tmp * x[3]) +
                             t_varargout_2_tmp * x[3]) +
                            hb_varargout_2_tmp * x[1]) +
                           x[0] * W[138] * x[2] * 2.0) +
                          ib_varargout_2_tmp * x[1]) +
                         u_varargout_2_tmp * x[3]) +
                        x[1] * W[154] * x[2] * 2.0) +
                       x[0] * W[162] * x[2] * 2.0) +
                      v_varargout_2_tmp * x[3]) +
                     x[1] * W[166] * x[2] * 2.0) +
                    uc_varargout_2_tmp * x[3]) +
                   x[2] * W[174] * x[3] * 2.0) +
                  vc_varargout_2_tmp * x[3]) +
                 x[2] * W[186] * x[3] * 2.0) +
                wc_varargout_2_tmp * x[3]) +
               xc_varargout_2_tmp * x[3]) +
              yc_varargout_2_tmp * x[3]) +
             ad_varargout_2_tmp * x[3]) +
            x[2] * W[234] * x[3] * 2.0;
    db_varargout_2_tmp = x[0] * W[47];
    eb_varargout_2_tmp = x[0] * W[59];
    fb_varargout_2_tmp = x[1] * W[111];
    gb_varargout_2_tmp = x[1] * W[123];
    hb_varargout_2_tmp = x[0] * W[143];
    ib_varargout_2_tmp = x[1] * W[159];
    bd_varargout_2_tmp = x[0] * W[179];
    cd_varargout_2_tmp = x[1] * W[183];
    dd_varargout_2_tmp = x[0] * W[203];
    ed_varargout_2_tmp = x[1] * W[219];
    fd_varargout_2_tmp = x[0] * W[227];
    gd_varargout_2_tmp = x[1] * W[231];
    varargout_2[11] =
            ((((((((((((((((((((((((((((((((((((-Q[11] + W[11] * d_varargout_2_tmp) +
                                               W[35] * d_varargout_2_tmp) +
                                              W[91] * varargout_2_tmp) +
                                             W[103] * varargout_2_tmp) +
                                            W[131] * d_varargout_2_tmp) +
                                           W[151] * varargout_2_tmp) +
                                          W[171] * b_varargout_2_tmp * 3.0) +
                                         W[191] * c_varargout_2_tmp) +
                                        W[239] * c_varargout_2_tmp) +
                                       W[251] * c_varargout_2_tmp) +
                                      pb_varargout_2_tmp * x[1]) +
                                     qb_varargout_2_tmp * x[1]) +
                                    x[0] * W[43] * x[2] * 2.0) +
                                   db_varargout_2_tmp * x[3]) +
                                  eb_varargout_2_tmp * x[3]) +
                                 rb_varargout_2_tmp * x[1]) +
                                sb_varargout_2_tmp * x[1]) +
                               x[1] * W[107] * x[2] * 2.0) +
                              fb_varargout_2_tmp * x[3]) +
                             gb_varargout_2_tmp * x[3]) +
                            tb_varargout_2_tmp * x[1]) +
                           x[0] * W[139] * x[2] * 2.0) +
                          ub_varargout_2_tmp * x[1]) +
                         hb_varargout_2_tmp * x[3]) +
                        x[1] * W[155] * x[2] * 2.0) +
                       x[0] * W[163] * x[2] * 2.0) +
                      ib_varargout_2_tmp * x[3]) +
                     x[1] * W[167] * x[2] * 2.0) +
                    bd_varargout_2_tmp * x[3]) +
                   x[2] * W[175] * x[3] * 2.0) +
                  cd_varargout_2_tmp * x[3]) +
                 x[2] * W[187] * x[3] * 2.0) +
                dd_varargout_2_tmp * x[3]) +
               ed_varargout_2_tmp * x[3]) +
              fd_varargout_2_tmp * x[3]) +
             gd_varargout_2_tmp * x[3]) +
            x[2] * W[235] * x[3] * 2.0;
    varargout_2[12] =
            ((((((((((((((((((((((((((((((((((((-Q[12] + W[12] * d_varargout_2_tmp) +
                                               W[48] * d_varargout_2_tmp) +
                                              W[92] * varargout_2_tmp) +
                                             W[116] * varargout_2_tmp) +
                                            W[172] * b_varargout_2_tmp) +
                                           W[184] * b_varargout_2_tmp) +
                                          W[192] * d_varargout_2_tmp) +
                                         W[212] * varargout_2_tmp) +
                                        W[232] * b_varargout_2_tmp) +
                                       W[252] * c_varargout_2_tmp * 3.0) +
                                      k_varargout_2_tmp * x[1]) +
                                     cc_varargout_2_tmp * x[2]) +
                                    l_varargout_2_tmp * x[1]) +
                                   dc_varargout_2_tmp * x[2]) +
                                  x[0] * W[60] * x[3] * 2.0) +
                                 m_varargout_2_tmp * x[1]) +
                                n_varargout_2_tmp * x[1]) +
                               ec_varargout_2_tmp * x[2]) +
                              fc_varargout_2_tmp * x[2]) +
                             x[1] * W[124] * x[3] * 2.0) +
                            gc_varargout_2_tmp * x[2]) +
                           hc_varargout_2_tmp * x[2]) +
                          ic_varargout_2_tmp * x[2]) +
                         jc_varargout_2_tmp * x[2]) +
                        o_varargout_2_tmp * x[1]) +
                       x[2] * W[188] * x[3] * 2.0) +
                      kc_varargout_2_tmp * x[2]) +
                     p_varargout_2_tmp * x[1]) +
                    x[0] * W[204] * x[3] * 2.0) +
                   lc_varargout_2_tmp * x[2]) +
                  mc_varargout_2_tmp * x[2]) +
                 x[1] * W[220] * x[3] * 2.0) +
                nc_varargout_2_tmp * x[2]) +
               x[0] * W[240] * x[3] * 2.0) +
              x[2] * W[236] * x[3] * 2.0) +
             x[1] * W[244] * x[3] * 2.0) +
            x[2] * W[248] * x[3] * 2.0;
    varargout_2[13] =
            ((((((((((((((((((((((((((((((((((((-Q[13] + W[13] * d_varargout_2_tmp) +
                                               W[49] * d_varargout_2_tmp) +
                                              W[93] * varargout_2_tmp) +
                                             W[117] * varargout_2_tmp) +
                                            W[173] * b_varargout_2_tmp) +
                                           W[185] * b_varargout_2_tmp) +
                                          W[193] * d_varargout_2_tmp) +
                                         W[213] * varargout_2_tmp) +
                                        W[233] * b_varargout_2_tmp) +
                                       W[253] * c_varargout_2_tmp * 3.0) +
                                      w_varargout_2_tmp * x[1]) +
                                     e_varargout_2_tmp * x[2]) +
                                    x_varargout_2_tmp * x[1]) +
                                   f_varargout_2_tmp * x[2]) +
                                  x[0] * W[61] * x[3] * 2.0) +
                                 y_varargout_2_tmp * x[1]) +
                                ab_varargout_2_tmp * x[1]) +
                               g_varargout_2_tmp * x[2]) +
                              h_varargout_2_tmp * x[2]) +
                             x[1] * W[125] * x[3] * 2.0) +
                            i_varargout_2_tmp * x[2]) +
                           j_varargout_2_tmp * x[2]) +
                          oc_varargout_2_tmp * x[2]) +
                         pc_varargout_2_tmp * x[2]) +
                        bb_varargout_2_tmp * x[1]) +
                       x[2] * W[189] * x[3] * 2.0) +
                      qc_varargout_2_tmp * x[2]) +
                     cb_varargout_2_tmp * x[1]) +
                    x[0] * W[205] * x[3] * 2.0) +
                   rc_varargout_2_tmp * x[2]) +
                  sc_varargout_2_tmp * x[2]) +
                 x[1] * W[221] * x[3] * 2.0) +
                tc_varargout_2_tmp * x[2]) +
               x[0] * W[241] * x[3] * 2.0) +
              x[2] * W[237] * x[3] * 2.0) +
             x[1] * W[245] * x[3] * 2.0) +
            x[2] * W[249] * x[3] * 2.0;
    varargout_2[14] =
            ((((((((((((((((((((((((((((((((((((-Q[14] + W[14] * d_varargout_2_tmp) +
                                               W[50] * d_varargout_2_tmp) +
                                              W[94] * varargout_2_tmp) +
                                             W[118] * varargout_2_tmp) +
                                            W[174] * b_varargout_2_tmp) +
                                           W[186] * b_varargout_2_tmp) +
                                          W[194] * d_varargout_2_tmp) +
                                         W[214] * varargout_2_tmp) +
                                        W[234] * b_varargout_2_tmp) +
                                       W[254] * c_varargout_2_tmp * 3.0) +
                                      jb_varargout_2_tmp * x[1]) +
                                     q_varargout_2_tmp * x[2]) +
                                    kb_varargout_2_tmp * x[1]) +
                                   r_varargout_2_tmp * x[2]) +
                                  x[0] * W[62] * x[3] * 2.0) +
                                 lb_varargout_2_tmp * x[1]) +
                                mb_varargout_2_tmp * x[1]) +
                               s_varargout_2_tmp * x[2]) +
                              t_varargout_2_tmp * x[2]) +
                             x[1] * W[126] * x[3] * 2.0) +
                            u_varargout_2_tmp * x[2]) +
                           v_varargout_2_tmp * x[2]) +
                          uc_varargout_2_tmp * x[2]) +
                         vc_varargout_2_tmp * x[2]) +
                        nb_varargout_2_tmp * x[1]) +
                       x[2] * W[190] * x[3] * 2.0) +
                      wc_varargout_2_tmp * x[2]) +
                     ob_varargout_2_tmp * x[1]) +
                    x[0] * W[206] * x[3] * 2.0) +
                   xc_varargout_2_tmp * x[2]) +
                  yc_varargout_2_tmp * x[2]) +
                 x[1] * W[222] * x[3] * 2.0) +
                ad_varargout_2_tmp * x[2]) +
               x[0] * W[242] * x[3] * 2.0) +
              x[2] * W[238] * x[3] * 2.0) +
             x[1] * W[246] * x[3] * 2.0) +
            x[2] * W[250] * x[3] * 2.0;
    varargout_2[15] =
            (((((((((((((((((((((((((((((((((((((-Q[15] + x[4]) +
                                                W[15] * d_varargout_2_tmp) +
                                               W[51] * d_varargout_2_tmp) +
                                              W[95] * varargout_2_tmp) +
                                             W[119] * varargout_2_tmp) +
                                            W[175] * b_varargout_2_tmp) +
                                           W[187] * b_varargout_2_tmp) +
                                          W[195] * d_varargout_2_tmp) +
                                         W[215] * varargout_2_tmp) +
                                        W[235] * b_varargout_2_tmp) +
                                       W[255] * c_varargout_2_tmp * 3.0) +
                                      vb_varargout_2_tmp * x[1]) +
                                     db_varargout_2_tmp * x[2]) +
                                    wb_varargout_2_tmp * x[1]) +
                                   eb_varargout_2_tmp * x[2]) +
                                  x[0] * W[63] * x[3] * 2.0) +
                                 xb_varargout_2_tmp * x[1]) +
                                yb_varargout_2_tmp * x[1]) +
                               fb_varargout_2_tmp * x[2]) +
                              gb_varargout_2_tmp * x[2]) +
                             x[1] * W[127] * x[3] * 2.0) +
                            hb_varargout_2_tmp * x[2]) +
                           ib_varargout_2_tmp * x[2]) +
                          bd_varargout_2_tmp * x[2]) +
                         cd_varargout_2_tmp * x[2]) +
                        ac_varargout_2_tmp * x[1]) +
                       x[2] * W[191] * x[3] * 2.0) +
                      dd_varargout_2_tmp * x[2]) +
                     bc_varargout_2_tmp * x[1]) +
                    x[0] * W[207] * x[3] * 2.0) +
                   ed_varargout_2_tmp * x[2]) +
                  fd_varargout_2_tmp * x[2]) +
                 x[1] * W[223] * x[3] * 2.0) +
                gd_varargout_2_tmp * x[2]) +
               x[0] * W[243] * x[3] * 2.0) +
              x[2] * W[239] * x[3] * 2.0) +
             x[1] * W[247] * x[3] * 2.0) +
            x[2] * W[251] * x[3] * 2.0;
    varargout_2[16] = x[0];
    varargout_2[17] = x[1];
    varargout_2[18] = x[2];
    varargout_2[19] = x[3];
}

//
// Arguments    : const double W[256]
//                const double Q[16]
//                const double x_init[5]
//                double xx[5]
// Return Type  : void
//
void QPEP_fsolve_func(const double W[256], const double Q[16],
                      const double x_init[5], double xx[5])
{
    double augJacobian[45];
    double jacob[20];
    double varargout_2[20];
    double rhs[9];
    double dx[5];
    double gradf[5];
    double fval[4];
    double varargout_1[4];
    double absxk;
    double b_gamma;
    double b_scale;
    double c;
    double funDiff;
    double normGradF;
    double relFactor;
    double resnorm;
    double scale;
    double t;
    int exitflag;
    int funcCount;
    int i;
    int iter;
    int iy0;
    boolean_T exitg1;
    boolean_T stepSuccessful;
    funcCount = 1;
    funDiff = rtInf;
    iter = 0;
    anon(W, Q, x_init, varargout_1, jacob);
    fval[0] = varargout_1[0];
    fval[1] = varargout_1[1];
    fval[2] = varargout_1[2];
    fval[3] = varargout_1[3];
    for (i = 0; i < 5; i++) {
        dx[i] = rtInf;
        xx[i] = x_init[i];
        exitflag = i << 2;
        augJacobian[9 * i] = jacob[exitflag];
        augJacobian[9 * i + 1] = jacob[exitflag + 1];
        augJacobian[9 * i + 2] = jacob[exitflag + 2];
        augJacobian[9 * i + 3] = jacob[exitflag + 3];
    }
    resnorm =
            ((varargout_1[0] * varargout_1[0] + varargout_1[1] * varargout_1[1]) +
             varargout_1[2] * varargout_1[2]) +
            varargout_1[3] * varargout_1[3];
    b_gamma = 0.01;
    for (i = 0; i < 5; i++) {
        iy0 = 9 * (i + 1);
        for (exitflag = 0; exitflag < 5; exitflag++) {
            augJacobian[(iy0 + exitflag) - 5] = 0.0;
        }
        augJacobian[(i + 9 * i) + 4] = 0.1;
        exitflag = 9 * i;
        iy0 = i << 2;
        jacob[iy0] = augJacobian[exitflag];
        jacob[iy0 + 1] = augJacobian[exitflag + 1];
        jacob[iy0 + 2] = augJacobian[exitflag + 2];
        jacob[iy0 + 3] = augJacobian[exitflag + 3];
        gradf[i] = 0.0;
    }
    for (iy0 = 0; iy0 <= 16; iy0 += 4) {
        c = 0.0;
        i = iy0 + 4;
        for (exitflag = iy0 + 1; exitflag <= i; exitflag++) {
            c += jacob[exitflag - 1] * fval[(exitflag - iy0) - 1];
        }
        exitflag = iy0 >> 2;
        gradf[exitflag] += c;
    }
    c = 0.0;
    stepSuccessful = true;
    normGradF = 0.0;
    for (exitflag = 0; exitflag < 5; exitflag++) {
        scale = std::abs(gradf[exitflag]);
        if (std::isnan(scale) || (scale > c)) {
            c = scale;
        }
        if (std::isnan(scale) || (scale > normGradF)) {
            normGradF = scale;
        }
    }
    relFactor = std::fmax(c, 1.4901161193847656E-8);
    if (normGradF <= 1.0E-12 * relFactor) {
        exitflag = 1;
    } else {
        c = 0.0;
        scale = 3.3121686421112381E-170;
        normGradF = 0.0;
        b_scale = 3.3121686421112381E-170;
        for (exitflag = 0; exitflag < 5; exitflag++) {
            if (rtInf > scale) {
                c = c * 0.0 * 0.0 + 1.0;
                scale = rtInf;
            } else {
                c = rtNaN;
            }
            absxk = std::abs(x_init[exitflag]);
            if (absxk > b_scale) {
                t = b_scale / absxk;
                normGradF = normGradF * t * t + 1.0;
                b_scale = absxk;
            } else {
                t = absxk / b_scale;
                normGradF += t * t;
            }
        }
        c = scale * std::sqrt(c);
        normGradF = b_scale * std::sqrt(normGradF);
        if (c < 1.0E-6 * (normGradF + 1.4901161193847656E-8)) {
            exitflag = 4;
        } else {
            exitflag = -5;
        }
    }
    exitg1 = false;
    while ((!exitg1) && (exitflag == -5)) {
        boolean_T evalOK;
        boolean_T guard1{false};
        rhs[0] = -fval[0];
        rhs[1] = -fval[1];
        rhs[2] = -fval[2];
        rhs[3] = -fval[3];
        for (exitflag = 0; exitflag < 5; exitflag++) {
            rhs[exitflag + 4] = 0.0;
        }
        coder::optim::coder::levenbergMarquardt::linearLeastSquares(augJacobian,
                                                                    rhs, dx);
        for (i = 0; i < 5; i++) {
            gradf[i] = xx[i] + dx[i];
        }
        anon(W, Q, gradf, varargout_1, varargout_2);
        for (i = 0; i < 5; i++) {
            exitflag = i << 2;
            augJacobian[9 * i] = varargout_2[exitflag];
            augJacobian[9 * i + 1] = varargout_2[exitflag + 1];
            augJacobian[9 * i + 2] = varargout_2[exitflag + 2];
            augJacobian[9 * i + 3] = varargout_2[exitflag + 3];
        }
        evalOK = ((!std::isinf(varargout_1[0])) && (!std::isnan(varargout_1[0])));
        if ((!evalOK) ||
            (std::isinf(varargout_1[1]) || std::isnan(varargout_1[1]))) {
            evalOK = false;
        }
        if ((!evalOK) ||
            (std::isinf(varargout_1[2]) || std::isnan(varargout_1[2]))) {
            evalOK = false;
        }
        c = ((varargout_1[0] * varargout_1[0] + varargout_1[1] * varargout_1[1]) +
             varargout_1[2] * varargout_1[2]) +
            varargout_1[3] * varargout_1[3];
        if ((!evalOK) ||
            (std::isinf(varargout_1[3]) || std::isnan(varargout_1[3]))) {
            evalOK = false;
        }
        funcCount++;
        guard1 = false;
        if ((c < resnorm) && evalOK) {
            iter++;
            funDiff = std::abs(c - resnorm) / resnorm;
            fval[0] = varargout_1[0];
            fval[1] = varargout_1[1];
            fval[2] = varargout_1[2];
            fval[3] = varargout_1[3];
            resnorm = c;
            for (i = 0; i < 5; i++) {
                exitflag = 9 * i;
                iy0 = i << 2;
                jacob[iy0] = augJacobian[exitflag];
                jacob[iy0 + 1] = augJacobian[exitflag + 1];
                jacob[iy0 + 2] = augJacobian[exitflag + 2];
                jacob[iy0 + 3] = augJacobian[exitflag + 3];
            }
            evalOK = true;
            for (i = 0; i < 20; i++) {
                if ((!evalOK) || (std::isinf(jacob[i]) || std::isnan(jacob[i]))) {
                    evalOK = false;
                }
            }
            if (evalOK) {
                for (i = 0; i < 5; i++) {
                    xx[i] = gradf[i];
                }
                if (stepSuccessful) {
                    b_gamma *= 0.1;
                }
                stepSuccessful = true;
                guard1 = true;
            } else {
                exitg1 = true;
            }
        } else {
            b_gamma *= 10.0;
            stepSuccessful = false;
            for (i = 0; i < 5; i++) {
                exitflag = i << 2;
                augJacobian[9 * i] = jacob[exitflag];
                augJacobian[9 * i + 1] = jacob[exitflag + 1];
                augJacobian[9 * i + 2] = jacob[exitflag + 2];
                augJacobian[9 * i + 3] = jacob[exitflag + 3];
            }
            guard1 = true;
        }
        if (guard1) {
            c = std::sqrt(b_gamma);
            for (i = 0; i < 5; i++) {
                iy0 = 9 * (i + 1);
                for (exitflag = 0; exitflag < 5; exitflag++) {
                    augJacobian[(iy0 + exitflag) - 5] = 0.0;
                }
                augJacobian[(i + 9 * i) + 4] = c;
                gradf[i] = 0.0;
            }
            for (iy0 = 0; iy0 <= 16; iy0 += 4) {
                c = 0.0;
                i = iy0 + 4;
                for (exitflag = iy0 + 1; exitflag <= i; exitflag++) {
                    c += jacob[exitflag - 1] * fval[(exitflag - iy0) - 1];
                }
                exitflag = iy0 >> 2;
                gradf[exitflag] += c;
            }
            normGradF = 0.0;
            for (exitflag = 0; exitflag < 5; exitflag++) {
                c = std::abs(gradf[exitflag]);
                if (std::isnan(c) || (c > normGradF)) {
                    normGradF = c;
                }
            }
            if (normGradF <= 1.0E-12 * relFactor) {
                exitflag = 1;
            } else if (funcCount >= 10000) {
                exitflag = 0;
            } else if (iter >= 10000) {
                exitflag = 0;
            } else {
                c = 0.0;
                scale = 3.3121686421112381E-170;
                normGradF = 0.0;
                b_scale = 3.3121686421112381E-170;
                for (exitflag = 0; exitflag < 5; exitflag++) {
                    absxk = std::abs(dx[exitflag]);
                    if (absxk > scale) {
                        t = scale / absxk;
                        c = c * t * t + 1.0;
                        scale = absxk;
                    } else {
                        t = absxk / scale;
                        c += t * t;
                    }
                    absxk = std::abs(xx[exitflag]);
                    if (absxk > b_scale) {
                        t = b_scale / absxk;
                        normGradF = normGradF * t * t + 1.0;
                        b_scale = absxk;
                    } else {
                        t = absxk / b_scale;
                        normGradF += t * t;
                    }
                }
                c = scale * std::sqrt(c);
                normGradF = b_scale * std::sqrt(normGradF);
                if (c < 1.0E-6 * (normGradF + 1.4901161193847656E-8)) {
                    exitflag = 4;
                    if (!stepSuccessful) {
                        iter++;
                    }
                } else if (funDiff <= 1.0E-8) {
                    exitflag = 3;
                } else {
                    exitflag = -5;
                }
            }
            if (exitflag != -5) {
                exitg1 = true;
            }
        }
    }
}


Eigen::VectorXd absvec(const Eigen::VectorXd vec)
{
    Eigen::VectorXd t(vec);
    for(int i = 0; i < vec.size(); ++i)
    {
        t(i) = std::fabs(t(i));
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
        if(err < 1e-16 || (j > 0 && grad.norm() < 1e-16))
        {
            break;
        }
        else if(err >= 1e-16 && err < last_err / 2.0)
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



struct QPEP_runtime QPEP_lm_fsolve(Eigen::Matrix3d& R,
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
                                   const Eigen::MatrixXd& W,
                                   const Eigen::MatrixXd& Q,
                                   const struct QPEP_runtime& stat_)
{
#ifdef __linux__
    return QPEP_lm_single(R, t, X, q0, max_iter, mu, eq_Jacob_func, t_func, coef_f_q_sym, coefs_tq, pinvG, stat_);
#endif
    clock_t time1 = clock();
    struct QPEP_runtime stat = stat_;
    double x_init[5] = {};
    x_init[0] = q0(0);
    x_init[1] = q0(1);
    x_init[2] = q0(2);
    x_init[3] = q0(3);
    double W_array[256], Q_array[16], xx[5];
    for(int i = 0; i < 4; ++i)
        for(int j = 0; j < 64; ++j)
            W_array[(i - 1) * 64 + j] = W(i, j);
    for(int i = 0; i < 4; ++i)
        for(int j = 0; j < 4; ++j)
            Q_array[(i - 1) * 4 + j] = Q(i, j);
    QPEP_fsolve_func(W_array, Q_array, x_init, xx);
    Eigen::Vector4d qq0;
    qq0(0) = xx[0];
    qq0(1) = xx[1];
    qq0(2) = xx[2];
    qq0(3) = xx[3];
    Eigen::Quaterniond qq(qq0);
    Eigen::Matrix3d aaa;
    Eigen::Matrix<double, 3, 10> bbb;
    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
            aaa(i, j) = pinvG(i, j);

    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 10; ++j)
            bbb(i, j) = coefs_tq(i, j);
    t = t_func(aaa, bbb, qq);
    R = q2R(qq);
    X << R, t, Eigen::Vector3d::Zero(3).transpose(), 1.0;
    clock_t time2 = clock();
    stat.timeLM = (time2 - time1) / double(CLOCKS_PER_SEC);
    return stat;
}

