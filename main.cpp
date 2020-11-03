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
#include <opencv2/core/types.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "QPEP.h"
#include <fstream>
//#include <GL/glut.h>
#include "omp.h"




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

void readPnPdata(std::string filename,
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


void test_generateProjectedPoints()
{
    Eigen::Matrix3d R0;
    Eigen::Vector3d t0;
    Eigen::Matrix3d K;
    K.setZero();
    std::vector<Eigen::Vector3d> world_pt0;
    std::vector<Eigen::Vector2d> image_pt0;
    const std::string filename = "/Users/zarathustra/git/LibQPEP/data/pnp_data1.txt";
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

Eigen::Matrix3d R0;
Eigen::Vector3d t0;
Eigen::Matrix3d K;
std::vector<Eigen::Vector3d> world_pt0;
std::vector<Eigen::Vector2d> image_pt0;
const std::string filename = "/Users/zarathustra/git/LibQPEP/data/pnp_data-500pt-1.txt";

void test_pnp_WQD()
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

    std::cout << "True X: " << std::endl << XX << std::endl;
    std::cout << "QPEP X: " << std::endl << X << std::endl << std::endl;



//    std::vector<cv::Point2d> image_pt0_cv = Vector2dToPoint2d(image_pt0);
//    std::vector<cv::Point3d> world_pt0_cv = Vector3dToPoint3d(world_pt0);
//    cv::Mat intrinsics, distCoeffs;
//    distCoeffs = cv::Mat::zeros(4, 1, CV_64FC1);
//    for (int i = 0; i < 3; i++)
//        distCoeffs.at<double>(i, 0) = 0.0;
//    intrinsics.create(3, 3, CV_64FC1);
//    intrinsics.at<double>(0, 0) = K(0, 0);
//    intrinsics.at<double>(1, 1) = K(1, 1);
//    intrinsics.at<double>(0, 2) = K(0, 2);
//    intrinsics.at<double>(1, 2) = K(1, 2);
//    intrinsics.at<double>(2, 2) = 1;

//    cv::Mat rvec, tvec;
//    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_EPNP);
//    cv::Mat rotation_matrix = cv::Mat(3, 3, CV_64FC1, cv::Scalar::all(0));
//    cv::Rodrigues(rvec, rotation_matrix);
//    Eigen::Matrix3d Rot;
//    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
//            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
//            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
//    Eigen::Vector3d Trans;
//    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
//    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
//    std::cout << "Opencv EPnP X: " << std::endl << XX << std::endl;
//
//
//    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_P3P);
//    cv::Rodrigues(rvec, rotation_matrix);
//    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
//            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
//            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
//    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
//    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
//    std::cout << "Opencv P3P X: " << std::endl << XX << std::endl;
//
//
//    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_DLS);
//    cv::Rodrigues(rvec, rotation_matrix);
//    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
//            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
//            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
//    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
//    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
//    std::cout << "Opencv DLS X: " << std::endl << XX << std::endl;
//
//
//    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_UPNP);
//    cv::Rodrigues(rvec, rotation_matrix);
//    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
//            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
//            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
//    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
//    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
//    std::cout << "Opencv UPnP X: " << std::endl << XX << std::endl;
//
//
//    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_AP3P);
//    cv::Rodrigues(rvec, rotation_matrix);
//    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
//            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
//            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
//    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
//    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
//    std::cout << "Opencv AP3P X: " << std::endl << XX << std::endl;
//
//
//    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_IPPE);
//    cv::Rodrigues(rvec, rotation_matrix);
//    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
//            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
//            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
//    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
//    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
//    std::cout << "Opencv IPPE X: " << std::endl << XX << std::endl;
//
//
//    cv::solvePnP(cv::Mat(world_pt0_cv), cv::Mat(image_pt0_cv), intrinsics, distCoeffs, rvec, tvec, false, cv::SOLVEPNP_IPPE_SQUARE);
//    cv::Rodrigues(rvec, rotation_matrix);
//    Rot << rotation_matrix.at<double>(0, 0), rotation_matrix.at<double>(0, 1), rotation_matrix.at<double>(0, 2),
//            rotation_matrix.at<double>(1, 0), rotation_matrix.at<double>(1, 1), rotation_matrix.at<double>(1, 2),
//            rotation_matrix.at<double>(2, 0), rotation_matrix.at<double>(2, 1), rotation_matrix.at<double>(2, 2);
//    Trans << tvec.at<double>(0), tvec.at<double>(1), tvec.at<double>(2);
//    XX << Rot, Trans, Eigen::Vector3d::Zero(3).transpose(), 1.0;
//    std::cout << "Opencv IPPE SQUARE X: " << std::endl << XX << std::endl;
}



void readpTopdata(std::string filename,
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



std::vector<Eigen::Vector3d> rr0;
std::vector<Eigen::Vector3d> bb0;
std::vector<Eigen::Vector3d> nv0;
const std::string filename2 = "/Users/zarathustra/git/LibQPEP/data/pTop_data-2500pt-1.txt";

void test_pTop_WQD() {
    if(rr0.size() < 3) {
        readpTopdata(filename2, R0, t0, rr0, bb0, nv0);
    }
    Eigen::Matrix4d XX;
    XX << R0, t0, Eigen::Vector3d::Zero(3).transpose(), 1.0;

    Eigen::Matrix<double, 4, 64> W;
    Eigen::Matrix<double, 4, 4> Q;
    Eigen::Matrix<double, 3, 28> D;
    Eigen::Matrix<double, 4, 24> coef_f_q_sym;
    Eigen::Matrix<double, 1, 85> coef_J_pure;
    Eigen::Matrix<double, 3, 11> coefs_tq;
    Eigen::Matrix<double, 3, 3> pinvG;
    pTop_WQD(W, Q, D, coef_f_q_sym, coef_J_pure, coefs_tq, pinvG, rr0, bb0, nv0);

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

//    Eigen::Vector4d q0 = R2q(R);
//
//    stat = QPEP_lm_single(R, t, X, q0, 100, 5e-2,
//                   reinterpret_cast<eq_Jacob_func_handle>(eq_Jacob_pTop_func),
//                   reinterpret_cast<t_func_handle>(t_pTop_func),
//                   coef_f_q_sym, coefs_tq, pinvG, stat);
//
//    std::cout << "True X: " << std::endl << XX << std::endl;
//    std::cout << "QPEP X: " << std::endl << X << std::endl << std::endl;
}


int plot();

int main(int argc,char ** argv) {
//    glutInit(&argc, argv);
//    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
//    glutInitWindowSize(400, 400);
//    glutInitWindowPosition(200, 200);
//
//    plot();


    std::cout.precision(16);
//    test_generateProjectedPoints();

    double time = 0.0;
    clock_t time1 = clock();
    double loops = 10.0;

    for(int i = 0; i < (int) loops; ++i)
    {
//        test_pnp_WQD();
        test_pTop_WQD();
    }
    clock_t time2 = clock();
    time = time2 - time1;
    std::cout << "Time: " << time / loops / double(CLOCKS_PER_SEC) << std::endl;


    return 0;
}






































//GLfloat ctrlpoints[5][5][3] = {{{-3,0,0},{-1,1,0},{0,0,0},{1,-1,0},{3,0,0}},
//                               {{-3,0,-1},{-1,1,-1},{0,0,-1},{1,-1,-1},{3,0,-1}},
//                               {{-3,0,-3},{-1,1,-3},{0,0,-3},{1,-1,-3},{3,0,-3}},
//                               {{-3,0,-3},{-1,1,-3},{0,0,-3},{1,-1,-3},{3,0,-3}},
//                               {{-3,0,-4},{-1,1,-4},{0,0,-4},{1,-1,-4},{3,0,-4}}};
//
//GLfloat mat_ambient[] = {0.1,0.1,0.1,1.0};
//GLfloat mat_diffuse[] = {1.0,0.6,0.0,1.0};
//GLfloat mat_specular[] = {1.0,1.0,1.0,1.0};
//
//GLfloat light_ambient[] = {0.1,0.1,0.1,1.0};
////GLfloat light_diffuse[] = {1.0,1.0,1.0,0.0};
//GLfloat light_diffuse[] = {0.2,1.0,0.2,0.0};
//GLfloat light_specular[] = {1.0,1.0,1.0,0.0};
//GLfloat light_position[] = {2.0,23.0,-4.0,1.0};
//
//void myInit(void)
//{
//    glClearColor(1.0,1.0,1.0,0.0);//设置背景色
//
//    /*为光照模型指定材质参数*/
//    glMaterialfv(GL_FRONT,GL_AMBIENT,mat_ambient);
//    glMaterialfv(GL_FRONT,GL_DIFFUSE,mat_diffuse);
//    glMaterialfv(GL_FRONT,GL_SPECULAR,mat_specular);
//    glMaterialf(GL_FRONT,GL_SHININESS,60.0);
//
//    /*设置光源参数*/
//    glLightfv(GL_LIGHT0,GL_AMBIENT,light_ambient);
//    glLightfv(GL_LIGHT0,GL_DIFFUSE,light_diffuse);
//    glLightfv(GL_LIGHT0,GL_SPECULAR,light_specular);
//    glLightfv(GL_LIGHT0,GL_POSITION,light_position);
//
//    glEnable(GL_LIGHTING);
//    glEnable(GL_LIGHT0);
//
//    /*enable depth comparisons and update the depth buffer*/
//    glEnable(GL_DEPTH_TEST);
//    /*设置特殊效果*/
//    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
//    glHint(GL_LINE_SMOOTH_HINT,GL_DONT_CARE);
//    glEnable(GL_BLEND);
//
//    glEnable(GL_AUTO_NORMAL);
//    glEnable(GL_NORMALIZE);
//    glFrontFace(GL_CW);
//    glShadeModel(GL_SMOOTH);
//    glEnable(GL_LINE_SMOOTH);
//
//}
//
//void myDisplay(void)
//{
//    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
//    glColor3f(0.4,0.4,0.4);
//    glTranslatef(0.0,-1.0,0.0);
//    glRotatef(50.0,1.0,0.0,0.0);
//    glPushMatrix();
//    /*绘制曲面*/
//    glEnable(GL_MAP2_VERTEX_3);
//    glMap2f(GL_MAP2_VERTEX_3,0,1,3,5,0,1,15,5,&ctrlpoints[0][0][0]);
//    glMapGrid2f(10.0,0.0,1.0,10.0,0.0,1.0);
//    glEvalMesh2(GL_FILL,0,10.0,0,10.0);
//    glPopMatrix();
//    glutSwapBuffers();
//}
//
//void myReshape(GLsizei w,GLsizei h) {
//    glViewport(0, 0, w, h);
//    glMatrixMode(GL_PROJECTION);
//    glLoadIdentity();
//    gluPerspective(60.0, (GLfloat) w / (GLfloat) h, 1.0, 100.0);
//    glMatrixMode(GL_MODELVIEW);
//    glLoadIdentity();
//    glTranslatef(0.0, 0.0, -5.0);
//}
//
//int plot() {
//    /*初始化*/
//
//
//    /*创建窗口*/
//    glutCreateWindow("lighted Bezier surface");
//
//    /*绘制与显示*/
//    myInit();
//    glutReshapeFunc(myReshape);
//    glutDisplayFunc(myDisplay);
//
//    /*进入GLUT事件处理循环*/
//    glutMainLoop();
//    return (0);
//}
