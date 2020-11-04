#include "CVLib.h"
#include "utils.h"

#ifdef USE_OPENCV

cv::RotatedRect getErrorEllipse(const double& chisquare_val,
                                const cv::Point2f& means,
                                const cv::Mat& covmat){

    //Get the eigenvalues and eigenvectors
    cv::Mat eigenvalues, eigenvectors;
    cv::eigen(covmat, eigenvalues, eigenvectors);

    //Calculate the angle between the largest eigenvector and the x-axis
    double angle = atan2(eigenvectors.at<double>(0,1), eigenvectors.at<double>(0,0));

    //Shift the angle to the [0, 2pi] interval instead of [-pi, pi]
    if(angle < 0)
        angle += 6.28318530718;

    //Conver to degrees instead of radians
    angle = 180*angle/3.14159265359;

    //Calculate the size of the minor and major axes
    double halfmajoraxissize= chisquare_val*sqrt(eigenvalues.at<double>(0));
    double halfminoraxissize= chisquare_val*sqrt(eigenvalues.at<double>(1));

    //Return the oriented ellipse
    //The -angle is used because OpenCV defines the angle clockwise instead of anti-clockwise
    return cv::RotatedRect(means, cv::Size2f(halfmajoraxissize, halfminoraxissize), angle);

}

void plotDataPoint(cv::Mat& img,
                   const int& x,
                   const int& y,
                   const int& r,
                   const cv::Scalar& color)
{
    cv::circle(img, cv::Point(x, y), r, color, -1);
    return;
}

void plotCov4d(cv::Mat& img,
               double& x_min,
               double& x_max,
               double& y_min,
               double& y_max,
               const Eigen::Matrix4d& cov,
               const std::vector<Eigen::Vector4d>& data,
               const Eigen::Vector4d& mean,
               const int& a,
               const int& b,
               const int& ellipse_size,
               const cv::Scalar& ellipse_color,
               const int& point_size,
               const cv::Scalar& point_color,
               const double& size,
               const int& linestyle,
               const Eigen::Vector2d& bias)
{
    cv::Mat cov_cv = (cv::Mat_<double>(2, 2) << cov(a - 1, a - 1), cov(a - 1, b - 1), cov(b - 1, a - 1), cov(b - 1, b - 1));
    cv::Point2d mean_cv(mean(a - 1), mean(b - 1));
    double scale = 5e6;


    std::vector<double> demeaned_x(data.size());
    std::vector<double> demeaned_y(data.size());
    for(int i = 0; i < data.size(); ++i)
    {
        demeaned_x[i] = data[i](a - 1) - mean_cv.x;
        demeaned_y[i] = data[i](b - 1) - mean_cv.y;
    }
    std::sort(demeaned_x.begin(), demeaned_x.end());
    std::sort(demeaned_y.begin(), demeaned_y.end());
    x_min = (demeaned_x[0]);
    x_max = (demeaned_x[data.size() - 1]);
    y_min = (demeaned_y[0]);
    y_max = (demeaned_y[data.size() - 1]);
    scale = 1.0 / MAX4(fabs(x_min), fabs(x_max), fabs(y_min), fabs(y_max));
    scale = scale * size * 0.9;

    double chiqaured = 6.636 * scale;
    cv::RotatedRect ellipse = getErrorEllipse(chiqaured, cv::Point2f(bias.x(), bias.y()), cov_cv);
    cv::ellipse(img, ellipse, ellipse_color, 2, linestyle);

    for(int i = 0; i < data.size(); ++i)
    {
        plotDataPoint(img,
                      (int)((data[i](a - 1) - mean_cv.x) * scale + bias.x()),
                      (int)((data[i](b - 1) - mean_cv.y) * scale + bias.y()),
                      point_size,
                      point_color);
    }
}

void plotQuatCov(cv::Mat& img,
                 const Eigen::Matrix4d& cov1,
                 const Eigen::Matrix4d& cov2,
                 const std::vector<Eigen::Vector4d>& qs,
                 const Eigen::Vector4d& mean_q)
{
    int num_data = qs.size();
    int margin = img.rows / 50;
    double size = img.rows / 6.0 - margin * 2;
    double spacing = 5.0;
    double textlen = 200.0;

    Eigen::Vector2d bias;
    double x_min, x_max, y_min, y_max;
    std::string str_min, str_max;
    char str[512];

    for(int i = 0; i < 4; ++i)
    {
        for(int j = i + 1; j < 4; ++j)
        {
            bias(0) = (2.0 * j - 1.0) * img.rows / 6.0;
            bias(1) = (2.0 * i + 1.0) * img.cols / 6.0;
            plotCov4d(img, x_min, x_max, y_min, y_max,
                      cov1, qs, mean_q, i + 1, j + 1,
                      12, cv::Scalar(0, 0, 0) ,
                      2, cv::Scalar(1, 1, 1) * 120, size, cv::LINE_8, bias);

            double rect_xx, rect_xy;
            double rect_yx, rect_yy;
            rect_xx = (j - 1) * img.cols / 3.0;
            rect_yx = j * img.cols / 3.0;
            rect_xy = i * img.cols / 3.0;
            rect_yy = (i + 1) * img.cols / 3.0;

            cv::rectangle( img,
                           cv::Point( rect_xx + margin, rect_xy + margin ),
                           cv::Point( rect_yx - margin, rect_yy - margin),
                           cv::Scalar( 0, 0, 0 ), 2);

            std::sprintf(str, "q%d", i);
            str_min = std::string(str);
            putText(img, str_min, cv::Point(rect_xx, (rect_xy + rect_yy) / 2.0),
                    cv::FONT_HERSHEY_COMPLEX, 1, cv::Scalar(0, 0, 0), 2);
            std::sprintf(str, "q%d", j);
            str_min = std::string(str);
            putText(img, str_min, cv::Point((rect_yx + rect_xx) / 2.0, rect_yy - margin + 5 * spacing),
                    cv::FONT_HERSHEY_COMPLEX, 1, cv::Scalar(0, 0, 0), 2);


            std::sprintf(str, "%.3e", y_min);
            str_min = std::string(str);
            std::sprintf(str, "%.3e", y_max);
            str_max = std::string(str);
            putText(img, str_min, cv::Point(rect_xx + margin + 1 * spacing, rect_xy + margin + 5 * spacing),
                    cv::FONT_HERSHEY_COMPLEX, 1, cv::Scalar(0, 0, 0), 2);

            putText(img, str_max, cv::Point(rect_xx + margin + 1 * spacing, rect_yy - margin - 2 * spacing),
                    cv::FONT_HERSHEY_COMPLEX, 1, cv::Scalar(0, 0, 0), 2);

            std::sprintf(str, "%.3e", x_min);
            str_min = std::string(str);
            std::sprintf(str, "%.3e", x_max);
            str_max = std::string(str);
            putText(img, str_min, cv::Point(rect_xx + spacing, rect_yy - margin + 5 * spacing),
                    cv::FONT_HERSHEY_COMPLEX, 1, cv::Scalar(0, 0, 0), 2);

            putText(img, str_max, cv::Point(rect_yx - textlen, rect_yy - margin + 5 * spacing),
                    cv::FONT_HERSHEY_COMPLEX, 1, cv::Scalar(0, 0, 0), 2);
        }
    }

    margin = img.rows / 200;
    plotDataPoint(img,
                  0 * img.cols / 3 + 5 * margin, 2 * img.rows / 3 + 3 * margin,
                  8, cv::Scalar(1, 1, 1) * 120);

    putText(img, "Data Points", cv::Point(0 * img.cols / 3 + 10 * margin, 2 * img.rows / 3 + 5 * margin),
            cv::FONT_HERSHEY_COMPLEX, 1, cv::Scalar(0, 0, 0), 2);

    cv::line( img,
              cv::Point(0 * img.cols / 3 + 1 * margin, 2 * img.rows / 3 + 12 * margin),
              cv::Point(0 * img.cols / 3 + 8 * margin, 2 * img.rows / 3 + 12 * margin),
              cv::Scalar( 0, 0, 0 ), 2);

    putText(img, "Stat Quaternion Covariance", cv::Point(0 * img.cols / 3 + 10 * margin, 2 * img.rows / 3 + 14 * margin),
            cv::FONT_HERSHEY_COMPLEX, 1, cv::Scalar(0, 0, 0), 2);
}

#endif
