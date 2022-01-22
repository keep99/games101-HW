#include <chrono>
#include <iostream>
#include <opencv2/opencv.hpp>

std::vector<cv::Point2f> control_points;

void mouse_handler(int event, int x, int y, int flags, void *userdata) 
{
    if (event == cv::EVENT_LBUTTONDOWN && control_points.size() < 5) 
    {
        std::cout << "Left button of the mouse is clicked - position (" << x << ", "
        << y << ")" << '\n';
        control_points.emplace_back(x, y);
    }     
}

void naive_bezier(const std::vector<cv::Point2f> &points, cv::Mat &window) 
{
    auto &p_0 = points[0];
    auto &p_1 = points[1];
    auto &p_2 = points[2];
    auto &p_3 = points[3];
    auto &p_4 = points[4];  // 5

    for (double t = 0.0; t <= 1.0; t += 0.001) 
    {
        // auto point = std::pow(1 - t, 3) * p_0 + 3 * t * std::pow(1 - t, 2) * p_1 +
        //          3 * std::pow(t, 2) * (1 - t) * p_2 + std::pow(t, 3) * p_3;
        
        auto point = std::pow(1 - t, 4) * p_0 + 4 * t * std::pow(1 - t, 3) * p_1 + 
                     6 * std::pow(t, 2) * std::pow((1- t), 2) * p_2 +
                     4 * std::pow(t, 3) * (1 - t) * p_3 + 
                     std::pow(t, 4) * p_4;

        window.at<cv::Vec3b>(point.y, point.x)[2] = 255;
    }
}

cv::Point2f recursive_bezier(const std::vector<cv::Point2f> &control_points, float t) 
{
    // Iteration
    // // TODO: Implement de Casteljau's algorithm
    // auto &p_0 = control_points[0];
    // auto &p_1 = control_points[1];
    // auto &p_2 = control_points[2];
    // auto &p_3 = control_points[3];

    // auto p_01 = cv::Point2f(p_0.x + (p_1.x - p_0.x) * t, p_0.y + (p_1.y - p_0.y) * t);
    // auto p_11 = cv::Point2f(p_1.x + (p_2.x - p_1.x) * t, p_1.y + (p_2.y - p_1.y) * t);
    // auto p_21 = cv::Point2f(p_2.x + (p_3.x - p_2.x) * t, p_2.y + (p_3.y - p_2.y) * t);

    // auto p_02 = cv::Point2f(p_01.x + (p_11.x - p_01.x) * t, p_01.y + (p_11.y - p_01.y) * t);
    // auto p_12 = cv::Point2f(p_11.x + (p_21.x - p_11.x) * t, p_11.y + (p_21.y - p_11.y) * t);
    
    // // last
    // auto p_03 = cv::Point2f(p_02.x + (p_12.x - p_02.x) * t, p_02.y + (p_12.y - p_02.y) * t);
    // return p_03;

    // Recursion
    if (control_points.size() == 2)
    {
        auto p_0 = control_points[0];
        auto p_1 = control_points[1];
        return cv::Point2f(p_0.x + (p_1.x - p_0.x) * t, p_0.y + (p_1.y - p_0.y) * t);
    }

    std::vector<cv::Point2f> tmp;
    int size = control_points.size();

    for (int i = 0; i < size - 1; ++i)
    {
        tmp.push_back(cv::Point2f(control_points[i].x + (control_points[i + 1].x - control_points[i].x) * t,
                                  control_points[i].y + (control_points[i + 1].y - control_points[i].y) * t));
    }

    return recursive_bezier(tmp, t);
}

void bezier(const std::vector<cv::Point2f> &control_points, cv::Mat &window) 
{
    // TODO: Iterate through all t = 0 to t = 1 with small steps, and call de Casteljau's 
    // recursive Bezier algorithm.
    for (double t = 0.0; t <= 1.0; t += 0.001)
    {
        auto return_point = recursive_bezier(control_points, t);
        window.at<cv::Vec3b>(return_point.y, return_point.x)[1] = 255;

        // TODO
        // Anti-aliasing
        float x = return_point.x;
        float y = return_point.y;
        float dx = x - std::floor(x);
        float dy = y - std::floor(y);

        cv::Point2f nearest;
        cv::Point2f p00;
        cv::Point2f p01;
        cv::Point2f p10;
        cv::Point2f p11;

        p00 = cv::Point2f(std::floor(x), std::floor(y));
        p01 = cv::Point2f(std::floor(x), std::floor(y) + 1);
        p10 = cv::Point2f(std::floor(x) + 1, std::floor(y));
        p11 = cv::Point2f(std::floor(x) + 1, std::floor(y) + 1);

        if (dx <= 0.5f && dy <= 0.5f)
        {
            nearest = p00;
        }
        if (dx <= 0.5f && dy > 0.5f)
        {
            nearest = p01;
        }
        if (dx > 0.5f && dy <= 0.5f)
        {
            nearest = p10;
        }
        if (dx > 0.5f && dy > 0.5f)
        {
            nearest = p11;
        }

        std::vector<cv::Point2f> vec{p00, p01, p10, p11};
        float distance = std::sqrt(std::pow((return_point - nearest).x, 2) + 
                                   std::pow((return_point - nearest).y, 2));  // nearest 

        for (auto a : vec)
        {
            float d = std::sqrt(std::pow((a - return_point).x, 2) + std::pow((a - return_point).y, 2));
            float ratio = distance / d;  // more far, more less
            cv::Vec3d color = window.at<cv::Vec3b>(return_point.y, return_point.x);
            color[1] = std::max(color[1], double(255 * ratio));
            window.at<cv::Vec3b>(a.y, a.x) = color;
        }
    }
}

int main() 
{
    cv::Mat window = cv::Mat(700, 700, CV_8UC3, cv::Scalar(0));
    cv::cvtColor(window, window, cv::COLOR_BGR2RGB);
    cv::namedWindow("Bezier Curve", cv::WINDOW_AUTOSIZE);

    cv::setMouseCallback("Bezier Curve", mouse_handler, nullptr);

    int key = -1;
    while (key != 27) 
    {
        for (auto &point : control_points) 
        {
            cv::circle(window, point, 3, {255, 255, 255}, 3);
        }

        if (control_points.size() == 5) 
        {
            naive_bezier(control_points, window);
            bezier(control_points, window);

            cv::imshow("Bezier Curve", window);
            cv::imwrite("my_bezier_curve.png", window);
            key = cv::waitKey(0);

            return 0;
        }

        cv::imshow("Bezier Curve", window);
        key = cv::waitKey(20);
    }

return 0;
}
