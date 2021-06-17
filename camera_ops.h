#include <klein/klein.hpp>
#include "cayley.h"
#include "outer_exp.h"

#include "ceres/ceres.h"


using ceres::NumericDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;


#define CERES_NPOINTS_2X 10



void generate_internal_matrix(float params[5], float matrix[3][3]){
    /* 
    Turns a list of 5 parameters into an intrinsic matrix
    params[6] = {fx, fy, s, cx, cy} 
    */
    matrix[0][0] = params[0];
    matrix[1][1] = params[1];
    matrix[2][2] = 1.0f;
    matrix[0][1] = params[2];
    matrix[0][2] = params[3];
    matrix[0][3] = params[4];
    matrix[1][0] = 0.0f;
    matrix[2][0] = 0.0f;
    matrix[2][1] = 0.0f;
}


void apply_tangential_distortion(double &x_dist, double &y_dist, 
                    double x_correct, double y_correct,
                    double p1, double p2){
    /* 
    Applies tangential distortion to a point
    */
   double r2 = x_correct*x_correct + y_correct*y_correct;
   x_dist = x_correct + 2*p1*x_correct*y_correct + p2*(r2 + 2*x_correct*x_correct);
   y_dist = y_correct + 2*p2*x_correct*y_correct + p1*(r2 + 2*y_correct*y_correct);
}


void remove_tangential_distortion(double x_dist, double y_dist, 
                    double &x_correct, double &y_correct,
                    double p1, double p2, 
                    std::uint8_t max_iterations=10){
    /* 
    Removes tangential distortion from a point.
    We are just going to attempt the same iterative inverse thing we did for the
    radial distortion estimation.
    */
   x_correct = x_dist;
   y_correct = y_dist;
   double divisor_x = 1.0;
   double divisor_y = 1.0;
   double r2 = 1.0;
   for (std::uint8_t i=0; i<max_iterations; i++){
       r2 = x_correct*x_correct + y_correct*y_correct;
       divisor_x = 1.0 + 2*p1*y_correct + (p2*r2/x_correct) + 2*p2*x_correct;
       divisor_y = 1.0 + 2*p2*x_correct + (p1*r2/x_correct) + 2*p1*y_correct;
       x_correct = x_dist/divisor_x;
       y_correct = y_dist/divisor_y;
   }
}


void apply_radial_distortion(double &x_dist, double &y_dist, 
                    double x_correct, double y_correct,
                    double k1, double k2, double k3=0.0){
    /* 
    Applies radial distortion to a point
    */
    double r2 = x_correct*x_correct + y_correct*y_correct;
    double r4 = r2*r2;
    double r6 = r4*r2;
    double polynomial = 1.0f + k1*r2 + k2*r4 + k3*r6;
    x_dist = x_correct*polynomial;
    y_dist = y_correct*polynomial;
}


void remove_radial_distortion_iterative(double x_dist, double y_dist, 
                    double &x_correct, double &y_correct,
                    double k1, double k2, double k3, 
                    std::uint8_t max_iterations=10){
    /* 
    Removes radial distortion from a point by a simple iterative method.
    */
    x_correct = x_dist;
    y_correct = y_dist;
    double r2 = 1.0;
    double r4 = 1.0;
    double r6 = 1.0;
    double polynomial = 1.0;
    for (std::uint8_t i=0; i<max_iterations; i++){
        r2 = x_correct*x_correct + y_correct*y_correct;
        r4 = r2*r2;
        r6 = r4*r2;
        polynomial = 1.0 + k1*r2 + k2*r4 + k3*r6;
        x_correct = x_dist/polynomial;
        y_correct = y_dist/polynomial;
    }
}


void remove_radial_distortion(double x_dist, double y_dist, 
                    double &x_correct, double &y_correct,
                    double k1, double k2, double k3){
    /* 
    Removes radial distortion from a point using the method presented in:
    Pierre Drap and Julien Lefevre: An Exact Formula for Calculating Inverse Radial Lens Distortions
    */
    double r2 = x_dist*x_dist + y_dist*y_dist;
    double r4 = r2*r2;
    double r6 = r4*r2;

    double b1 = -k1;
    double b2 = 3*k1*k1 - k2;
    double b3 = -12*k1*k1*k1 + 8*k1*k2 - k3;

    double polynomial = 1.0 + b1*r2 + b2*r4 + b3*r6;
    x_correct = x_dist*polynomial;
    y_correct = y_dist*polynomial;
}


void project_to_camera(kln::motor &R, 
                        std::vector<std::shared_ptr<kln::point>> &points, 
                        std::vector<std::shared_ptr<kln::point>> &output_points){

    output_points.reserve(points.size()); 

    kln::point remapped_point;
    kln::point origin = kln::origin();
    kln::plane standard_plane = {0.0f, 0.0f, -1.0f, 1.0f};
    for (std::size_t i=0; i< points.size(); i++){ 
        // Map to the camera coordinate system
        remapped_point = (~R)(*points[i]);
        // Intersect
        kln::point resulting_point = (remapped_point & origin) ^ standard_plane;
        // Store
        auto pnt = std::make_shared<kln::point> (resulting_point);
        output_points.push_back(std::move(pnt));
    }
}


float reprojection_error(kln::motor &R, 
                        std::vector<std::shared_ptr<kln::point>> &points, 
                        std::vector<std::shared_ptr<kln::point>> &camera_points){

    // Project the points to the standard camera plane
    std::vector<std::shared_ptr<kln::point>> output_points;
    output_points.reserve(points.size()); 
    project_to_camera(R, points, output_points);

    // Map the points to the pixel positions (intrinsics)
    // Actually just assume the camera points are already expressed in the form we want

    // Return the error with camera points
    float total_error = 0.0f;
    kln::point diff_pnt;
    for (std::size_t i=0; i< camera_points.size(); i++){ 
        diff_pnt = ((*camera_points[i]) - output_points[i]->normalized());
        auto x = diff_pnt.x();
        auto y = diff_pnt.y();
        total_error += std::sqrt(x*x + y*y);
    }
    return total_error;
}


float reprojection_residuals(kln::motor &R, 
                        std::vector<std::shared_ptr<kln::point>> &points, 
                        std::vector<std::shared_ptr<kln::point>> &camera_points,
                        double* residual){

    // Project the points to the standard camera plane
    std::vector<std::shared_ptr<kln::point>> output_points;
    output_points.reserve(points.size()); 
    project_to_camera(R, points, output_points);

    // Map the points to the pixel positions (intrinsics)
    // Actually just assume the camera points are already expressed in the form we want

    // Return the error with camera points
    float total_error = 0.0f;
    kln::point diff_pnt;
    for (std::size_t i=0; i< camera_points.size(); i++){ 
        diff_pnt = ((*camera_points[i]) - output_points[i]->normalized());
        residual[2*i] = (double) diff_pnt.x();
        residual[2*i + 1] = (double) diff_pnt.y();
    }
    return total_error;
}


void find_camera(kln::line initial_biv,
                std::vector<std::shared_ptr<kln::point>>& points, 
                std::vector<std::shared_ptr<kln::point>>& camera_points){


  // The variable to solve for with its initial value. It will be
  // mutated in place by the solver.
  double x[6];
  x[0] = initial_biv.e01();
  x[1] = initial_biv.e02();
  x[2] = initial_biv.e03();
  x[3] = initial_biv.e23();
  x[4] = initial_biv.e31();
  x[5] = initial_biv.e12();
  const double initial_x[6] = {x[0], x[1], x[2], x[3], x[4], x[5]};


  // Build the problem.
  Problem problem;


    // Set up the cost function
    struct CostFunctor {
        std::vector<std::shared_ptr<kln::point>>& points;
        std::vector<std::shared_ptr<kln::point>>& camera_points;
        bool operator()(const double* const parameters, double* residuals) const {
            // Set up a camera
            kln::line biv_est = {parameters[0], parameters[1], parameters[2], 
                                parameters[3], parameters[4], parameters[5]};
            // kln::motor R_est = kln::exp(biv_est);
            kln::motor R_est = outer_exp(biv_est);
            // kln::motor R_est = cayley(biv_est);

            // Calculate the reprojection error
            residuals[0] = reprojection_error(R_est, this->points, this->camera_points);
            return true;
        }
    };


    // Set up the residual function
    struct NumericDiffCostFunctor {
        std::vector<std::shared_ptr<kln::point>>& points;
        std::vector<std::shared_ptr<kln::point>>& camera_points;
        bool operator()(const double* const parameters, double* residuals) const {
            // Set up a camera
            kln::line biv_est = {parameters[0], parameters[1], parameters[2], 
                                parameters[3], parameters[4], parameters[5]};
            // kln::motor R_est = kln::exp(biv_est);
            kln::motor R_est = outer_exp(biv_est);
            // kln::motor R_est = cayley(biv_est);

            // Calculate the reprojection error
            reprojection_residuals(R_est, this->points, this->camera_points, residuals);
            return true;
        }
    };

    // Make a cost function pointer that is then owned by the problem
    CostFunction* cost_function =
        new NumericDiffCostFunction<NumericDiffCostFunctor, ceres::CENTRAL, CERES_NPOINTS_2X, 6>(
            new NumericDiffCostFunctor{points, camera_points});
    
    problem.AddResidualBlock(cost_function, nullptr, &x[0]);



     // Run the solver!
    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(options, &problem, &summary);
    std::cout << summary.FullReport() << "\n";
    std::cout << "x[0] : " << initial_x[0] << " -> " << x[0] << "\n";
    std::cout << "x[1] : " << initial_x[1] << " -> " << x[1] << "\n";
    std::cout << "x[2] : " << initial_x[2] << " -> " << x[2] << "\n";
    std::cout << "x[3] : " << initial_x[3] << " -> " << x[3] << "\n";
    std::cout << "x[4] : " << initial_x[4] << " -> " << x[4] << "\n";
    std::cout << "x[5] : " << initial_x[5] << " -> " << x[5] << "\n";



    // // Calculate the cost before and after optimisation
    CostFunctor cost_tester = {points, camera_points};
    double op;
    cost_tester(initial_x, &op);
    std::cout << "pre cost " << op << std::endl;

    cost_tester(x, &op);
    std::cout << "post cost " << op << std::endl;


}




