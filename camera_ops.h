
#include <klein/klein.hpp>

#include "ceres/ceres.h"


using ceres::NumericDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;


#define CERES_NPOINTS_2X 10


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

    // // Check that the cost function evaluates the same way
    // std::vector<const double*> parameters(2);
    // parameters[0] = &initial_x[0];
    // parameters[1] = &initial_x[1];
    // cost_function->Evaluate(parameters.data(), &op, NULL);
    // std::cout << "init cost " << op << std::endl;

    // parameters[0] = &x[0];
    // parameters[1] = &x[1];
    // cost_function->Evaluate(parameters.data(), &op, NULL);
    // std::cout << "end cost " << op << std::endl;


}





