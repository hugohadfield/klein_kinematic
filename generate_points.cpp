
#include <memory>
#include <random>
#include <fstream>
#include <iostream>
#include "klein/klein.hpp"
#include "klein_ops.h"
#include "cayley.h"
#include "outer_exp.h"
#include "camera_ops.h"


#include "ceres/ceres.h"
using ceres::NumericDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;


// Set up the random number generator
std::default_random_engine default_generator;
std::normal_distribution<float> coordinate_distribution(0.0, 1.0);


void random_point(float array[3], std::default_random_engine& generator=default_generator){
    array[0] = coordinate_distribution(generator);
    array[1] = coordinate_distribution(generator);
    array[2] = coordinate_distribution(generator);
}


void generate_random_points(std::vector<std::shared_ptr<kln::point>>& points, 
        unsigned int npoints, 
        std::default_random_engine& generator=default_generator){
    float array[3];
    for (std::size_t i=0; i < npoints; i++){
        random_point(array, generator);
        std::shared_ptr<kln::point> my_shared_pointer(new kln::point(array[0], array[1], array[2]));
        points.push_back(std::move(my_shared_pointer));
    }
}


void write_points(std::ofstream& file, std::vector<std::shared_ptr<kln::point>>& points){
    std::cout << points.size() << std::endl;
    for (std::size_t i=0; i < points.size(); i++){
        file << points[i]->x() <<','<< points[i]->y() <<','<< points[i]->z() << std::endl;
    }
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
    struct NumericDiffCostFunctor {
        std::vector<std::shared_ptr<kln::point>>& points;
        std::vector<std::shared_ptr<kln::point>>& camera_points;
        bool operator()(const double* const parameters, double* residuals) const {
            // Set up a camera
            kln::line biv_est = {parameters[0], parameters[1], parameters[2], 
                                parameters[3], parameters[4], parameters[5]};
            kln::motor R_est = kln::exp(biv_est);

            // Calculate the reprojection error
            residuals[0] = (double) reprojection_error(R_est, this->points, this->camera_points);
            return true;
        }
    };

    // Instantiate the struct to capture the local data
    NumericDiffCostFunctor instantiated = {points, camera_points};


    CostFunction* cost_function =
        new NumericDiffCostFunction<NumericDiffCostFunctor, ceres::CENTRAL, 1, 6>(
            &instantiated);
    
    problem.AddResidualBlock(cost_function, nullptr, &x[0]);

     // Run the solver!
    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    Solve(options, &problem, &summary);
    std::cout << summary.BriefReport() << "\n";
    std::cout << "x[0] : " << initial_x[0] << " -> " << x[0] << "\n";
    std::cout << "x[1] : " << initial_x[1] << " -> " << x[1] << "\n";
    std::cout << "x[2] : " << initial_x[2] << " -> " << x[2] << "\n";
    std::cout << "x[3] : " << initial_x[3] << " -> " << x[3] << "\n";
    std::cout << "x[4] : " << initial_x[4] << " -> " << x[4] << "\n";
    std::cout << "x[5] : " << initial_x[5] << " -> " << x[5] << "\n";

    // Calculate the cost before and after optimisation
    double op;
    instantiated(initial_x, &op);
    std::cout << "pre cost " << op << std::endl;

    instantiated(x, &op);
    std::cout << "post cost " << op << std::endl;


}






int main(){
    std::ofstream output_stream;
    output_stream.open("100.txt");

    // Generate a load of points
    std::vector<std::shared_ptr<kln::point>> points;
    points.reserve(100);
    generate_random_points(points, 100);

    std::cout << points.size() << std::endl;

    // Write the points to file
    write_points(output_stream, points);
    output_stream.close();



    // Set up a true camera
    kln::line biv_cam{0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 3.0f};
    kln::motor R = kln::exp(biv_cam);

    // Project the points into that camera
    std::vector<std::shared_ptr<kln::point>> camera_points;
    camera_points.reserve(points.size()); 
    project_to_camera(R, points, camera_points);
    for (std::size_t i=0; i<camera_points.size(); i++){
        camera_points[i]->normalize();
    }

    // Assert that there is no reprojection error with the true camera
    auto output = reprojection_error(R, points, camera_points);
    std::cout << output << std::endl;
    


    // Set up a false camera
    kln::line biv_est{0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 3.3f};
    kln::motor R_est = kln::exp(biv_est);

    // Calculate the reprojection error
    output = reprojection_error(R_est, points, camera_points);
    std::cout << output << std::endl;


    find_camera(biv_est, points, camera_points);

}

