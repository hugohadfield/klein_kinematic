
#include <memory>
#include <random>
#include <fstream>
#include <iostream>
#include "klein/klein.hpp"
#include "klein_ops.h"
#include "cayley.h"
#include "outer_exp.h"
#include "camera_ops.h"



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



int main(){

    const unsigned int npoints = CERES_NPOINTS_2X/2;

    std::ofstream output_stream;
    output_stream.open("200.txt");

    // Generate a load of points
    std::vector<std::shared_ptr<kln::point>> points;
    points.reserve(npoints);
    generate_random_points(points, npoints);

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
    kln::line biv_est{0.0f, 0.1f, 0.0f, 0.0f, 0.0f, 3.8f};
    kln::motor R_est = kln::exp(biv_est);

    // Calculate the reprojection error
    output = reprojection_error(R_est, points, camera_points);
    std::cout << output << std::endl;

    kln::line biv_est_outer = outer_log(R_est);
    kln::line biv_est_cayley = cayley(R_est);

    find_camera(biv_est_outer, points, camera_points); 

}

