
#include <klein/klein.hpp>


// void apply_internal(kln::mat3x4 &internal, 
//                     std::vector<kln::point> &points_in, 
//                     std::vector<float[2]> &points_out){
//     float tmp_point[4];
//     for (std::size_t i=0; i< points_in.length(); i++){ 
//         internal(points_in[i]).store(tmp_point);
//         points_out[i][0] = tmp_point[0]/tmp_point[2];
//         points_out[i][1] = tmp_point[1];
//     }
// }



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






