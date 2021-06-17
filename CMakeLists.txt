cmake_minimum_required(VERSION 3.15)

include(FetchContent)

# This tracks the latest commit on the master branch of the Klein
# repository. Instead of `origin/master`, you can specify a specific
# commit hash or tag.
FetchContent_Declare(
    klein
    GIT_REPOSITORY https://github.com/jeremyong/Klein.git
    GIT_TAG origin/master
)
FetchContent_MakeAvailable(klein)

# Now, you can use target_link_libraries(your_lib PUBLIC klein::klein)
# If you can target SSE4.1 (~97% market penetration), you can link against
# the target klein::klein_sse42 instead.


project(kinematic_klein LANGUAGES C CXX)

add_executable(cayley cayley.cpp)
target_link_libraries(cayley PRIVATE klein::klein_sse42)

# Stuff from klein repo
target_compile_options(cayley
    PRIVATE
    -fno-omit-frame-pointer
    -fsanitize=address
    -Wall
    -Wno-comment # Needed for doxygen
)
target_link_options(cayley PRIVATE -fno-omit-frame-pointer -fsanitize=address)


add_executable(outer_exp outer_exp.cpp)
target_link_libraries(outer_exp PRIVATE klein::klein_sse42)

# Stuff from klein repo
target_compile_options(outer_exp
    PRIVATE
    -fno-omit-frame-pointer
    -fsanitize=address
    -Wall
    -Wno-comment # Needed for doxygen
)
target_link_options(outer_exp PRIVATE -fno-omit-frame-pointer -fsanitize=address)

find_package(Ceres REQUIRED)

add_executable(generate_points generate_points.cpp)
target_link_libraries(generate_points klein::klein_sse42 Ceres::ceres)

target_compile_options(generate_points
    PRIVATE
    -fno-omit-frame-pointer
    -fsanitize=address
    -Wall
    -Wno-comment # Needed for doxygen
)
target_link_options(generate_points PRIVATE -fno-omit-frame-pointer -fsanitize=address)


add_executable(test_ops test_ops.cpp)
target_link_libraries(test_ops klein::klein_sse42 Ceres::ceres)

target_compile_options(test_ops
    PRIVATE
    -fno-omit-frame-pointer
    -fsanitize=address
    -Wall
    -Wno-comment # Needed for doxygen
)
target_link_options(test_ops PRIVATE -fno-omit-frame-pointer -fsanitize=address)
