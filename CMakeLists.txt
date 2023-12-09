cmake_minimum_required(VERSION 3.8)
project(ndt_2d_karto)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(Boost REQUIRED COMPONENTS thread)
find_package(Eigen3 REQUIRED)
find_package(ndt_2d REQUIRED)
find_package(pluginlib REQUIRED)
find_package(rclcpp REQUIRED)

include_directories(include)

add_library(scan_matcher_karto SHARED
  src/Karto.cpp
  src/Mapper.cpp
  src/scan_matcher_karto.cpp
)
target_link_libraries(scan_matcher_karto Eigen3::Eigen ${Boost_LIBRARIES})
ament_target_dependencies(scan_matcher_karto
  ndt_2d
  pluginlib
  rclcpp
)

if(BUILD_TESTING)
  find_package(ament_cmake_cpplint REQUIRED)
  ament_cpplint(FILTERS "-whitespace/braces" "-whitespace/newline")
endif()

install(
  TARGETS
    scan_matcher_karto
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION lib/${PROJECT_NAME}
)

pluginlib_export_plugin_description_file(ndt_2d plugins.xml)
ament_export_dependencies(rosidl_default_runtime)
ament_package()
