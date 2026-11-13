cmake_minimum_required(VERSION 3.16)
project(cartesian_limit_adapter)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
add_compile_options(-Wall -Wextra -Wpedantic)

find_package(ament_cmake REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(moveit_core REQUIRED)
find_package(pluginlib REQUIRED)
find_package(rclcpp REQUIRED)

# The algorithm alone (Eigen only), unit-tested without MoveIt.
add_library(${PROJECT_NAME}_core SHARED src/toppra_cartesian.cpp)
target_include_directories(${PROJECT_NAME}_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include/${PROJECT_NAME}>)
target_link_libraries(${PROJECT_NAME}_core PUBLIC Eigen3::Eigen)
target_compile_features(${PROJECT_NAME}_core PUBLIC cxx_std_17)

# The MoveIt plugin.
add_library(${PROJECT_NAME} SHARED src/limit_cartesian_speed.cpp)
target_link_libraries(${PROJECT_NAME}
  ${PROJECT_NAME}_core
  moveit_core::moveit_planning_interface
  moveit_core::moveit_robot_state
  moveit_core::moveit_robot_trajectory
  pluginlib::pluginlib
  rclcpp::rclcpp)

pluginlib_export_plugin_description_file(moveit_core planning_response_adapters.xml)

install(TARGETS ${PROJECT_NAME}_core ${PROJECT_NAME}
  EXPORT export_${PROJECT_NAME}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin)
install(DIRECTORY include/ DESTINATION include/${PROJECT_NAME})

if(BUILD_TESTING)
  find_package(ament_cmake_gtest REQUIRED)
  ament_add_gtest(test_toppra_cartesian test/test_toppra_cartesian.cpp)
  target_link_libraries(test_toppra_cartesian ${PROJECT_NAME}_core)
endif()

ament_export_targets(export_${PROJECT_NAME} HAS_LIBRARY_TARGET)
ament_export_dependencies(Eigen3 moveit_core pluginlib rclcpp)
ament_package()
