// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Ignacio Gavier and LaurusAg
//
// MoveIt planning response adapter: retimes the planned path with TOPP-RA so one
// link stays under a linear speed and a linear acceleration, on top of the joint
// velocity/acceleration limits of the robot model (URDF + joint_limits.yaml, the
// same ones TOTG uses).
//
// Put it after AddTimeOptimalParameterization and before ValidateSolution, so the
// trajectory that is validated and executed is the limited one.
//
// Parameters, under <pipeline>.cartesian_limit (read on every request, so
// `ros2 param set` takes effect at once):
//   enabled                  bool    true
//   link                     string  ""    (empty: last link of the group)
//   max_linear_speed         double  0.0   m/s  (<= 0: no bound)
//   max_linear_acceleration  double  0.0   m/s2 (<= 0: no bound)
//   resample_dt              double  0.02  s
//   grid_points, waypoint_subdivisions, polygon_sides, refine_rounds, refine_tolerance
//
// A request with max_cartesian_speed > 0 overrides max_linear_speed, and its
// cartesian_speed_limited_link (if set) overrides link.
#include <chrono>
#include <memory>
#include <string>

#include <moveit/planning_interface/planning_interface.hpp>
#include <moveit/planning_interface/planning_response_adapter.hpp>
#include <moveit/planning_scene/planning_scene.hpp>
#include <moveit/robot_state/robot_state.hpp>
#include <moveit/robot_trajectory/robot_trajectory.hpp>
#include <pluginlib/class_list_macros.hpp>
#include <rclcpp/rclcpp.hpp>

#include "cartesian_limit_adapter/toppra_cartesian.hpp"

namespace cartesian_limit_adapter
{
class LimitCartesianSpeed : public planning_interface::PlanningResponseAdapter
{
public:
  void initialize(const rclcpp::Node::SharedPtr& node, const std::string& parameter_namespace) override
  {
    node_ = node;
    logger_ = node->get_logger().get_child("limit_cartesian_speed");  // also goes to /rosout
    prefix_ = (parameter_namespace.empty() ? std::string() : parameter_namespace + ".") + "cartesian_limit.";
    declare("enabled", true);
    declare("link", std::string());
    declare("max_linear_speed", 0.0);
    declare("max_linear_acceleration", 0.0);
    declare("resample_dt", 0.02);
    declare("grid_points", 100);
    declare("waypoint_subdivisions", 2);
    declare("polygon_sides", 12);
    declare("refine_rounds", 3);
    declare("refine_tolerance", 0.005);
  }

  [[nodiscard]] std::string getDescription() const override
  {
    return "Limit Cartesian speed (TOPP-RA)";
  }

  void adapt(const planning_scene::PlanningSceneConstPtr& planning_scene,
             const planning_interface::MotionPlanRequest& req,
             planning_interface::MotionPlanResponse& res) const override
  {
    if (!res.error_code || !res.trajectory || res.trajectory->getWayPointCount() < 2 ||
        !get<bool>("enabled"))
      return;
    const auto started = std::chrono::steady_clock::now();
    const moveit::core::RobotModelConstPtr& model = planning_scene->getRobotModel();
    const moveit::core::JointModelGroup* group = res.trajectory->getGroup();
    if (group == nullptr)
      group = model->getJointModelGroup(req.group_name);
    if (group == nullptr)
      return fail(res, "no hay grupo para '" + req.group_name + "'");

    std::string link_name = get<std::string>("link");
    if (req.max_cartesian_speed > 0.0 && !req.cartesian_speed_limited_link.empty())
      link_name = req.cartesian_speed_limited_link;
    const moveit::core::LinkModel* link =
        link_name.empty() ? group->getLinkModels().back() : model->getLinkModel(link_name);
    if (link == nullptr)
      return fail(res, "no existe el link '" + link_name + "'");

    Limits limits;
    limits.max_linear_speed =
        req.max_cartesian_speed > 0.0 ? req.max_cartesian_speed : get<double>("max_linear_speed");
    limits.max_linear_acceleration = get<double>("max_linear_acceleration");
    if (limits.max_linear_speed <= 0.0 && limits.max_linear_acceleration <= 0.0)
      return;
    if (!jointLimits(group, req, limits))
      return fail(res, "el grupo tiene joints sin límites de velocidad");

    // The path: the joint positions of the planned waypoints.
    const std::size_t n_wp = res.trajectory->getWayPointCount();
    const std::size_t dof = group->getVariableCount();
    Eigen::MatrixXd wp(n_wp, dof);
    Eigen::VectorXd q;
    for (std::size_t k = 0; k < n_wp; ++k)
    {
      res.trajectory->getWayPoint(k).copyJointGroupPositions(group, q);
      wp.row(k) = q.transpose();
    }
    const CubicPath path(wp);
    if (!path.valid())
      return;  // already at the goal: nothing to limit

    moveit::core::RobotState state(res.trajectory->getFirstWayPoint());
    Eigen::MatrixXd J, Jp, Jm;
    auto jacobian = [&](const Eigen::VectorXd& qq, Eigen::MatrixXd& out) {
      state.setJointGroupPositions(group, qq);
      state.updateLinkTransforms();
      state.getJacobian(group, link, Eigen::Vector3d::Zero(), out, false);
    };
    const LinkDerivatives derivatives = [&](const Eigen::VectorXd& qq, const Eigen::VectorXd& d1,
                                            const Eigen::VectorXd& d2, Eigen::Vector3d& p1,
                                            Eigen::Vector3d& p2) {
      constexpr double h = 1e-6;
      jacobian(qq, J);
      jacobian(qq + h * d1, Jp);
      jacobian(qq - h * d1, Jm);
      p1 = J.topRows<3>() * d1;
      p2 = J.topRows<3>() * d2 + ((Jp.topRows<3>() - Jm.topRows<3>()) / (2.0 * h)) * d1;
    };

    Options opt;
    opt.grid_points = std::max(2, get<int>("grid_points"));
    opt.waypoint_subdivisions = std::max(1, get<int>("waypoint_subdivisions"));
    opt.polygon_sides = std::max(4, get<int>("polygon_sides"));
    opt.refine_rounds = std::max(0, get<int>("refine_rounds"));
    opt.refine_tolerance = get<double>("refine_tolerance");
    const Result result = parameterize(path, derivatives, limits, opt);
    if (!result.ok)
      return fail(res, "TOPP-RA: " + result.message);

    std::vector<double> times;
    std::vector<Eigen::VectorXd> pos, vel, acc;
    result.sample(path, std::max(1e-3, get<double>("resample_dt")), times, pos, vel, acc);
    auto out = std::make_shared<robot_trajectory::RobotTrajectory>(model, group);
    moveit::core::RobotState point(res.trajectory->getFirstWayPoint());
    for (std::size_t k = 0; k < times.size(); ++k)
    {
      point.setJointGroupPositions(group, pos[k]);
      point.setJointGroupVelocities(group, vel[k]);
      point.setJointGroupAccelerations(group, acc[k]);
      point.update();
      out->addSuffixWayPoint(point, k == 0 ? 0.0 : times[k] - times[k - 1]);
    }
    const double before = res.trajectory->getDuration();
    res.trajectory = out;
    const double ms =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - started).count();
    RCLCPP_INFO(logger_, "%s <= %.2f m/s, %.2f m/s2: %.2fs -> %.2fs (%.1f ms, %zu puntos, %d rondas)",
                link->getName().c_str(), limits.max_linear_speed, limits.max_linear_acceleration,
                before, result.duration, ms, result.s.size(), result.rounds);
  }

private:
  template <typename T>
  void declare(const std::string& name, const T& value)
  {
    if (!node_->has_parameter(prefix_ + name))
      node_->declare_parameter(prefix_ + name, value);
  }

  template <typename T>
  T get(const std::string& name) const
  {
    return node_->get_parameter(prefix_ + name).get_value<T>();
  }

  // Joint velocity/acceleration bounds of the model (URDF + joint_limits.yaml),
  // scaled like TOTG does with the request scaling factors.
  static bool jointLimits(const moveit::core::JointModelGroup* group,
                          const planning_interface::MotionPlanRequest& req, Limits& limits)
  {
    const auto scale = [](double f) { return (f > 0.0 && f <= 1.0) ? f : 1.0; };
    const double sv = scale(req.max_velocity_scaling_factor);
    const double sa = scale(req.max_acceleration_scaling_factor);
    const std::vector<std::string>& names = group->getVariableNames();
    limits.max_velocity.resize(names.size());
    limits.max_acceleration.resize(names.size());
    const moveit::core::RobotModel& model = group->getParentModel();
    for (std::size_t j = 0; j < names.size(); ++j)
    {
      const moveit::core::VariableBounds& b = model.getVariableBounds(names[j]);
      if (!b.velocity_bounded_)
        return false;
      limits.max_velocity[j] = std::min(std::abs(b.max_velocity_), std::abs(b.min_velocity_)) * sv;
      // TOTG uses 1.0 for joints without an acceleration bound.
      const double a = b.acceleration_bounded_ ?
                           std::min(std::abs(b.max_acceleration_), std::abs(b.min_acceleration_)) :
                           1.0;
      limits.max_acceleration[j] = a * sa;
    }
    return true;
  }

  void fail(planning_interface::MotionPlanResponse& res, const std::string& why) const
  {
    RCLCPP_ERROR(logger_, "%s", why.c_str());
    res.error_code = moveit::core::MoveItErrorCode(moveit_msgs::msg::MoveItErrorCodes::FAILURE, why,
                                                   "LimitCartesianSpeed");
  }

  rclcpp::Node::SharedPtr node_;
  rclcpp::Logger logger_ = rclcpp::get_logger("limit_cartesian_speed");
  std::string prefix_;
};

}  // namespace cartesian_limit_adapter

PLUGINLIB_EXPORT_CLASS(cartesian_limit_adapter::LimitCartesianSpeed,
                       planning_interface::PlanningResponseAdapter)
