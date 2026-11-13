// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Ignacio Gavier and LaurusAg
//
// Time-optimal path parameterization (TOPP-RA, Pham & Pham 2018) with joint
// velocity/acceleration limits plus a bound on the linear speed and linear
// acceleration of one link.
//
// The path q(s), s in [0, L], is a natural cubic spline through the waypoints,
// parameterized by joint-space chord length. With x = sdot^2 and u = sddot:
//
//   joints   |q'_j| sqrt(x) <= v_j            |q'_j u + q''_j x| <= a_j
//   link     |p'| sqrt(x) <= v_lin            |p' u + p'' x| <= a_lin
//
// where p' = J q' and p'' = J q'' + Jdot q'. p' u + p'' x lies in the plane of
// {p', p''}, so its norm is bounded by a K-gon inscribed in the circle of radius
// a_lin in that plane (never exceeds it; gives up at most 1 - cos(pi/K)).
//
// Each interval enforces the acceleration rows at both ends (interpolation
// scheme). The grid is uniform plus the waypoints (subdivided), and intervals
// whose midpoints exceed a limit are split and solved again.
#pragma once

#include <Eigen/Core>
#include <functional>
#include <string>
#include <vector>

namespace cartesian_limit_adapter
{

struct Limits
{
  Eigen::VectorXd max_velocity;      // per joint
  Eigen::VectorXd max_acceleration;  // per joint
  double max_linear_speed = 0.0;         // <= 0: no bound
  double max_linear_acceleration = 0.0;  // <= 0: no bound
};

struct Options
{
  int grid_points = 100;         // uniform part of the grid
  int waypoint_subdivisions = 2; // points per waypoint interval
  int polygon_sides = 12;        // K for the norm bound
  int refine_rounds = 3;
  double refine_tolerance = 0.005;  // fraction of a limit
};

class CubicPath
{
public:
  // waypoints: one row per waypoint. Consecutive duplicates (< eps) are dropped.
  explicit CubicPath(const Eigen::MatrixXd& waypoints, double eps = 1e-6);

  double length() const { return knots_.back(); }
  const std::vector<double>& knots() const { return knots_; }
  bool valid() const { return valid_; }

  Eigen::VectorXd q(double s) const;
  Eigen::VectorXd dq(double s) const;
  Eigen::VectorXd ddq(double s) const;

private:
  std::size_t segment(double s) const;
  std::vector<double> knots_;
  Eigen::MatrixXd points_;  // knots x dof
  Eigen::MatrixXd second_;  // spline second derivatives at the knots
  bool valid_ = false;
};

// Linear velocity / acceleration of the limited link along the path:
// given q, q' and q'', returns p' = J q' and p'' = J q'' + Jdot q'.
using LinkDerivatives = std::function<void(const Eigen::VectorXd& q, const Eigen::VectorXd& dq,
                                           const Eigen::VectorXd& ddq, Eigen::Vector3d& p1,
                                           Eigen::Vector3d& p2)>;

struct Result
{
  bool ok = false;
  std::string message;
  std::vector<double> s;  // grid
  std::vector<double> x;  // sdot^2 on the grid
  std::vector<double> u;  // sddot per interval
  std::vector<double> t;  // time on the grid
  double duration = 0.0;
  int rounds = 0;

  // Samples every dt (the last sample lands exactly on the duration).
  void sample(const CubicPath& path, double dt, std::vector<double>& times,
              std::vector<Eigen::VectorXd>& pos, std::vector<Eigen::VectorXd>& vel,
              std::vector<Eigen::VectorXd>& acc) const;
};

Result parameterize(const CubicPath& path, const LinkDerivatives& link, const Limits& limits,
                    const Options& options = Options());

}  // namespace cartesian_limit_adapter
