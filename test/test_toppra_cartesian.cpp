// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Ignacio Gavier and LaurusAg
//
#include <gtest/gtest.h>

#include <cmath>

#include "cartesian_limit_adapter/toppra_cartesian.hpp"

using cartesian_limit_adapter::CubicPath;
using cartesian_limit_adapter::LinkDerivatives;
using cartesian_limit_adapter::Limits;
using cartesian_limit_adapter::parameterize;

namespace
{
const LinkDerivatives kNoLink = [](const Eigen::VectorXd&, const Eigen::VectorXd&, const Eigen::VectorXd&,
                                   Eigen::Vector3d& p1, Eigen::Vector3d& p2) {
  p1.setZero();
  p2.setZero();
};

// Planar 2R arm, links of 0.5 m: tip position in the plane z = 0.
Eigen::Vector3d tip(const Eigen::VectorXd& q)
{
  return { 0.5 * std::cos(q[0]) + 0.5 * std::cos(q[0] + q[1]),
           0.5 * std::sin(q[0]) + 0.5 * std::sin(q[0] + q[1]), 0.0 };
}

Eigen::Matrix<double, 3, 2> jac(const Eigen::VectorXd& q)
{
  Eigen::Matrix<double, 3, 2> J = Eigen::Matrix<double, 3, 2>::Zero();
  J(0, 0) = -0.5 * std::sin(q[0]) - 0.5 * std::sin(q[0] + q[1]);
  J(0, 1) = -0.5 * std::sin(q[0] + q[1]);
  J(1, 0) = 0.5 * std::cos(q[0]) + 0.5 * std::cos(q[0] + q[1]);
  J(1, 1) = 0.5 * std::cos(q[0] + q[1]);
  return J;
}

const LinkDerivatives kArm = [](const Eigen::VectorXd& q, const Eigen::VectorXd& d1,
                                const Eigen::VectorXd& d2, Eigen::Vector3d& p1, Eigen::Vector3d& p2) {
  const double h = 1e-6;
  p1 = jac(q) * d1;
  p2 = jac(q) * d2 + (jac(q + h * d1) - jac(q - h * d1)) / (2 * h) * d1;
};
}  // namespace

TEST(CubicPath, DropsDuplicatesAndInterpolates)
{
  Eigen::MatrixXd wp(4, 2);
  wp << 0, 0, 0, 0, 1, 0, 1, 1;
  const CubicPath path(wp);
  ASSERT_TRUE(path.valid());
  EXPECT_NEAR(path.length(), 2.0, 1e-12);
  EXPECT_TRUE(path.q(1.0).isApprox(Eigen::Vector2d(1, 0), 1e-12));
  EXPECT_TRUE(path.q(2.0).isApprox(Eigen::Vector2d(1, 1), 1e-12));
}

TEST(CubicPath, ZeroLengthIsInvalid)
{
  Eigen::MatrixXd wp(3, 2);
  wp << 0.3, 0.1, 0.3, 0.1, 0.3, 0.1;
  EXPECT_FALSE(CubicPath(wp).valid());
}

TEST(Toppra, TrapezoidMatchesClosedForm)
{
  // One joint, 1 rad at v = 1 rad/s, a = 2 rad/s2: 0.5 s up, 0.5 s cruise, 0.5 s down.
  Eigen::MatrixXd wp(3, 1);
  wp << 0.0, 0.5, 1.0;
  Limits lim;
  lim.max_velocity = Eigen::VectorXd::Constant(1, 1.0);
  lim.max_acceleration = Eigen::VectorXd::Constant(1, 2.0);
  const auto r = parameterize(CubicPath(wp), kNoLink, lim);
  ASSERT_TRUE(r.ok) << r.message;
  EXPECT_NEAR(r.duration, 1.5, 0.01);
}

TEST(Toppra, RespectsCartesianAndJointLimits)
{
  // Curved path of a 2R arm.
  Eigen::MatrixXd wp(5, 2);
  wp << 0.0, 0.3, 0.6, 0.9, 1.4, 1.2, 2.0, 0.8, 2.4, 0.2;
  Limits lim;
  lim.max_velocity = Eigen::Vector2d(3.0, 3.0);
  lim.max_acceleration = Eigen::Vector2d(10.0, 10.0);
  lim.max_linear_speed = 0.5;
  lim.max_linear_acceleration = 1.0;
  const CubicPath path(wp);
  const auto r = parameterize(path, kArm, lim);
  ASSERT_TRUE(r.ok) << r.message;

  std::vector<double> t;
  std::vector<Eigen::VectorXd> q, dq, ddq;
  r.sample(path, 0.002, t, q, dq, ddq);
  double v_max = 0, a_max = 0, jv = 0, ja = 0;
  for (std::size_t k = 0; k < t.size(); ++k)
  {
    Eigen::Vector3d v, a;
    kArm(q[k], dq[k], ddq[k], v, a);  // with time derivatives: v = J qdot, a = J qddot + Jdot qdot
    v_max = std::max(v_max, v.norm());
    a_max = std::max(a_max, a.norm());
    jv = std::max(jv, dq[k].cwiseAbs().cwiseQuotient(lim.max_velocity).maxCoeff());
    ja = std::max(ja, ddq[k].cwiseAbs().cwiseQuotient(lim.max_acceleration).maxCoeff());
  }
  EXPECT_LE(v_max, 0.5 * 1.01);
  EXPECT_LE(a_max, 1.0 * 1.01);
  EXPECT_LE(jv, 1.01);
  EXPECT_LE(ja, 1.01);
  EXPECT_GT(v_max, 0.5 * 0.95);  // and it is not needlessly slow
  EXPECT_NEAR((tip(q.back()) - tip(path.q(path.length()))).norm(), 0.0, 1e-9);
}

TEST(Toppra, LooseCartesianLimitChangesNothing)
{
  Eigen::MatrixXd wp(4, 2);
  wp << 0.0, 0.3, 0.6, 0.9, 1.4, 1.2, 2.0, 0.8;
  Limits joints;
  joints.max_velocity = Eigen::Vector2d(1.0, 1.0);
  joints.max_acceleration = Eigen::Vector2d(2.0, 2.0);
  Limits loose = joints;
  loose.max_linear_speed = 100.0;
  loose.max_linear_acceleration = 100.0;
  const CubicPath path(wp);
  const auto a = parameterize(path, kNoLink, joints);
  const auto b = parameterize(path, kArm, loose);
  ASSERT_TRUE(a.ok && b.ok);
  EXPECT_NEAR(a.duration, b.duration, 1e-3 * a.duration);
}
