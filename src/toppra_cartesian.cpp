// SPDX-License-Identifier: BSD-3-Clause
// Copyright (c) 2026, Ignacio Gavier and LaurusAg
//
#include "cartesian_limit_adapter/toppra_cartesian.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <limits>

namespace cartesian_limit_adapter
{
namespace
{
constexpr double kInf = std::numeric_limits<double>::infinity();
constexpr double kTiny = 1e-12;

// a_u * u + a_x * x <= b
struct Row
{
  double au, ax, b;
};

// In (x, y = x_next): y <= g - h x (upper) or y >= g - h x (lower).
struct Line
{
  double g, h;
};

struct Interval
{
  double x_low = 0.0, x_high = kInf;  // bounds that do not depend on the next set
  std::vector<Line> upper, lower;
};

struct PointData
{
  Eigen::VectorXd q1, q2;
  Eigen::Vector3d p1 = Eigen::Vector3d::Zero(), p2 = Eigen::Vector3d::Zero();
};
}  // namespace

// ----------------------------------------------------------------------------
// Natural cubic spline, one per joint, on chord-length knots.

CubicPath::CubicPath(const Eigen::MatrixXd& waypoints, double eps)
{
  std::vector<int> keep;
  for (int k = 0; k < waypoints.rows(); ++k)
  {
    if (keep.empty() || (waypoints.row(k) - waypoints.row(keep.back())).norm() > eps)
      keep.push_back(k);
  }
  if (!keep.empty() && keep.back() != waypoints.rows() - 1)
    keep.back() = static_cast<int>(waypoints.rows()) - 1;  // the goal is exact
  if (keep.size() < 2)
    return;
  const int n = static_cast<int>(keep.size());
  points_.resize(n, waypoints.cols());
  for (int k = 0; k < n; ++k)
    points_.row(k) = waypoints.row(keep[k]);
  knots_.assign(n, 0.0);
  for (int k = 1; k < n; ++k)
    knots_[k] = knots_[k - 1] + (points_.row(k) - points_.row(k - 1)).norm();

  // Tridiagonal system for the second derivatives (natural: M_0 = M_{n-1} = 0).
  second_ = Eigen::MatrixXd::Zero(n, points_.cols());
  if (n > 2)
  {
    const int m = n - 2;
    std::vector<double> sub(m), diag(m), sup(m);
    Eigen::MatrixXd rhs(m, points_.cols());
    for (int i = 1; i <= m; ++i)
    {
      const double h0 = knots_[i] - knots_[i - 1];
      const double h1 = knots_[i + 1] - knots_[i];
      sub[i - 1] = h0;
      diag[i - 1] = 2.0 * (h0 + h1);
      sup[i - 1] = h1;
      rhs.row(i - 1) = 6.0 * ((points_.row(i + 1) - points_.row(i)) / h1 -
                              (points_.row(i) - points_.row(i - 1)) / h0);
    }
    for (int i = 1; i < m; ++i)  // Thomas algorithm
    {
      const double w = sub[i] / diag[i - 1];
      diag[i] -= w * sup[i - 1];
      rhs.row(i) -= w * rhs.row(i - 1);
    }
    second_.row(m) = rhs.row(m - 1) / diag[m - 1];
    for (int i = m - 2; i >= 0; --i)
      second_.row(i + 1) = (rhs.row(i) - sup[i] * second_.row(i + 2)) / diag[i];
  }
  valid_ = knots_.back() > 0.0;
}

std::size_t CubicPath::segment(double s) const
{
  auto it = std::upper_bound(knots_.begin(), knots_.end(), s);
  std::size_t i = (it == knots_.begin()) ? 0 : static_cast<std::size_t>(it - knots_.begin()) - 1;
  return std::min(i, knots_.size() - 2);
}

Eigen::VectorXd CubicPath::q(double s) const
{
  s = std::clamp(s, 0.0, length());
  const std::size_t i = segment(s);
  const double h = knots_[i + 1] - knots_[i];
  const double a = (knots_[i + 1] - s) / h, b = (s - knots_[i]) / h;
  return (a * points_.row(i) + b * points_.row(i + 1) +
          ((a * a * a - a) * second_.row(i) + (b * b * b - b) * second_.row(i + 1)) * h * h / 6.0)
      .transpose();
}

Eigen::VectorXd CubicPath::dq(double s) const
{
  s = std::clamp(s, 0.0, length());
  const std::size_t i = segment(s);
  const double h = knots_[i + 1] - knots_[i];
  const double a = (knots_[i + 1] - s) / h, b = (s - knots_[i]) / h;
  return ((points_.row(i + 1) - points_.row(i)) / h +
          (-(3 * a * a - 1) * second_.row(i) + (3 * b * b - 1) * second_.row(i + 1)) * h / 6.0)
      .transpose();
}

Eigen::VectorXd CubicPath::ddq(double s) const
{
  s = std::clamp(s, 0.0, length());
  const std::size_t i = segment(s);
  const double h = knots_[i + 1] - knots_[i];
  const double a = (knots_[i + 1] - s) / h, b = (s - knots_[i]) / h;
  return (a * second_.row(i) + b * second_.row(i + 1)).transpose();
}

// ----------------------------------------------------------------------------

namespace
{
bool usesLink(const Limits& lim)
{
  return lim.max_linear_speed > 0.0 || lim.max_linear_acceleration > 0.0;
}

PointData evaluate(const CubicPath& path, const LinkDerivatives& link, const Limits& lim, double s)
{
  PointData d;
  d.q1 = path.dq(s);
  d.q2 = path.ddq(s);
  if (usesLink(lim))
    link(path.q(s), d.q1, d.q2, d.p1, d.p2);
  return d;
}

double maxX(const PointData& d, const Limits& lim)
{
  double x = 1e8;
  for (int j = 0; j < d.q1.size(); ++j)
  {
    const double a = std::abs(d.q1[j]);
    if (a > kTiny)
      x = std::min(x, std::pow(lim.max_velocity[j] / a, 2));
  }
  const double n = d.p1.norm();
  if (lim.max_linear_speed > 0.0 && n > kTiny)
    x = std::min(x, std::pow(lim.max_linear_speed / n, 2));
  return x;
}

void accelerationRows(const PointData& d, const Limits& lim, int sides, std::vector<Row>& rows)
{
  rows.clear();
  for (int j = 0; j < d.q1.size(); ++j)
  {
    rows.push_back({ d.q1[j], d.q2[j], lim.max_acceleration[j] });
    rows.push_back({ -d.q1[j], -d.q2[j], lim.max_acceleration[j] });
  }
  if (lim.max_linear_acceleration > 0.0)
  {
    Eigen::Vector3d e1 = d.p1.norm() > kTiny ? d.p1 : d.p2;
    if (e1.norm() > kTiny)
    {
      e1.normalize();
      Eigen::Vector3d r = d.p2 - d.p2.dot(e1) * e1;
      const bool planar = r.norm() > 1e-9 * std::max(1.0, d.p2.norm());
      const Eigen::Vector3d e2 = planar ? Eigen::Vector3d(r.normalized()) : Eigen::Vector3d::Zero();
      const double b = lim.max_linear_acceleration * std::cos(M_PI / sides);
      for (int k = 0; k < sides; ++k)
      {
        const double phi = 2.0 * M_PI * k / sides;
        const Eigen::Vector3d dir = std::cos(phi) * e1 + std::sin(phi) * e2;
        rows.push_back({ dir.dot(d.p1), dir.dot(d.p2), b });
      }
    }
  }
  for (Row& r : rows)  // normalize so tolerances mean the same for every row
  {
    const double n = std::hypot(r.au, r.ax);
    if (n > kTiny)
      r = { r.au / n, r.ax / n, r.b / n };
  }
}

// Interval between grid points i and i+1 in (x, y = x_{i+1}), with u = (y - x) / 2D.
bool buildInterval(const std::vector<Row>& a0, const std::vector<Row>& a1, double D, double xmax,
                   Interval& out)
{
  out = Interval();
  out.x_high = xmax;
  const double k2 = 1.0 / (2.0 * D);
  auto add = [&](double cx, double cy, double b) {
    const double n = std::hypot(cx, cy);
    if (n < kTiny)
      return b >= -1e-9;
    cx /= n, cy /= n, b /= n;
    if (cy > 1e-9)
      out.upper.push_back({ b / cy, cx / cy });
    else if (cy < -1e-9)
      out.lower.push_back({ b / cy, cx / cy });
    else if (cx > 1e-9)
      out.x_high = std::min(out.x_high, b / cx);
    else if (cx < -1e-9)
      out.x_low = std::max(out.x_low, b / cx);
    else if (b < -1e-9)
      return false;
    return true;
  };
  for (const Row& r : a0)  // at s_i: a_u (y - x) k2 + a_x x <= b
    if (!add(r.ax - r.au * k2, r.au * k2, r.b))
      return false;
  for (const Row& r : a1)  // at s_{i+1}: a_u (y - x) k2 + a_x y <= b
    if (!add(-r.au * k2, r.ax + r.au * k2, r.b))
      return false;
  // Fourier-Motzkin on y: every (upper, lower) pair bounds x.
  for (const Line& up : out.upper)
  {
    for (const Line& lo : out.lower)
    {
      const double dh = up.h - lo.h, dg = up.g - lo.g;  // dh x <= dg
      if (dh > kTiny)
        out.x_high = std::min(out.x_high, dg / dh);
      else if (dh < -kTiny)
        out.x_low = std::max(out.x_low, dg / dh);
      else if (dg < -1e-9)
        return false;
    }
  }
  return out.x_low <= out.x_high + 1e-9;
}

Result solve(const CubicPath& path, const LinkDerivatives& link, const Limits& lim,
             const std::vector<double>& s, const std::vector<PointData>& pts, int sides)
{
  Result res;
  const std::size_t G = s.size();
  std::vector<Interval> iv(G - 1);
  std::vector<Row> a0, a1;
  accelerationRows(pts[0], lim, sides, a0);
  for (std::size_t i = 0; i + 1 < G; ++i)
  {
    accelerationRows(pts[i + 1], lim, sides, a1);
    if (!buildInterval(a0, a1, s[i + 1] - s[i], maxX(pts[i], lim), iv[i]))
    {
      res.message = "intervalo " + std::to_string(i) + " sin región factible";
      return res;
    }
    a0.swap(a1);
  }

  // Backward pass: controllable sets [lo_i, hi_i], ending at rest.
  std::vector<double> lo(G, 0.0), hi(G, 0.0);
  for (std::size_t k = G - 1; k-- > 0;)
  {
    const Interval& I = iv[k];
    double xl = I.x_low, xh = I.x_high;
    for (const Line& l : I.upper)  // lo' <= g - h x
    {
      if (std::abs(l.h) <= kTiny)
      {
        if (l.g < lo[k + 1] - 1e-12)
          xh = -1.0;
      }
      else if (l.h > 0)
        xh = std::min(xh, (l.g - lo[k + 1]) / l.h);
      else
        xl = std::max(xl, (l.g - lo[k + 1]) / l.h);
    }
    for (const Line& l : I.lower)  // g - h x <= hi'
    {
      if (std::abs(l.h) <= kTiny)
      {
        if (l.g > hi[k + 1] + 1e-12)
          xh = -1.0;
      }
      else if (l.h < 0)
        xh = std::min(xh, (l.g - hi[k + 1]) / l.h);
      else
        xl = std::max(xl, (l.g - hi[k + 1]) / l.h);
    }
    if (xl > xh + 1e-10 * (1.0 + std::abs(xh)))
    {
      res.message = "etapa " + std::to_string(k) + " no controlable";
      return res;
    }
    lo[k] = std::max(xl, 0.0);
    hi[k] = std::max(xh, 0.0);
  }
  if (lo[0] > 1e-9)
  {
    res.message = "no se puede arrancar en reposo";
    return res;
  }

  // Forward pass: the largest feasible next x.
  res.x.assign(G, 0.0);
  for (std::size_t i = 0; i + 1 < G; ++i)
  {
    double yh = hi[i + 1], yl = lo[i + 1];
    for (const Line& l : iv[i].upper)
      yh = std::min(yh, l.g - l.h * res.x[i]);
    for (const Line& l : iv[i].lower)
      yl = std::max(yl, l.g - l.h * res.x[i]);
    res.x[i + 1] = std::max(yh >= yl ? yh : yl, 0.0);
  }
  res.x.back() = 0.0;
  res.s = s;
  res.u.resize(G - 1);
  res.t.assign(G, 0.0);
  for (std::size_t i = 0; i + 1 < G; ++i)
  {
    const double D = s[i + 1] - s[i];
    res.u[i] = (res.x[i + 1] - res.x[i]) / (2.0 * D);
    const double v = std::sqrt(std::max(res.x[i], 0.0)) + std::sqrt(std::max(res.x[i + 1], 0.0));
    res.t[i + 1] = res.t[i] + 2.0 * D / std::max(v, kTiny);
  }
  res.duration = res.t.back();
  res.ok = true;
  (void)path;
  (void)link;
  return res;
}

// Largest fraction of a limit used inside interval i (at 1/4, 1/2 and 3/4).
double excess(const CubicPath& path, const LinkDerivatives& link, const Limits& lim,
              const Result& r, std::size_t i)
{
  double worst = 0.0;
  for (double f : { 0.25, 0.5, 0.75 })
  {
    const double s = r.s[i] + f * (r.s[i + 1] - r.s[i]);
    const double x = std::max(r.x[i] + f * (r.x[i + 1] - r.x[i]), 0.0);
    const PointData d = evaluate(path, link, lim, s);
    for (int j = 0; j < d.q1.size(); ++j)
    {
      worst = std::max(worst, std::abs(d.q1[j] * r.u[i] + d.q2[j] * x) / lim.max_acceleration[j]);
      worst = std::max(worst, std::abs(d.q1[j]) * std::sqrt(x) / lim.max_velocity[j]);
    }
    if (lim.max_linear_acceleration > 0.0)
      worst = std::max(worst, (d.p1 * r.u[i] + d.p2 * x).norm() / lim.max_linear_acceleration);
    if (lim.max_linear_speed > 0.0)
      worst = std::max(worst, d.p1.norm() * std::sqrt(x) / lim.max_linear_speed);
  }
  return worst;
}
}  // namespace

Result parameterize(const CubicPath& path, const LinkDerivatives& link, const Limits& lim,
                    const Options& opt)
{
  Result res;
  if (!path.valid())
  {
    res.message = "camino sin longitud";
    return res;
  }
  const double L = path.length();
  std::vector<double> s;
  for (int k = 0; k <= opt.grid_points; ++k)
    s.push_back(L * k / opt.grid_points);
  const auto& w = path.knots();
  for (std::size_t k = 0; k + 1 < w.size(); ++k)
    for (int j = 0; j < opt.waypoint_subdivisions; ++j)
      s.push_back(w[k] + (w[k + 1] - w[k]) * j / opt.waypoint_subdivisions);
  std::sort(s.begin(), s.end());
  std::vector<double> grid;
  for (double v : s)
    if (grid.empty() || v - grid.back() > 1e-9 * std::max(L, 1.0))
      grid.push_back(v);
  grid.back() = L;

  std::vector<PointData> pts;
  pts.reserve(grid.size());
  for (double v : grid)
    pts.push_back(evaluate(path, link, lim, v));

  for (int round = 0;; ++round)
  {
    res = solve(path, link, lim, grid, pts, opt.polygon_sides);
    res.rounds = round;
    if (!res.ok || round >= opt.refine_rounds)
      return res;
    std::vector<std::size_t> bad;
    for (std::size_t i = 0; i + 1 < grid.size(); ++i)
      if (excess(path, link, lim, res, i) > 1.0 + opt.refine_tolerance)
        bad.push_back(i);
    if (bad.empty())
      return res;
    std::vector<double> new_grid;
    std::vector<PointData> new_pts;
    std::size_t b = 0;
    for (std::size_t i = 0; i < grid.size(); ++i)
    {
      new_grid.push_back(grid[i]);
      new_pts.push_back(pts[i]);
      if (b < bad.size() && bad[b] == i)
      {
        for (double f : { 0.25, 0.5, 0.75 })
        {
          const double v = grid[i] + f * (grid[i + 1] - grid[i]);
          new_grid.push_back(v);
          new_pts.push_back(evaluate(path, link, lim, v));
        }
        ++b;
      }
    }
    grid.swap(new_grid);
    pts.swap(new_pts);
  }
}

void Result::sample(const CubicPath& path, double dt, std::vector<double>& times,
                    std::vector<Eigen::VectorXd>& pos, std::vector<Eigen::VectorXd>& vel,
                    std::vector<Eigen::VectorXd>& acc) const
{
  times.clear();
  pos.clear();
  vel.clear();
  acc.clear();
  for (double tt = 0.0; tt < duration - 1e-9; tt += dt)
    times.push_back(tt);
  times.push_back(duration);
  std::size_t k = 0;
  for (double tt : times)
  {
    while (k + 1 < u.size() && t[k + 1] <= tt)
      ++k;
    const double tau = tt - t[k];
    const double sd0 = std::sqrt(std::max(x[k], 0.0));
    const double ss = std::clamp(s[k] + sd0 * tau + 0.5 * u[k] * tau * tau, 0.0, path.length());
    const double sd = std::max(sd0 + u[k] * tau, 0.0);
    const Eigen::VectorXd d1 = path.dq(ss);
    pos.push_back(path.q(ss));
    vel.push_back(d1 * sd);
    acc.push_back(path.ddq(ss) * sd * sd + d1 * u[k]);
  }
  pos.back() = path.q(path.length());
  vel.back().setZero();
}

}  // namespace cartesian_limit_adapter
