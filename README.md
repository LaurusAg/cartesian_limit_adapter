# cartesian_limit_adapter

A MoveIt 2 planning response adapter that limits the **linear speed and linear acceleration of one link** (usually the tool) on every planned trajectory. It keeps the joint velocity and acceleration limits too.

MoveIt's `MotionPlanRequest` has `max_cartesian_speed` and `cartesian_speed_limited_link` fields, but no MoveIt 2 component applies them (checked on Jazzy, 2.12). This adapter applies them, and it can also take its limits from parameters.

## How it works

The adapter takes the path of the planned trajectory, the waypoints left by the planner and `AddTimeOptimalParameterization`, and retimes it with TOPP-RA (Pham & Pham, 2018). The geometric path does not change; only the timing does.

- **Path:** a natural cubic spline through the waypoints, parameterized by joint-space chord length.
- **Joint limits:** velocity and acceleration of every joint, read from the robot model (URDF + `joint_limits.yaml`) and scaled by `max_velocity_scaling_factor` / `max_acceleration_scaling_factor`, the same way TOTG does.
- **Link speed:** `‖J_v q'‖ ṡ ≤ v_max`.
- **Link acceleration:** `‖p' s̈ + p'' ṡ²‖ ≤ a_max`, with `p' = J_v q'` and `p'' = J_v q'' + J̇_v q'`.
  - That vector lies in the plane of `{p', p''}`, so its norm is bounded exactly by a polygon inscribed in that plane.
  - The bound is conservative: it never exceeds the limit and gives up at most `1 − cos(π/K)` of it (3.4% with the default K = 12).
- **Discretization:** constraints are enforced at both ends of every grid interval. The grid is uniform plus the waypoints. After each solve the adapter checks every interval at ¼, ½ and ¾ and splits the ones that exceed a limit by more than `refine_tolerance`, then solves again.
- **Output:** the trajectory is resampled every `resample_dt`, with positions, velocities and accelerations.

Computation takes a few milliseconds per trajectory (5–20 ms on a Jetson Orin Nano for a 7-joint arm).

## Usage

Add the adapter to the response adapters of your planning pipeline, after `AddTimeOptimalParameterization` and before `ValidateSolution`, and set its parameters under `cartesian_limit`. For example, in `ompl_planning.yaml`:

```yaml
response_adapters:
  - default_planning_response_adapters/AddTimeOptimalParameterization
  - cartesian_limit_adapter/LimitCartesianSpeed
  - default_planning_response_adapters/ValidateSolution
  - default_planning_response_adapters/DisplayMotionPath
cartesian_limit:
  link: tool0
  max_linear_speed: 2.0         # m/s
  max_linear_acceleration: 4.0  # m/s^2
```

The parameters live in the pipeline's namespace, for example `ompl.cartesian_limit.max_linear_speed`. They are read on every request, so they can be changed at runtime:

```bash
ros2 param set /move_group ompl.cartesian_limit.max_linear_speed 0.5
ros2 param set /move_group ompl.cartesian_limit.enabled false
```

A request with `max_cartesian_speed > 0` overrides `max_linear_speed`. A non-empty `cartesian_speed_limited_link` overrides `link`.

### Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `enabled` | `true` | Turn the adapter off without removing it from the pipeline. |
| `link` | `""` | Link to limit. Empty: the last link of the planning group. |
| `max_linear_speed` | `0.0` | m/s. `<= 0` means no bound. |
| `max_linear_acceleration` | `0.0` | m/s². `<= 0` means no bound. |
| `resample_dt` | `0.02` | Seconds between output waypoints. |
| `grid_points` | `100` | Uniform part of the grid. |
| `waypoint_subdivisions` | `2` | Grid points per waypoint interval. |
| `polygon_sides` | `12` | Sides of the polygon that bounds the acceleration norm. |
| `refine_rounds` | `3` | Maximum grid refinements. |
| `refine_tolerance` | `0.005` | Allowed excess (fraction of a limit) before an interval is split. |

With both limits at `0` (the default) and no `max_cartesian_speed` in the request, the adapter leaves the trajectory untouched.

## Behavior

- If the trajectory does not move (already at the goal), it is passed through.
- If TOPP-RA finds no feasible timing, the plan fails with `FAILURE`. The faster trajectory is never passed on in its place.
- Joints without a velocity bound make the plan fail. Joints without an acceleration bound use 1.0, like TOTG.

## Limitations

- Only the linear part is limited. Angular speed and acceleration are not implemented yet.
- Jerk is not bounded: acceleration is piecewise constant, as with TOTG.
- The retimed path follows the TOTG waypoints through a spline. Keep `ValidateSolution` after this adapter so the result is collision-checked.

## Install

```bash
sudo apt install ros-jazzy-cartesian-limit-adapter
```

(Once the package is released to the ROS buildfarm. Until then, build it from source.)

## Build and test

```bash
colcon build --packages-select cartesian_limit_adapter
colcon test --packages-select cartesian_limit_adapter && colcon test-result --verbose
```

Tested on ROS 2 Jazzy with MoveIt 2.12.

## License

BSD-3-Clause. Copyright (c) 2026, Ignacio Gavier and LaurusAg.
