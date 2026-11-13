^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Changelog for package cartesian_limit_adapter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

0.1.0 (2026-11-13)
------------------
* First release.
* ``cartesian_limit_adapter/LimitCartesianSpeed`` planning response adapter: TOPP-RA
  retiming with joint velocity/acceleration limits plus a bound on the linear speed and
  linear acceleration of one link.
* Honors ``MotionPlanRequest.max_cartesian_speed`` and ``cartesian_speed_limited_link``.
