- Before changes
  ```
  @time generateXRay(nav, world, 1e6, 1)
  272.065976 seconds (2.91 G allocations: 72.223 GiB, 2.17% gc time)
  ```
- With changes from Philippe G.
  ```
  @time generateXRay(nav, world, 1e6, 1)
  223.448555 seconds (10 allocations: 7.632 MiB)
  ```

./test/XRayBenchmarkFromROOTFile ../Geom4hep/examples/cms2018.gdml MB x 100 --novoxel
generateXRay(nav, world, 100, 1)


| program |  options            |  time | allocations | # steps | # distanceToIn | Comments |
|---------|---------------------|-------|-------------|---------|----------------|----------|
| VecGeom |  NewSimpleNavigator | 8.1 s |             |  1.8 M  | 88 M           | |
| VecGeom |  default            | 1.27 s|             |  1.8 M  | 334 k ???      | |
| VecGeom |  --use-bvh-navigator | 1.04 s|            |  1.8 M  | 2.3 M          | |
| Geom4Hep | TrivialNavigator   | 4.6 s | 10          | 1.8 M |  | |   
| Geom4Hep | BVHNavigator       | 1.6 s |  10         |  1.8 M  | 3.48 M         | |




