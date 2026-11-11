using Geom4hep

const T = Float64

box = Box{T}(1,2,3)
benchmark_shape(box, npoints=10240, repetitions=1000)

trd = Trd{T}(3, 4, 5, 6, 7)
benchmark_shape(trd, npoints=10240, repetitions=1000)

trap = Trap{T}(3, 4, 5, 6, 7)
benchmark_shape(trap, npoints=10240, repetitions=1000)

cone1 = Geom4hep.Cone{T}(5, 10, 7, 15, 10, 0, 2π)
benchmark_shape(cone1, npoints=10240, repetitions=1000)

cone2 = Geom4hep.Cone{T}(5, 10, 7, 15, 10, π, π/2)
benchmark_shape(cone2, npoints=10240, repetitions=1000)

tube1 = Tube{T}(0, 5, 10, 0,2π)
benchmark_shape(tube1, npoints=10240, repetitions=1000)

tube2 = Tube{T}(2, 5, 10, 0,2π)
benchmark_shape(tube2, npoints=10240, repetitions=1000)

tube3 = Tube{T}(2, 5, 10, 0, 2π/3)
benchmark_shape(tube3, npoints=10240, repetitions=1000)

pcone = Polycone{T}([0.1, 0., 0., 0.2], [1., 2., 2., 1.5], [-1, -0.5, 0.5, 10], 0., 2π)
benchmark_shape(pcone, npoints=10240, repetitions=1000)
