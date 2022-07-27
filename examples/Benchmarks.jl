using Revise
using Geom4hep

const T = Float64

box = Box{T}(1,2,3)
benchmarkShape(box, Npoints=10240, Repetitions=1000)

trd = Trd{T}(3, 4, 5, 6, 7)
benchmarkShape(trd, Npoints=10240, Repetitions=1000)

trap = Trap{T}(3, 4, 5, 6, 7)
benchmarkShape(trap, Npoints=10240, Repetitions=1000)

cone1 = Cone{T}(5, 10, 7, 15, 10, 0, 2π)
benchmarkShape(cone1, Npoints=10240, Repetitions=1000)

cone2 = Cone{T}(5, 10, 7, 15, 10, π, π/2)
benchmarkShape(cone2, Npoints=10240, Repetitions=1000)

tube1 = Tube{T}(0, 5, 10, 0,2π)
benchmarkShape(tube1, Npoints=10240, Repetitions=1000)

tube2 = Tube{T}(2, 5, 10, 0,2π)
benchmarkShape(tube2, Npoints=10240, Repetitions=1000)

tube3 = Tube{T}(2, 5, 10, 0, 2π/3)
benchmarkShape(tube3, Npoints=10240, Repetitions=1000)

pcone = Polycone{T}([0.1, 0., 0., 0.2], [1., 2., 2., 1.5], [-1, -0.5, 0.5, 10], 0., 2π)
benchmarkShape(pcone, Npoints=10240, Repetitions=1000)
