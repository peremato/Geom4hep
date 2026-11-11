"""
    Geom4hep.Units

System of units (Geant4 conventions): lengths in mm, angles in radians.
Use as `using Geom4hep.Units; Box{Float64}(1cm, 2cm, 3cm)`.
"""
module Units

export mm, cm, m, km, um, nm, radian, rad, mrad, degree, deg

#--- Length (internal unit: mm)
const mm = 1.0
const cm = 10mm
const m  = 1000mm
const km = 1000m
const um = 1e-3mm
const nm = 1e-6mm

#--- Angle (internal unit: radian)
const radian = 1.0
const rad = radian
const mrad = 1e-3rad
const degree = π / 180.0 * radian
const deg = degree

end # module Units
