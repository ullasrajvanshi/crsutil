# crsutil
Coordinate Reference System UTILities (crsutil.py)

# ----------------------------------------------------------------------------
# crsutil.py    CRSUTIL Coordinate and Time Reference System Toolbox.
# Version 1.0 (1 November 2019).
# Created by: Hans van der Maarel and Ullas Rajvanshi
# Date:       13 Nov 2019
# Modified:   -
#
#   Copyright: Hans van der Marel, Ullas Rajvanshi, Delft University of Technology
#   Email:     H.vanderMarel@tudelft.nl ; U.Rajvanshi@student.tudelft.nl
#   Github:    www.github.com/ullasrajvanshi
# ----------------------------------------------------------------------------
# Functions created:
# Coordinate transformations (ECEF reference frame)
#   inqell      - Semi-major axis, flattening and GM for various ellipsoids
#   xyz2plh     - Cartesian Coordinates to Ellipsoidal coordinates
#   plh2xyz     - Ellipsoidal coordinates to Cartesian Coordinates
#   xyz2neu     - North, East, Up (dN, dE, dU) to Cartesian delta's (dX, dY, dZ)
#   neu2xyz     - Cartesian delta's (dX, dY, dZ) to North, East, Up (dN, dE, dU)
#   plh2neu     - Ellipsoidal (Lat,Lon,Hgt) to North,East,Up (dN, dE, dU)
#   xyz2zas     - Cartesian coordinates to Zenith angle, azimuth and distance
#   zas2xyz     - Zenith angle, azimuth and distance to cartesian coordinates
#
# UT1 to GMST, and ECI/ECEF, conversions
#
#   ut2gmst    - Compute Greenwich Mean Siderial Time from UT1
#   ecef2eci   - Convert position and velocity from ECEF to ECI reference frame
#   eci2ecef   - Convert position and velocity from ECI to ECEF reference frame
#
# Keplerian elements
#
#   vec2orb     - Convert inertial state vector into Keplerian elements
#   orb2vec     - Convert Keplerian elements into inertial state vector
#   kepler      - Compute mean anomaly from eccentric anomaly (Kepler's equation)
#   keplernu    - Compute mean anomaly from true anomaly (Kepler's equation)
#   keplerm     - Compute eccentric/true from mean anomaly solving Kepler's eqn
#
# ----------------------------------------------------------------------------
