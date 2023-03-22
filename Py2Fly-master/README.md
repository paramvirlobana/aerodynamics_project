# Welcome to Py2Fly!
This is a library of some useful routines for predicting airfoil and wing performance, generating airfoil shapes, and predicting aircraft performance.

You can run this interactive in the cloud without having to install a thing! Click below to run it in Binder:
https://mybinder.org/v2/gl/bvermeir%2FPy2Fly/HEAD

### Py2Fly Performance: A library for predicting general aircraft performance characteristics including:
 - Lift to drag ratio vs. airspeed
 - Thrust required vs. airspeed
 - Power required vs. airspeed
 - Rate of climb vs. airspeed for multiple throttle settings
 - Specific fuel consumption vs. airspeed
 - Specific range vs. airspeed

### Py2Fly Airfoils:
 - Airfoil coordinate and camberline generators for standard airfoil types

### Py2Fly Aerodynamics: A library for aerodynamic analysis using thin airfoil and vortex panel methods
 - Two-dimensional vortex panel method for predicting airfoil section performance including:
   - Lift coefficient vs. angle of attack
   - Pressure coefficient distribution vs. chord-wise position
 - Thin airfoil solver for predicting airfoil section performance including:
   - Lift coefficient vs. angle of attack
   - Moment coefficient vs. angle of attack
   - Center of pressure vs. angle of attack
   - Vortex line strength distribution along the camberline
 - Finite span wing solver for predicting lift distribution and induced drag
   - Lift coefficient, induced drag coefficient, and Oswald efficiency factor
   - Lift coefficient vs. angle of attack
   - Induced drag coefficient vs. angle of attack
   - Drag polar plot
   - Circulation vs. span