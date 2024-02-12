import FLOWPanel as pnl

run_name = "rectangle"
save_path = run_name

paraview=true

AoA=5.0
magVinf = 10.0
Vinf = magVinf*[cos(AoA*pi/180), 0, sin(AoA*pi/180)]

rho = 1.225

b = 10.0
c = 4.0


bodytype = pnl.RigidWakeBody{pnl.VortexRing}
