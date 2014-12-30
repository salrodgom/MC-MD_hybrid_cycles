cutp 12.0 1.0
integrator leapfrog verlet
ensemble npt 0.05 0.05
temperature  TEMPERATURE
mdmaxtemp    1000.0
tscale       PRODUEQUI 0.0002 ps
equil        EQUILIBRATION  ps
produ        PRODUCTION ps
timestep     TIMESTEP ps
sample       0.001 ps
write        0.001 ps
output cif        md.cif
dump   every   1  md.res
output trajectory md.trg
output movie arc  md
