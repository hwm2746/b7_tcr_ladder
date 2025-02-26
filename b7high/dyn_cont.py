from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from datetime import datetime

iprefx = '1bd2i_incr8_1' # output prefix
oprefx = '1bd2i_incr8_2' # output prefix

simtype='nvt' # choose between cpt & nvt
# nose is preferred. Avoid andersen or langevin
thermo0='nose' # andersen, nose, or langevin

istep=25000000
nsavc=10000

# psf/crd should precede params declaration
psf = CharmmPsfFile('1bd2i.psf') 

# set chmdir for your system
chmdir='/scratch/user/acchangg12/charmm/c47a11/toppar/'
params = CharmmParameterSet(chmdir+'par_all36_prot.prm',chmdir+'top_all36_prot.rtf',chmdir+'toppar_water_ions.str')

# setBox: in nm unit. should precede createSystem. 
# For ortho box, use Lx, Ly, Lz separately. From *d0*.rst 
# from shift.inp 
Lx = 0.1 * 227
Ly = 0.1 * 88
Lz = 0.1 * 81.6426699990669

psf.setBox(Lx, Ly, Lz)

# 1.2nm: Cutoff used in charmm
system = psf.createSystem(params, nonbondedMethod=PME,
                          nonbondedCutoff=1.2*nanometer, switchDistance=1.0*nanometer,
                          constraints=HBonds, ewaldErrorTolerance=1.e-4)

#######################################################
# CONS HARM . apply constraints to end residues about a position
# hf: 1.0 kcal/[mol*A^2] ; 1 kcal/[mol*A^2] = 418.40 kJ/[mol*nm^2]

force = CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"); # harmonic

force.addPerParticleParameter("x0")
force.addPerParticleParameter("y0")
force.addPerParticleParameter("z0")

hf0 = 1.0 # spring constant in charmm
hf = hf0*KJPerKcal*(1/(NmPerAngstrom*NmPerAngstrom))
print(hf)  # in job output file 
force.addGlobalParameter("k", hf); # spring const

# index for CA atoms (counts from 0)
# the following from dcd/out_dyn0.dat 
mhc_calpha=4463 # pull_atom, hold_atom CA 
tcr_calpha=[9613, 13563] # hold_atom [0]==alpha, [1]==beta, hmcm 

d0 = [113.500000 ,  44.000000 , 40.821335] # "translation vector" from out_shift.dat
d0[0]  = NmPerAngstrom * d0[0]
d0[1]  = NmPerAngstrom * d0[1]
d0[2]  = NmPerAngstrom * d0[2]
mhc_x0 = NmPerAngstrom * -95.021356 # start position from out_get_pos.dat
tcr_x0 = NmPerAngstrom * 95.2036753 # centroid pos tcr hold_atom from out_get_pos.dat

mhc_r0 = [ d0[0]+mhc_x0, d0[1], d0[2] ] 
tcr_r0 = [ d0[0]+tcr_x0, d0[1], d0[2] ] 
print(mhc_r0)
print(tcr_r0) 

system.addForce(force) # comment to remove cons harm 
force.addParticle(mhc_calpha,[])
force.setParticleParameters(0,mhc_calpha,mhc_r0)

##########
## add hmcm to tcr end
hmcm0 = 'k*((x1-x0)^2+(y1-y0)^2+(z1-z0)^2);'
hmcm0 += 'k  = %f;' % hf #418.4;' # spring const
# hmcm0 += 'k = 418.4;' # spring const
hmcm0 += 'x0 = %f;' % tcr_r0[0]
hmcm0 += 'y0 = %f;' % tcr_r0[1]
hmcm0 += 'z0 = %f;' % tcr_r0[2]

hmcm = CustomCentroidBondForce(1, hmcm0)
hmcm.addGroup(tcr_calpha) # particle mass used as weights 
hmcm.addBond([0])

system.addForce(hmcm) 

#######################################################
# RESD
# distance = 10A between tcr_calpha residues 
# add to the tcr_calpha residues

ca0=[tcr_calpha[0]]; ca1=[tcr_calpha[1]];

d1 = NmPerAngstrom * 10  # resd distance

# max: apply force only when distance is > r0
resd = CustomCentroidBondForce(2,"k*max(0,distance(g1,g2)-r0)^2")
resd.addPerBondParameter("k");
resd.addPerBondParameter("r0");
resd.addGroup(ca0)
resd.addGroup(ca1)
bondGroups= [0,1]
bondParameters= [ hf, d1]  # [k, r0]
resd.addBond(bondGroups, bondParameters);

system.addForce(resd) 

####################################
platform = Platform.getPlatformByName('CUDA') # for CUDA
properties = {'CudaPrecision': 'mixed','DeterministicForces':'true'} # for double precision

####################################
# Choose between CPT and NVT
# http://docs.openmm.org/latest/userguide/application.html#integrators
#################################

if simtype=='cpt':
    system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))

if thermo0=='andersen':
#    system.addForce(AndersenThermostat(300*kelvin, 25/picosecond)) # for testing
    system.addForce(AndersenThermostat(300*kelvin, 1/picosecond))
    integrator = VerletIntegrator(0.002*picoseconds)
elif thermo0 == 'nose':
    integrator = NoseHooverIntegrator(300*kelvin, 1/picosecond,
                                0.002*picoseconds);    
else: # langevin
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)

odcd = 'dcd/'+oprefx+'.dcd'
ofn =  'out/'+oprefx+'.dat'

chkpt0 =  'chkpt/'+iprefx+'.dat' # in openmm, rst is in xml format
chkpt =  'chkpt/'+oprefx+'.dat' # in openmm, rst is in xml format

####################################
simulation = Simulation(psf.topology, system, integrator,platform,properties)
simulation.loadCheckpoint(chkpt0)

simulation.reporters.append(DCDReporter(odcd, nsavc))
simulation.reporters.append(StateDataReporter(ofn, nsavc, step=True,
    time=True,potentialEnergy=True, temperature=True, totalEnergy=True,\
    volume=True,totalSteps=istep,kineticEnergy=True,separator=' '))

t0=datetime.now() # measure elapsed time
simulation.step(istep)
t1=datetime.now()

simulation.saveCheckpoint(chkpt) # saving state: Sec 3.6.15 of user guide
ffo=open(ofn,'a')
ffo.write('# Elapsed time  (hh:mm:ss.ms) {}'.format(t1-t0))
ffo.close()
