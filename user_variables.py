# This version for large networks. Fixed network size and loop over 
# values of d and pin.
# User variables:

#Simulation parameters
nodes_to_use=30 #Determines maximum number jobs to place in the queue
qsub_dir='~/psojobs' #Directory into which qsub places its output files (must exist)
working_dir='/home/ucgbgbi/data/BPN_v0.5/Cnet/BPNPSO/' #Directory in which these files are located
output_dir='/home/ucgbgbi/data/BPN_v0.5/Cnet/BPNPSO/output/' #Directory to save files to
mpi_dir='/opt/SUNWhpc/HPC8.2.1/gnu/bin/' #Location of MPI
script_base='pso_config_' #Base for mathematica file names
seconds=144000# seconds
mem=5 # Memory requested to run job (Gb)
numprocs=2 #Number of processors per Swarm
pe='openmpi' #-pe setting for MPI

#Parameters to be directly placed in PSO configuration file
prec=200
timestps=1000000000
particles=10
muini=20
mumax=100
gammaini=20
gammamax=100
dini=1
dmax=20
pini=0.005
pmax=0.995
pinc=0.005
omega=1
rp=0.1
rg=0.1
statp=0.7
iterations=5

#Scan over cost parameters
deltaini=0.1
deltamax=0.2
deltainc=0.1
xini=1
xmax=1
xinc=1
yini=1
ymax=1
yinc=1
