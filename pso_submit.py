# Script to submit a batch of particle 
# swarm optimisation jobs.
# Script will keep the number of jobs in the 
# queue around a target number. 
# Guy Billings, UCL, 2011.

import subprocess
import math
import sys
import os.path
import time
import user_variables as uv

def qsize(jobstr):
  # Returns number of queued jobs  
  os.system('qstat -r  | grep '+jobstr+' | wc -l > qj_'+jobstr)
  time.sleep(1)
  f=open('qj_'+jobstr,'r')
  numq=f.read()
  return int(numq)
  
def output(message):
  sys.stdout.write(message+'\n')
  sys.stdout.flush()

now = time.localtime(time.time())
tstr= time.strftime("%Y-%m-%d.%H-%M-%S", now)

output('Submitting cost model jobs at '+tstr)


# Construct index over all requested jobs
index=[]

deltadex=int((uv.deltamax-uv.deltaini)/uv.deltainc)+1
xdex=int((uv.xmax-uv.xini)/uv.xinc)+1
ydex=int((uv.ymax-uv.xini)/uv.yinc)+1

for d_i in range(deltadex):
  delta=uv.deltaini+d_i*uv.deltainc
  for x_i in range(xdex):
    x=uv.xini+x_i*uv.xinc
    for y_i in range(ydex):
      y=uv.yini+y_i*uv.yinc
      index.append([x,y,delta])
njobs=len(index)
output('There are '+str(njobs)+' jobs.')

# Submit jobs

if njobs > uv.nodes_to_use:
  to_submit=uv.nodes_to_use
else:
  to_submit=njobs
  
absub=0
all_jobs_not_submitted=True

while all_jobs_not_submitted:
  
  remaining=njobs-absub
  if to_submit>remaining:
    to_submit=remaining
    
  for sub in range(to_submit):
    #Build configuration file
    config_file=uv.output_dir+uv.script_base+tstr+'_'+str(index[absub][0])+'_'+str(index[absub][1])+\
                '_'+str(index[absub][2])+'.m'
    cf=open(config_file,mode='a')
    cf.write(str(uv.prec)+'\n')
    cf.write(str(uv.timestps)+'\n')
    cf.write(str(uv.particles)+'\n')
    cf.write(str(uv.muini)+'\n')
    cf.write(str(uv.mumax)+'\n')
    cf.write(str(uv.gammaini)+'\n')
    cf.write(str(uv.gammamax)+'\n')
    cf.write(str(uv.dini)+'\n')
    cf.write(str(uv.dmax)+'\n')
    cf.write(str(uv.pini)+'\n')
    cf.write(str(uv.pinc)+'\n')
    cf.write(str(uv.pmax)+'\n')
    cf.write(str(index[absub][0])+'\n')
    cf.write(str(index[absub][1])+'\n')
    cf.write(str(index[absub][2])+'\n')
    cf.write(str(uv.omega)+'\n')
    cf.write(str(uv.rp)+'\n')
    cf.write(str(uv.rg)+'\n')
    cf.write(str(uv.statp)+'\n')
    cf.write(str(uv.iterations)+'\n')
    cf.write(uv.output_dir+'\n')
    cf.write(tstr+'_'+str(index[absub][0])+'_'+str(index[absub][1])+'_'+str(index[absub][2])+'\n')
    cf.close()
    
    #Submit job
    comstr="--mca btl_tcp_if_exclude lo,eth1"
    scrfile=uv.working_dir+'pso_'+tstr+'_'+str(index[absub][0])+'_'+str(index[absub][1])+'_'+str(index[absub][2])+'.sh'
    sc=open(scrfile,mode='w')
    sc.write(uv.mpi_dir+'mpirun -np '+str(slots)+' '+comstr+' '+uv.working_dir+'bpnpso '+config_file+'\n')
    sc.close()
    substr='qsub -v LD_LIBRARY_PATH -e '+uv.qsub_dir+' -o '+uv.qsub_dir+' -l fc=yes -pe '+uv.pe+' '+str(slots)+' -l h_rt='+\
            str(uv.seconds)+' -R y'+' -l h_vmem='+str(uv.mem)+'G'+' '+scrfile
    output(substr)                        
    subprocess.call(substr,shell=True)          
    absub=absub+1
    
  output(str(absub)+' Jobs have been submitted.')  
  time.sleep(10)    
  #Status check
  qs=qsize('pso_'+tstr)
  while qs>=uv.nodes_to_use:
    qs=qsize('pso_'+tstr)  
    time.sleep(2)
  to_submit=uv.nodes_to_use-qs
  
  if absub>=njobs:
    all_jobs_not_submitted=False
  

  


  
  
