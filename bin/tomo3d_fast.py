#!/usr/bin/env python
#
# NOTE: modifying to the python script.
#
############################################
# Script for running multi-parameter
# tomography program
############################################
#
import os,sys
import commands
import pfm3d_fast

# number of CPU used to run the tomography
ncpu = 32



########################################################
# Program and files for solving the forward problem with
# the Fast Marching Method
########################################################
#
# Name of program for calculating FMM traveltimes
fmm='fm3d'
# Name of program for computing frechet.in file
frech='frechgen'
# Name of file containing FMM traveltimes
ttim='arrivals.dat'
# Name of file containing current velocity grid
cvg='vgrids.in'
# Name of file containing current interface grid
cig='interfaces.in'
# Name of file containing current source locations
csl='sources.in'
# Output diagnostics file
diagnos='fm3dlog.out'
#
#########################################################
# Program and files for solving the inverse problem using
# subspace inversion
#########################################################
#
# Name of program for performing inversion
inv='invert3d'
# Name of file containing iteration number
itn='inviter.in'
# Name of file containing current model traveltimes
mtrav='mtimes.dat'
# Name of file containing reference model traveltimes
rtrav='rtimes.dat'
# Name of file containing initial velocity grid
ivg='vgridsref.in'
# Name of file containing initial interface grid
iig='interfacesref.in'
# Name of file containing initial source locations
isl='sourcesref.in'
#
#########################################################
# Program and files for calculating traveltime
# residuals
#########################################################
#
# Name of program for calculating traveltime residuals
resid='residuals'
# Name of output file for calculating traveltime residuals
resout='residuals.dat'
#
###############################
###############################
# cd=chdir
ifile='tomo3d.in'
#
###############################
###############################
#
fp = open(ifile, 'r')
lines = fp.readlines()
NI = int(lines[0].strip())
BGITER = int(lines[1].strip())
BINV = int(lines[2].strip())
fp.close()
#
# If necessary, copy the initial velocity,
# interface and source files to the current files
# read by fm3d. Generate the frechet.in file and 
# set the iteration number to 1.
#
if BGITER == 0:
    os.system('cp %s %s' % (ivg, cvg))
    os.system('cp %s %s' % (iig, cig))
    os.system('cp %s %s' % (isl, csl))
    os.system(frech)
    ITER=1
    fp = open(itn, 'w')
    fp.write(str(ITER))
    fp.close()
    
    
        
#
#
# Run FMM once to generate traveltimes for current model if
# necessary
#
if BGITER==0 or (BGITER==1 and BINV==1 ):
    # os.system(fmm+' > '+diagnos)
    pfm3d_fast.pfm3d_fast(ncpu=ncpu, cmd=fmm)
    os.system('cp %s %s' % (ttim, mtrav))
   #
   #  Repeat this copy for reference traveltimes
   #
    os.system('cp %s %s' % (ttim, rtrav))
   #
   # Calculate initial residual
   #
    fp = open(resout, 'a')
    output = commands.getoutput(resid)
    fp.write(output+'\n')
    fp.close()

#
# Now begin a loop to iteratively apply subspace inversion
# and FMM
#
try:
    fp = open(itn, 'r')
    line = fp.readline()
    fp.close()
    ITER=int(line.strip())
except:
    ITER=1

while ITER<=NI:
    os.system(inv)
    
    for line in ['interfaces.in', 'vgrids.in', 'sources.in']:
        os.system('cp %s %s' % (line, line+'.'+str(ITER).zfill(3)))
    
    # os.system(fmm+' > '+diagnos)
    pfm3d_fast.pfm3d_fast(ncpu=ncpu, cmd=fmm)
    os.system('cp %s %s' % (ttim, mtrav))
#
#  Calculate traveltime residual
#
    fp = open(resout, 'a')
    output = commands.getoutput(resid)
    fp.write(output+'\n')
    fp.close()
    ITER=ITER+1
    COUNT=ITER
    
    fpitn = open(itn, 'w')
    fpitn.write(str(COUNT))
    fpitn.close()
    
for line in ['interfaces.in', 'vgrids.in', 'sources.in']:
    os.system('cp %s %s' % (line, line+'.'+str(ITER).zfill(3)))
    
    
    
    
    
    
    
    
    
    
    
    
    
