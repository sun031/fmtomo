#!/usr/bin/env python


import commands
import os, sys
import time
import numpy as np

def srclineno(isrc):
    # return int(9*isrc-8), int(9*isrc+1)
   return int(6*isrc-5), int(6*isrc+1)

def reclineno(irec):
    return int(4*irec-3), int(4*irec+1)

def pfm3d(ncpu=12):

    print "number of CPUs or Cores:\t",ncpu

    with open("sources.in", "r") as fpsrc:
        srclst = fpsrc.readlines()

    with open("receivers.in", "r") as fprec:
        reclst = fprec.readlines()

    nsrc = int(srclst[0])
    print "number of sources:\t", nsrc

    pnsrc = nsrc/ncpu+1

    if nsrc-pnsrc*ncpu > ncpu/2:
        pnsrc = pnsrc+1

    print "number of sources for parallel:\t", pnsrc

    # start and end shot number for each CPU core
    psrc1 = np.zeros(ncpu)
    psrc2 = np.zeros(ncpu)
    for i in range(ncpu):
        psrc1[i] = i*pnsrc+1
        psrc2[i] = (i+1)*pnsrc

    psrc2[-1]=nsrc

    print "start and end shot number for each CPU core"
    print psrc1
    print psrc2

    prec1 = np.zeros(ncpu)
    prec2 = np.zeros(ncpu)
    for i in range(3,len(reclst),4):
#         print i,int(reclst[i]), i-2, i+2

        srcno = int(reclst[i])

        flag1 = True
        for icpu in range(ncpu):
            if srcno==psrc2[icpu]:
                prec2[icpu] = int(i+2)
                break



    prec1[0] = int(1)
    for icpu in range(1,ncpu):
        prec1[icpu] = int(prec2[icpu-1])

#     print prec1
#     print prec2


#         if i>50:
#             os._exit(0)

    fpsave = open("cpu", "w")
    for i in range(ncpu):
        pnsrc = int(psrc2[i]-psrc1[i]+1)

        i2s = str(i).zfill(2)

        if not os.path.exists(i2s):
            os.system("mkdir %s" % i2s)

        index1, index2 = srclineno(psrc1[i])
        index3, index4 = srclineno(psrc2[i])

        srcstr = ""
        fpsrc = open(i2s+"/sources.in", "w")
        srcstr += str(pnsrc)+"\n"

        fpsrc.write("%s" % srcstr)
        fpsrc.writelines(srclst[index1:index4])

        fpsrc.close()

        fprec = open(i2s+"/receivers.in", "w")
        nrec = int((prec2[i]-prec1[i])/4)
        fprec.write(str(nrec)+"\n")

        index1 = int(prec1[i])
        index2 = int(prec2[i])
        #print index1, index2

        lst = reclst[index1:index2]
        for irec in range(2,len(lst),4):

            isrc = int(lst[irec])-psrc1[i]+1

            lst[irec] = "       "+str(int(isrc))+"\n"

        fprec.writelines(lst)

        fprec.close()

        fpsave.write("%d\t%d\t%d\t%d\n" % (psrc1[i],psrc2[i],prec1[i],prec2[i]))

        cpfiles = "propgrid.in interfaces.in vgrids.in mode_set.in frechet.in"

        os.system("cp %s %s" % (cpfiles, i2s))


    fpsave.close()
    pass

def pfm3d_run(ncpu=12, name="fm3d"):
    for i in range(ncpu):
        i2s = str(i).zfill(2)


        os.chdir("./"+i2s)
        print os.getcwd()
        os.system("nohup %s &" % name)
        os.chdir("../")
        print os.getcwd()

    pass

def pfm3d_assem(ncpu=12):

    with open("cpu", "r") as fp:
        lst = fp.readlines()

    psrc1 = np.zeros(len(lst), dtype=np.int)
    for i in range(len(lst)):
        row = lst[i].split()
        psrc1[i] = int(row[0])

    print psrc1

    fparr = open("arrivals.dat", "w")

    counter = 0
    for i in range(ncpu):

        i2s = str(i).zfill(2)
        with open(i2s+"/arrivals.dat") as fp:
            lst = fp.readlines()

        for line in lst:
            counter += 1
            row = line.split()

            isrc = int(row[1])+psrc1[i]-1

            fparr.write("%6d%6d%6s%6s%15s%5s%5s\n" %(counter, isrc, row[2], row[3], row[4], row[5], row[6]))

    fparr.close()


def pfm3d_assem_frechet(ncpu=12):
    """
    assemble frechet.dat
    :param ncpu:
    :return:
    """

    with open("cpu", "r") as fp:
        lst = fp.readlines()

    psrc1 = np.zeros(len(lst), dtype=np.int)
    for i in range(len(lst)):
        row = lst[i].split()
        psrc1[i] = int(row[0])

    print psrc1

    fpfrechet = open("frechet.dat", "w")

    counter = 0
    nrec = 0
    for i in range(ncpu):

        i2s = str(i).zfill(2)
        with open(i2s+"/frechet.dat") as fp:
            lst = fp.readlines()


        # iline = 0
        for line in lst:
            row = line.split()
            if len(row) == 5:
                idrec = int(row[0])
                idsrc = int(row[1])
                idray = int(row[2])
                npath = int(row[3])
                npdev = int(row[4])

                idsrc += psrc1[i] - 1
                idrec += nrec


                fpfrechet.write("%d\t%d\t%d\t%d\t%d\n" % (idrec, idsrc, idray, npath, npdev))
            else:
                fpfrechet.write("%s" % line)

        nrec = idrec


    fpfrechet.close()









def pfm3d_fast(ncpu=12, cmd="fm3d"):

    output = commands.getoutput("which "+cmd)
    output = output.strip()

    tt = time.strftime("%Y%m%d%H%M%S",time.localtime(time.time()))

    nickname = cmd+"_"+tt

    print "Alias name:", nickname

    # copy fm3d into current directory
    os.system("cp %s %s" % (output, nickname))
    os.system("chmod u+x %s" % nickname)

    pfm3d(ncpu=ncpu)
    pfm3d_run(ncpu=ncpu, name="../"+nickname)

    # monitor the process, after it finished then run assem
    while True:
        ps_string = os.popen('pidof %s' % nickname,'r').read()

        if len(ps_string)==0:
            print "No related process. EXIT!"
            break

        ps_strings = ps_string.strip().split('\n')

        print ps_strings

        if len(ps_strings)>0:
            # os.system('sudo halt')
            # print ps_strings
            # for pid in ps_strings:
                # os.system("kill -9 %s" % pid)
                # print "PID", pid, "killed"
            print "Process " + nickname + " exist..."
            print "Sleeping 10 sec and then check..."
            time.sleep(10)
            continue

        else:
            print "No related process. EXIT!"
            break

    # sleep 30 seconds in order to make sure the data writing done
    time.sleep(30)
    pfm3d_assem(ncpu=ncpu)
    pfm3d_assem_frechet(ncpu=ncpu)

    # os.system("rm %s" % nickname)




    pass


if __name__=="__main__":

    pfm3d_fast(ncpu=24, cmd="fm3d")

    pass
