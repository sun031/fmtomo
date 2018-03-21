#!/usr/bin/env python

ttr1 = "ts1101559.ttr"
ttr2 = "picks/ts1101559.ttr"

with open(ttr1, "r") as fp:
    lines = fp.readlines()


fp = open(ttr2, "w")

nsta = float(lines[0].strip())
row = lines[1].split()
evla = float(row[0])
evlo = float(row[1])
evdp = float(row[2])
ph = lines[6].strip()

nsta = 0
lst = []
for line in lines[8:]:
    row = line.split()

    if int(row[5])==0:
        continue

    nsta += 1
    stla = nsta
    stlo = nsta
    stel = nsta
    uncertainty = float(row[3])

    shift1 = float(row[4])
    shift2 = float(row[2])

    lst.append([stla, stlo, stel, shift1, uncertainty])

print nsta
print len(lst)

fp.write("%d\n" % nsta)
fp.write("%.4f\t%.4f\t%.4f\n" % (evla, evlo, evdp))
fp.write("%s\n" % ph)
for a in lst:
    fp.write("%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n" % (a[0], a[1], a[2], a[3], a[4]))
fp.close()