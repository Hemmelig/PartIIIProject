#!/usr/bin/env python3

import os.path, sys

event = sys.argv[1]

latmin = 48
latmax = 70
lonmin = -167
lonmax = -130

outfile = 'receivers_forcsem_' + event + '.dat'
f = open(outfile, 'w')

f.write('Number of stations is:\n')
f.write('%d \n' % ((latmax - latmin + 1) * (lonmax - lonmin + 1)))
f.write('     nw stn lat lon:\n')
for i in range((lonmax - lonmin + 1)):
    slon = lonmin + (i * 1)
    if i < 10:
        ii = str(0) + str(i)
    else:
        ii = str(i)
    for j in range((latmax - latmin + 1)):
        if j < 10:
            jj = str(0) + str(j)
        else:
            jj = str(j)
        sta  = ii + jj
        nw   = 'AA'
        slat = latmin + (j * 1)

        f.write('%s %s  %2.3f  %3.3f\n'% (nw, sta, slat, slon))

f.close()
