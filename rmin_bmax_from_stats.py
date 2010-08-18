#!/usr/bin/env python
#
# Steven Bickerton
# Dept. of Physics/Astronomy, McMaster University
# bick@physics.mcmaster.ca
# Made with makeScript, Tue Mar 18, 2008  16:05:11 DST
# Host: bender.astro.princeton.edu
# Working Directory: /Users/bick/working/elongation_test


import sys, re
import numpy as num
from scipy.optimize import leastsq

def Usage():
    print """Usage: %s statfile """ % sys.argv[0]
    sys.exit(1)
    
def invPrit(a, x0, prob, sr):
    return (1.0/a) * num.tan(num.pi*(prob-0.5)) + x0 + sr

def pritchet(p, x):
    wid, x0 = p
    return  (1.0/num.pi) * num.arctan( wid*(x - x0) ) + 0.5

def pritchet2d(p, r, b, rbreak, rwid, bwid):
    # introduced a linear term 's' to allow ip to increas with radius
    #  ... beginning at r=rbreak (pritchet acts as heaviside)
    A, r0, b0 = p
    maxr = num.max(r)
    maxb = num.max(b)

    Amin, Amax = 0.0, 1.0
    #rwidmin, rwidmax = 5.0, 200.0
    #bwidmin, bwidmax = 20.0, 500.0
    r0min, r0max     = 10.0, maxr
    b0min, b0max     = 0.0, maxb
    #rbreakmin, rbreakmax = 20.0, 2.0*maxr

    if A > Amax: A = Amax
    if A < Amin: A = Amin
    #if rwid < 1.0/rwidmax: rwid = 1.0/rwidmax
    #if rwid > 1.0/rwidmin: rwid = 1.0/rwidmin
    #if bwid < 1.0/rwidmax: bwid = 1.0/bwidmax
    #if bwid > 1.0/bwidmin: bwid = 1.0/bwidmin
    if r0 < r0min: r0 = r0min
    if r0 > r0max: r0 = r0max
    if b0 < b0min: b0 = b0min
    if b0 > b0max: b0 = b0max
    #if rbreak < rbreakmin: rbreak = rbreakmin
    #if rbreak > rbreakmax: rbreak = rbreakmax

    s = 1.0    
    b0tmp = b0 + s*(r-rbreak)* pritchet( (bwid,rbreak), r)
    return A * pritchet( (rwid,r0), r) * (1.0 - pritchet( (rwid, b0tmp), b))

def residuals(p, data, r, b, rbreak, rwid, bwid):
    return data - pritchet2d(p, r, b, rbreak, rwid, bwid)

def main():
    if len(sys.argv) != 2:
        Usage()
        
    try:
        statfile = sys.argv[1]
    except:
        Usage()

    try:
        lamb = num.float(sys.argv[2])
    except:
        lamb = 5.5e-7

    thresh = 0.8

    data = num.array([], dtype=float)
    r = num.array([], dtype=float)
    b = num.array([], dtype=float)

    i = 0
    fp = open(statfile, 'r')
    for line in fp.read().splitlines():

        if ( re.search("^#", line) ):
            continue
        Rstar, aa, AU, ip, nrec, nadd, nhit, vret = (line.split())[0:8]
        frac = num.float(nrec) / num.float(nadd)
        if (frac > 1.0): frac = 1.0
        data = num.append(data, frac)
        r = num.append(r, num.float(aa))
        b = num.append(b, num.float(ip))

        # fake a negative r to help the fitter
        #r = num.append(r, -num.float(aa))
        #b = num.append(b, num.float(ip))
        #data = num.append(data, 0.0)

        # if there are multiple distances, result is meaningless ... exit
        if (i>0):
            if (AU != AU0):
                print >> sys.stderr, "ERROR: multiple AUs in stats file: %s\n" % (statfile)
                sys.exit(1)
        AU0 = AU
        i += 1

    fp.close()
    fsu = num.sqrt( num.float(AU)*1.5e11 * lamb / 2.0)
    rbreak = 0.5*fsu
               
    # try a leastsq fit
    rwid = 1.0/20.0
    bwid = 1.0 / 100.0
    p0 = [1.0, 750.0, 2000.0]
    p = leastsq(residuals, p0, (data,r,b, rbreak, rwid, bwid)) 
    A, r0, b0 =  p[0]

    if ( num.max(r)-r0 < 100 ): A = 0.0

    s = 1.0
    r_th = invPrit(rwid, r0, thresh, 0)
    b_th = invPrit(bwid, b0, thresh, 
                   s*(r_th-rbreak)*pritchet( (bwid,rbreak), r_th))
    print "%.2f %.2f %.2f %.0f %.0f %.0f %.0f %.2f %.2f %.0f" % (r_th, b_th, num.float(vret), 1.0/rwid, 1.0/bwid, r0, b0, A, s, rbreak)
        
    
if __name__ == '__main__':
    main()
    
