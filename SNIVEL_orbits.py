#!/usr/bin/env python
import numpy
import georinex as gr
import math
from scipy import optimize
from scipy.interpolate import lagrange
import datetime
####################################################
#Constants
####################################################
mu = 3.986004418e14 #G*Mearth
omeg_E = 7.2921151467e-5 #Earth's rotation rate
c = 299792458.0 #Speed of light

def bcorbit(nav, svn, tsv, x0, y0, z0):
    ####################################################
    #Loading Nav File
    ####################################################
    tk_all = tsv - nav['Toe'].sel(sv=svn)
    tk_all_min = numpy.nanmin(numpy.abs(tk_all))
    [tk_ind] = numpy.where(tk_all_min == numpy.abs(tk_all))
    dind = tk_ind[0] 
    tk = tk_all[dind]
    ####################################################
    #Different Satellite parameters
    ####################################################
    toe = nav['Toe'].sel(sv=svn).data[dind]
    a = math.pow(nav['sqrtA'].sel(sv=svn).data[dind],2)
    e = nav['Eccentricity'].sel(sv=svn).data[dind]
    m0 = nav['M0'].sel(sv=svn).data[dind]
    w = nav['omega'].sel(sv=svn).data[dind]
    i0 = nav['Io'].sel(sv=svn).data[dind]
    omg0 = nav['Omega0'].sel(sv=svn).data[dind]
    dn = nav['DeltaN'].sel(sv=svn).data[dind]
    idot = nav['IDOT'].sel(sv=svn).data[dind]
    odot = nav['OmegaDot'].sel(sv=svn).data[dind]
    iode = nav['IODE'].sel(sv=svn).data[dind]
    af0 = nav['SVclockBias'].sel(sv=svn).data[dind]
    af1 = nav['SVclockDrift'].sel(sv=svn).data[dind]
    af2 = nav['SVclockDriftRate'].sel(sv=svn).data[dind]
    toc = nav['TransTime'].sel(sv=svn).data[dind]


    cuc = nav['Cuc'].sel(sv=svn).data[dind]
    cus = nav['Cus'].sel(sv=svn).data[dind]
    cic = nav['Cic'].sel(sv=svn).data[dind]
    cis = nav['Cis'].sel(sv=svn).data[dind]
    crc = nav['Crc'].sel(sv=svn).data[dind]
    crs = nav['Crs'].sel(sv=svn).data[dind]
    ####################################################
    #Computations...
    ####################################################
    dt = 0
    ti = tsv
    while dt < 3:
        #mean angular motion
        n0 = math.sqrt(mu/math.pow(a,3))
        n = n0+dn
        #Mean anomaly
        M = m0 + tk*n
        ####################################################
        #Keplers equation for eccentric anomaly estimate
        ####################################################
        def f(x):
            return x-e*math.sin(x)-M
        E=optimize.newton(f,1)
        ####################################################
        #Time Correction
        ####################################################
        F = -2*math.sqrt(mu)/math.pow(c,2)
        delta_tr = F*e*math.sqrt(a)*math.sin(E)
        delta_tsv = af0+af1*(ti-toe)+delta_tr
        prange_corr=delta_tsv*c
        t = ti-delta_tsv
        
        tk = t-toe
        M = m0+n*tk
        def f(x):
            return x-e*math.sin(x)-M
        E=optimize.newton(f,1)    
        ####################################################
        #True anomaly
        ####################################################
        v = math.atan2((math.sqrt(1-math.pow(e,2))*math.sin(E)),(math.cos(E)-e))
        ####################################################
        #Central body distance
        ####################################################
        phi = v+w
        rc_correction = crc*math.cos(2*phi) + crs*math.sin(2*phi)
        rc = a*(1-e*math.cos(E)) + rc_correction
        ####################################################
        #Argument of latitude
        ####################################################
        u_correction = cuc*math.cos(2*phi) + cus*math.sin(2*phi)
        u = phi + u_correction
        ####################################################
        #Inclination
        ####################################################
        i_correction = cic*math.cos(2*phi) + cis*math.sin(2*phi)
        i = i0 + idot*tk + i_correction
        ####################################################
        #Longitude of ascending node (note corrected formula)
        ####################################################
        Omega = omg0 + (odot - omeg_E)*tk - omeg_E*toe
        ####################################################
        #Positions
        ####################################################
        x = rc*(math.cos(u)*math.cos(Omega) - math.sin(u)*math.cos(i)*math.sin(Omega))
        y = rc*(math.cos(u)*math.sin(Omega) + math.sin(u)*math.cos(i)*math.cos(Omega))
        z = rc*(math.sin(u)*math.sin(i))

        rho = math.sqrt(math.pow(x-x0,2)+math.pow(y-y0,2)+math.pow(z-z0,2))
        ts = tsv - rho/c
        dt = dt+1
        ti = ts

    return (x,y,z,delta_tr,rho,delta_tsv)



####################################################
#Precise orbits - reader and interpolator
####################################################
#This code reads in an sp3 file and outputs matrices of the satellite positions, times, and
#PRN numbers
def readsp3(sp3file):
    k=0
    with open(sp3file, 'rt') as g:
        rows = (line.split() for line in g)
        for grow in rows:
            if (k == 0):
                nt =  int(grow[6])
            if (k == 2):
                numsat = int(grow[1])
            k=k+1

    xpos = numpy.nan*numpy.ones([nt,numsat])
    ypos = numpy.nan*numpy.ones([nt,numsat])
    zpos = numpy.nan*numpy.ones([nt,numsat])
    satclock = numpy.nan*numpy.ones([nt,numsat])
    PRN = numpy.nan*numpy.ones([1,numsat])
    gpst = numpy.nan*numpy.ones([nt,1])

    k=0
    n=0
    t=0
    with open(sp3file, 'rt') as g:
        rows = (line.split() for line in g)
        for grow in rows:
            if (k == 0):
                if (grow[0] == "*"):
                    year = grow[1]
                    month = grow[2]
                    day = grow[3]
                    hour = grow[4]
                    minute = grow[5]
                    second = grow[6]
                    dtime64 = year + '-' + month.zfill(2) + '-' + day.zfill(2) + 'T' + hour.zfill(2) + ':' + minute.zfill(2) + ':' + "{0:011.8f}".format(float(second))
                    gps_time = (numpy.datetime64(dtime64) - numpy.datetime64('1980-01-06T00:00:00'))/ numpy.timedelta64(1, 's')
                    gps_week = int(gps_time/604800)
                    gps_sow = gps_time - gps_week*604800
                    gpst[t,0] = gps_time
                    k = 1
                    t=t+1
            else:
                if (n < numsat-1):
                    sat = grow[0]
                    PRN[0,n] = int(sat.replace("PG",""))
                    xpos[t-1,n] = grow[1]
                    ypos[t-1,n] = grow[2]
                    zpos[t-1,n] = grow[3]
                    satclock[t-1,n] = grow[4]
                    n=n+1
                elif (n == numsat-1):
                    sat = grow[0]
                    PRN[0,n] = int(sat.replace("PG",""))
                    xpos[t-1,n] = grow[1]
                    ypos[t-1,n] = grow[2]
                    zpos[t-1,n] = grow[3]
                    satclock[t-1,n] = grow[4]
                    n=0
                    k=0
    return(PRN, gpst, xpos, ypos, zpos, satclock)

def readantex(antexfile):
    k=0
    with open(antexfile, 'rt') as g:
        rows = (line.split() for line in g)
        for grow in rows:
            if (grow[0] == "START"):
                k=k+1
    SVN = numpy.nan*numpy.ones([k,1])
    T1 = numpy.nan*numpy.ones([k,1])
    T2 = 1e12*numpy.ones([k,1])
    DX = numpy.nan*numpy.ones([k,1])
    DY = numpy.nan*numpy.ones([k,1])
    DZ = numpy.nan*numpy.ones([k,1])
    
    ksvn=0
    kt1=0
    kt2=0
    kl1=0

    with open(antexfile, 'rt') as g:
        rows = (line.split() for line in g)
        for grow in rows:
            if (grow[0] == "BLOCK"):
                sat = grow[2]

                SVN[ksvn,0] = int(sat.replace("G",""))
            if (len(grow) >= 7):
                if (grow[6] == "VALID" and grow[7] == "FROM"):
                    year = grow[0]
                    month = grow[1]
                    day = grow[2]
                    hour = grow[3]
                    minute = grow[4]
                    second = grow[5]
                    dtime64 = year + '-' + month.zfill(2) + '-' + day.zfill(2) + 'T' + hour.zfill(2) + ':' + minute.zfill(2) + ':' + "{0:011.8f}".format(float(second))
                    T1[kt1,0] = float((numpy.datetime64(dtime64) - numpy.datetime64('1980-01-06T00:00:00'))/ numpy.timedelta64(1, 's'))
                if (grow[6] == "VALID" and grow[7] == "UNTIL"):
                    year = grow[0]
                    month = grow[1]
                    day = grow[2]
                    hour = grow[3]
                    minute = grow[4]
                    second = grow[5]
                    dtime64 = year + '-' + month.zfill(2) + '-' + day.zfill(2) + 'T' + hour.zfill(2) + ':' + minute.zfill(2) + ':' + "{0:011.8f}".format(float(second))
                    T2[kt2,0] = float((numpy.datetime64(dtime64) - numpy.datetime64('1980-01-06T00:00:00'))/ numpy.timedelta64(1, 's'))
            if (len(grow) > 3):
                if (grow[3] == "NORTH"):
                    DX[kl1,0] = grow[0]
                    DY[kl1,0] = grow[1]
                    DZ[kl1,0] = grow[2]

            if (len(grow) == 3):
                if (grow[0] == "END" and grow[2] == "ANTENNA"):
                    ksvn=ksvn+1
                    kt1=kt1+1
                    kt2=kt2+1
                    kl1=kl1+1
    return (SVN,T1,T2,DX,DY,DZ)



#Interpolate the precise orbits using 11-point Lagrangian interpolation
#Relativistic clock correction from broadcast orbits is used to correct transmission time
#as well as precise satellite clocks
def sp3interp(tt, satnum, PRN,  gpst, xpos, ypos, zpos, satclock, delta_tr, x0, y0, z0, SVN_ant, T1_ant, T2_ant, DX_ant, DY_ant, DZ_ant):
    a1 = numpy.argmin(abs(tt-gpst))
    b1 = numpy.where(int(satnum.replace("G","")) == PRN)[1]
    bant = (int(satnum.replace("G","")) == SVN_ant) & (tt >= T1_ant) & (tt < T2_ant)
    bant1 = numpy.where(bant == True)[0]

    dxant = DX_ant[bant1,0]/1000
    dyant = DY_ant[bant1,0]/1000
    dzant = DZ_ant[bant1,0]/1000

    if (len(b1) > 0):
        xinput = xpos[a1-5:a1+6,b1]
        yinput = ypos[a1-5:a1+6,b1]
        zinput = zpos[a1-5:a1+6,b1]
        cinput = satclock[a1-5:a1+6,b1]

        aclock = numpy.amax(cinput/1e-6)

        if (aclock > 99999):
            xnew = 999999.9
            ynew = 999999.9
            zpred = 999999.9
            cpred = 999999.9
            rclock = 999999.9
            rpath = 999999.9
            rho = 999999.9
        else:
            #light time equation correction for satellite broadcast time
            ti = tt
            dt = 0
            while dt < 4:
                tinput = gpst[a1-5:a1+6,0]-ti-delta_tr
                #Interpolate clocks first to estimate time of signal transmission
                pc = lagrange(tinput,cinput)
                cpred = numpy.polyval(pc,0)
                tinput = tinput-cpred
                px = lagrange(tinput,xinput)
                py = lagrange(tinput,yinput)
                pz = lagrange(tinput,zinput)
                xpred = numpy.polyval(px,0)+dxant
                ypred = numpy.polyval(py,0)+dyant
                zpred = numpy.polyval(pz,0)+dzant
                r = math.sqrt(math.pow(xpred,2)+math.pow(ypred,2)+math.pow(zpred,2))
                phi = -omeg_E*r/c
                xnew = math.cos(phi)*xpred-math.sin(phi)*ypred
                ynew = math.sin(phi)*xpred+math.cos(phi)*ypred
                rclock = (cpred+delta_tr)*c #relativistic clock correction with satellite clock correction
                rho = math.sqrt(math.pow(xnew-x0,2)+math.pow(ynew-y0,2)+math.pow(zpred-z0,2))
                rho_sta = math.sqrt(math.pow(x0,2)+math.pow(y0,2)+math.pow(z0,2))
                rpath = 2*mu/c/c*math.log((r+rho_sta+rho)/(r+rho_sta-rho)) #relativistic path correction
                ts = tt - rho/c
                dt = dt+1
                ti = ts
    else:
        xnew = 999999.9
        ynew = 999999.9
        zpred = 999999.9
        cpred = 999999.9
        rclock = 999999.9
        rpath = 999999.9
        rho = 999999.9
        

    return(xnew, ynew, zpred, cpred, rclock, rho, rpath)

