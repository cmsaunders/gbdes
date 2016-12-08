#!/usr/bin/env python
'''
Read FOF file and generate YAML-format specifications for differential chromatic refraction terms
for each exposure (DCR).  Also calculates an improved airmass value, placing this and parallactic
angle into columns of the EXPOSURE table in the FOF file.
'''
import numpy as np
import yaml
from gmbpy import dmsToDegrees

def parallactic(dec, ha, lat=-30.1716):
    '''Function will calculate airmass and
    parallactic angle (in radians) given
    input of source declination, hour angle,
    and observatory latitude (in degrees)
    '''
    # First get cosine of zenith angle
    dtor = np.pi / 180.
    d = dec*dtor
    if (ha>180.):
        h = (ha-360.)*dtor  # Make negative HA instead of large
    else:
        h = ha*dtor
    l = lat*dtor
    cosz = np.sin(l)*np.sin(d) + np.cos(l)*np.cos(d)*np.cos(h)
    sinz = np.sqrt(1-cosz*cosz)  # zenith angle always 0-180
    airmass = 1./cosz
    # Now the parallactic angle
    if sinz<=0:
        # at (anti-)zenith already, p is undefined, return 0
        return airmass,0.
    cosp = (np.sin(l)*np.cos(d)-np.cos(l)*np.sin(d)*np.cos(h))/sinz
    if np.abs(cosp)>1.:
        cosp = np.sign(cosp)
    p = np.arccos(cosp)* np.sign(h)
    return airmass,p

dcrConstant = {'g':45.0, 'r':8.4}  # DCR amplitude in mas/mag/tan(z)
latitude = -30.1716  # Observatory latitude, degrees
referenceColor = 0.61  # zeropoint for color terms

if __name__=='__main__':
    fits = sys.argv(1) # ?? use argparse
    dcrfile = sys.argv(2)

    pixmaps = {}
    
#    ff = pf.open(fits,'update')
    ff = pf.open(fits)

    # Get bands of all instruments, assign band to exposures
    bands = {}
    for hdu in ff[1:]:
        if hdu.header['EXTNAME']=='Instrument':
            bands[hdu.header['NUMBER'])] = hdu.header['BAND']

    exposures = ff['exposures'].data
    extns = ff['extensions'].data
    
    # Read exposure table
    for i,expo in enumerate(exposures['name']):
        iinst = exposures['instrumentnumber'][i]
        if iinst < 0:
            # Not a real instrument, no changes here
            continue
        name = "".join(expo)
        # Find all extensions with this exposure
        used = extns['exposure']==i
        # Get HA from them, do they agree?
        ha = np.unique(extensions['HA'][used])
        if len(ha) > 1:
            print "Disagreeing HA's for exposure",name,ha
        ha = dmsToDegree(ha) * 15.
        # get declination, derive parallactic angle
        dec = exposures['dec'][i]
        airmass, p = parallactic(dec, ha, lat = latitude)
        # replace airmass in table
        exposures['airmass'][i] = airmass
        # Write DCR to pixmaps
        b = bands[iinst]
        if b in dcrConstant.keys():
            # not a band that needs DCR
            pixmaps[name+'/dcr'] = {'Type':'Identity'}
        else:
            # Convert DCR amplitude from mas to degrees
            ampl = dcrConstant[b] / (3600.*1000.)
            dcry = ampl * np.cos(p)*np.sqrt(airmass*airmass-1)
            dcrx = ampl * np.sin(p)*np.sqrt(airmass*airmass-1)
            print name,dcrx,dcry
            # Specify the distortion
            pixmap = {'Type':'Color',
                      'Reference':referenceColor,
                      'Function':{'Type':'Constant',
                                  'Parameters':[dcrx,dcry]}}
            pixmaps[name+'/dcr'] = pixmap

    # write yaml
    fout = open(yamlFile,'w')
    yaml.dump(pixmaps,fout)
    fout.close()
    ff.close()
    
