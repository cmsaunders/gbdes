
from lsst.daf.butler import Butler
import wcsfit
import wcsfit_test
import numpy as np
import sip_tpv

#butler = Butler('/repo/main', collections="HSC/runs/RC2/w_2021_30/DM-31182")
butler = Butler('/repo/main', collections="HSC/runs/RC2/w_2021_34/DM-31524")

"""
print('test package:', wcsfit_test.__file__)
inputTables = '/home/csaunder/stack_projects/gbdes_tests/cosmos_pdr2_field/match_newgaia.fof'
"""
# GBDES defaults:
reserveFraction = 0.05
randomNumberSeed = 0
skipFile = ''
clipThresh = 5
maxError = 100
sysError = 2
referenceSysError = 2
freePM = True
pmEpoch = 2015.5
parallaxPrior = 10
pmPrior = 100
minMatches = 2
minFitExposures = 200
clipEntireMatch = False
chisqTolerance = 0.001
divideInPlace = False
purgeOutput = False
inputMaps = '/home/csaunder/stack_projects/gbdes_tests/polyExposure.astro'
fixMaps = ''
useInstruments = '.*'
skipExposures = ''
outCatalog = 'wcscat.fits'
outWcs = 'wcsfit.wcs'
starCatalog = 'starcat.fits'
colorExposures = ''
minColor = -10
maxColor = 10
verbose = 1

pvrange, tpvx, tpvy = sip_tpv.pvsiputils.sym_tpvexprs()


def pv_matrix_from_sip(sipx, sipy):
    # This function is a combination of sip_tpv.add_pv_keywords and ReadPV in gbdes/src/subs/TPVMap.cpp
    coeffs_size = max(pvrange) + 1

    pv1_vals = np.zeros(coeffs_size)
    pv2_vals = np.zeros(coeffs_size)
    for p in pvrange:
        # makes PV1_{p}, PV2_{p}
        pv1_vals[p] = float(sip_tpv.pvsiputils.calcpv(pvrange, 1, p, sipx, sipy, tpvx, tpvy).evalf())
        pv2_vals[p] = float(sip_tpv.pvsiputils.calcpv(pvrange, 2, p, sipx, sipy, tpvx, tpvy).evalf())
    
    coeffs_size = np.flatnonzero(pv1_vals != 0).max() + 1
    order = 1
    while ((order + 1) * (order + 2) / 2) < coeffs_size:
        order += 1

    pv1 = np.zeros((order + 1, order + 1))
    pv2 = np.zeros((order + 1, order + 1))
    i = 0
    j = 0
    for k in range(coeffs_size):
        pv1[i, j] = pv1_vals[k]
        pv2[i, j] = pv2_vals[k]
        if i == 0:
            i = j + 1
            j = 0
        else:
            i -= 1
            j += 1

    if 1 not in pvrange:
        pv1[1, 0] = 1.0
        pv2[1, 0] = 1.0

    pv1 = pv1.transpose()
    pv2 = pv2.transpose()
    return pv1, pv2

POLYSTEP = 1.0/3600.0
def convert_sip_to_tpv(butler_wcs, name=""):
    # This function is an adaptation of readTPV in gbdes/src/subs/TPVMap.cpp
    print('Start sip to tpv')
    # Format copied from sip_tpv.sip_pv function
    cd = butler_wcs.getCdMatrix()
    fits_metadata = butler_wcs.getFitsMetadata()

    
    if not ((fits_metadata.get('CTYPE1') == 'RA---TAN-SIP') and
            (fits_metadata.get('CTYPE2') == 'DEC--TAN-SIP')):
        raise ValueError(f'CTYPES {fits_metadata.get("CTYPE1")} and {fits_metadata.get("CTYPE2")}'
                         'do not match SIP convention')
    
    plin = np.zeros(6)
    crpix1 = fits_metadata.get('CRPIX1')
    crpix2 = fits_metadata.get('CRPIX2')
    plin[0] = - cd[0, 0] * crpix1 - cd[0, 1] * crpix2
    plin[1] = cd[0, 0]
    plin[2] = cd[0, 1]
    plin[3] = -cd[1, 0] * crpix1 - cd[1, 1] * crpix2
    plin[4] = cd[1, 0]
    plin[5] = cd[1, 1]
    linearMap = wcsfit.LinearMap(plin)

    pole = wcsfit.SphericalICRS(fits_metadata.get('CRVAL1') * wcsfit.DEGREE,
                                fits_metadata.get('CRVAL2') * wcsfit.DEGREE)
    orientIn = wcsfit.Orientation()
    orientIn.set(pole, 0.0)
    tp = wcsfit.Gnomonic(orientIn.getPole(), orientIn)

    a_order = int(fits_metadata.get('A_ORDER', 0))
    b_order = int(fits_metadata.get('B_ORDER', 0))
    ac = np.matrix(np.zeros((a_order+1, a_order+1), dtype=np.float64))
    bc = np.matrix(np.zeros((b_order+1, b_order+1), dtype=np.float64))
    for m in range(a_order+1):
        for n in range(0, a_order+1-m):
            ac[m, n] = fits_metadata.get('A_%d_%d' % (m, n), 0.0)
    for m in range(b_order+1):
        for n in range(0, b_order+1-m):
            bc[m, n] = fits_metadata.get('B_%d_%d' % (m, n), 0.0)
    sipx, sipy = sip_tpv.pvsiputils.real_sipexprs(cd, ac, bc)
    pv1, pv2 = pv_matrix_from_sip(sipx, sipy)
    p1 = wcsfit.Poly2d(pv1)
    p2 = wcsfit.Poly2d(pv2)

    pmlist = [linearMap]
    polyName = name + "_pv"
    if p1.nCoeffs() > 1:
        # Obtained a valid x polynomial.
        if p2.nCoeffs() > 1:
            # Also have valid y polynomial.  Make the map:
            pv = wcsfit.PolyMap(p1, p2, polyName, wcsfit.Bounds(-1.0, 1.0, -1.0, 1.0), POLYSTEP)
            
        else:
            # Did not find PV2's.  Install default:
            p2Identity = wcsfit.Poly2d(1)
            coeffs = p2Identity.getC()
            coeffs[p2Identity.vectorIndex(0, 1)] = 1.0
            p2Identity.setC(coeffs)
            pv = wcsfit.PolyMap(p1, p2Identity, polyName, wcsfit.Bounds(-1.0, 1.0, -1.0, 1.0), POLYSTEP)
        pmlist.append(pv)
    else:
        # Did not obtain any PV1's.  If there are PV2's, install identity map for x coeff:
        if p2.nCoeffs() > 1:
            p1Identity = wcsfit.Poly2d(1)
            coeffs = p1Identity.getC()
            coeffs[p1Identity.vectorIndex(1, 0)] = 1.0
            p1Identity.setC(coeffs)
            pv = wcsfit.PolyMap(p1Identity, p2, polyName, wcsfit.Bounds(-1.0, 1.0, -1.0, 1.0), POLYSTEP)
            pmlist.append(pv)

    #Create a SubMap that owns a duplicate of these one or two maps (no sharing):
    sm = wcsfit.SubMap(pmlist, name, False)
    
    # Return the Wcs, which again makes its own copy of the maps (no sharing):
    return wcsfit.Wcs(sm, tp, name, wcsfit.DEGREE, False)


def readExposures_wcsfof(refs):
    wcsList = []
    for v, visitSummaryRef in enumerate(refs):
        visitSummary = butler.get(visitSummaryRef)
        
        for row in visitSummary[:2]:
            print(f'Processing visit {row["visit"]}, detector {row["id"]}')
            calexp = butler.get('calexp', visit=row['visit'], detector=row['id'])
            butler_wcs = calexp.getWcs()
            wcs = convert_sip_to_tpv(butler_wcs)
            wcsList.append(wcs)
    return wcsList

def readExposuresExtensions(refs, fieldNumber=0, instrumentNumber=0, isReference=False):
    # TODO: add ExposureColorPriorities, ColorExtension / decide if necessary
    
    extensions = []
    exposures = []
    
    for v, visitSummaryRef in enumerate(refs):
        print(visitSummaryRef)        
        visitSummary = butler.get(visitSummaryRef)
        visInfo = visitSummary[0].getVisitInfo()
        
        # TODO: ra dec I think should be in radians based in line 482 in gbdes/src/FitSubroutines.cpp
        # -- need to double check
        ra = visInfo.getBoresightRaDec().getRa().asRadians()
        dec = visInfo.getBoresightRaDec().getDec().asRadians()
        gn = wcsfit.Gnomonic(wcsfit.Orientation(wcsfit.SphericalICRS(ra, dec)))
        # TODO: probably want a different id/name in future:
        exposure = wcsfit.Exposure(str(visInfo.getId()), gn)
        exposure.field = fieldNumber
        exposure.instrument = instrumentNumber
        airmass = visInfo.boresightAirmass
        exptime = visInfo.exposureTime
        exposure.airmass = airmass
        exposure.exptime = exptime
        exposure.mjd = visInfo.date.get(visInfo.date.DateSystem.MJD)
        if False:
            exposure.pmEpoch
        if False:
            exposure.apcorr
        # TODO: add astrometric covariance
        exposures.append(exposure)

        for row in visitSummary[:2]:
            extension = wcsfit.Extension()
            extension.exposure = v
            extension.device = row['id']
            extension.airmass = airmass
            extension.magshift = 2.5 * np.log10(exptime)
            if False:
                extension.apcorr
            # TODO: add WCS here!
            extensions.append(extension)
    return exposures, extensions

def readInstruments(filtername):
    
    instruments = []
    for i in list(range(9)) + list(range(10, 104)):
        instrument = wcsfit.Instrument()    
        instrument.band = filtername
        # TODO: replace with real bounds (~25 pixel clipping around edges?)
        instrument.addDevice(str(i), wcsfit.Bounds(0, 2047, 0, 4175))
        instruments.append(instrument)
        
    return instruments


def getMatchArray(catalog, allowSelfMatches=False):
    # Probably move this back to c++
    sequence = []
    extn = []
    obj = []
    matches = 0
    # Now loop through matches in this catalog
    print("Start loop")
    
    catalogVector = catalog.toVector()
    import ipdb
    ipdb.set_trace()
    for pt in catalogVector:
        print("get point")
        point = pt.toVector()
        # Skip any Match that is below minimum match size
        if len(point) < minMatches:
            print("not enough matches")
            continue
        inExposures = []
        selfMatch = False
        if not allowSelfMatches:
            print("in aSM")
            print(len(point))
            for detection in point:
                print(detection)
                if detection.exposureNumber in inExposures:
                    selfMatch = True
                    break
                inExposures.append(detection.exposureNumber)
            if selfMatch:
                continue
        print("now read")
        #Add elements of this match to the output vectors
        seq = 0
        matches += 1
        for detection in point:
            sequence.append(seq)
            extn.append(detection.extensionNumber)
            obj.append(detection.objectNumber)
            seq += 1
        1/0

fields = [9813]
bands = ['r']
usePM = False

skyMap = butler.get('skyMap')

fieldNames = []
fieldProjections = []
fieldObjects = []
for f, field in enumerate(fields):
    #fieldNames.append(wcsfit.spaceReplace(str(field)))
    fieldNames.append(str(field))
    # TODO: does tractinfo have a center?
    
    tractInfo = skyMap.generateTract(field)    
    skyOrigin = tractInfo.getWcs().getSkyOrigin()
    ra = skyOrigin.getRa().asRadians()
    dec = skyOrigin.getDec().asRadians()
    fieldProjection = wcsfit.Gnomonic(wcsfit.Orientation(wcsfit.SphericalICRS(ra, dec)))
    fieldProjections.append(fieldProjection)
    
    #expRefList = list(set(butler.registry.queryDatasets('calexp', dataId={"band": "r", "tract": field, 
    #                                                                      "patch": 25, 
    #                                                                      "skymap": "hsc_rings_v1"})))
    visitSummaryRefs = list(set(butler.registry.queryDatasets('visitSummary', dataId={"tract": field})))
    
    #exposures = readExposures(expRefList[:3], fieldNumber=f, instrumentNumber=0)
    
    ## WCSFoF setup:
    
    fieldObject = wcsfit.Field()
    fieldObject.name = str(field)
    fieldObject.projection = fieldProjection
    # TODO: matchProjection does not need to be a property of wcsfof
    fieldObject.matchRadius = 1
    fieldObject.extent = 3.0  # TODO: set something more intelligent here
    fieldObjects.append(fieldObject)
    
    instObject = wcsfit.Instr()
    instObject.name = 'Hyper Suprime-Cam'
    for d in list(range(9)) + list(range(10, 104)):
        device = wcsfit.Device()
        device.name = str(d)
        instObject.append(device)
    instrumentObjects = [instObject]
    
    exposureObjects = []
    ## TODO: reconsider whether we want to load with calexps of visitSummary tables
    for visitSummaryRef in visitSummaryRefs[:2]:
        visitSummary = butler.get(visitSummaryRef)
        visInfo = visitSummary[0].getVisitInfo()
        expo = wcsfit.Expo()
        ra = visInfo.getBoresightRaDec().getRa().asRadians()
        dec = visInfo.getBoresightRaDec().getDec().asRadians()
        expo.name = str(visInfo.getId())
        expo.field = f
        expo.instrument = 0
        expo.pointing = wcsfit.SphericalICRS(ra, dec)
        exposureObjects.append(expo)
    
    #wcsfof = wcsfit_test.WCSFoF(fieldObjects, instrumentObjects, exposureObjects)
    wcsfof = wcsfit.FoFClass()
    wcsfof.fields = fieldObjects
    wcsfof.instruments = instrumentObjects
    wcsfof.exposures = exposureObjects
    #import ipdb
    #ipdb.set_trace()
    
    iextn = 0
    for v, visitSummaryRef in enumerate(visitSummaryRefs[:2]):
        visitSummary = butler.get(visitSummaryRef)
        exposureNumber = v
        instrumentNumber = 0
        fieldNumber = f
        for row in visitSummary[:2]:
            thisAffinity = row['band'] # TODO: not sure if this is what is wanted for AFFINITY
            deviceNumber = row['id']
            srcCatalog = butler.get('src', dataId={"visit": row['visit'], "detector": row['id']})
            vx = srcCatalog['base_SdssCentroid_x']
            vy = srcCatalog['base_SdssCentroid_y']
            vid = np.arange(len(vx)) # TODO: anything better here?
            isStar = np.ones(len(vx))
            calexp = butler.get('calexp', dataId={"visit": row['visit'], "detector": row['id']})
            wcs = convert_sip_to_tpv(calexp.getWcs())
            # TODO: need reprojectWCS to fields[field].projection here
            wcsfof.addCatalog(wcs, thisAffinity, exposureNumber, fieldNumber, instrumentNumber, deviceNumber,
                              iextn, isStar, vx, vy, vid)
            print('Added catalog')
            iextn += 1  # TODO: replace with unique extn + visit number?
            
    # TODO: write out matches here, check against baseline result
    print(wcsfof.fields[0].catalogs.keys())
    catalog = wcsfof.fields[0].catalogs['STELLAR']
    print("getting matchArray")
    seqs = []
    extns = []
    objs = []
    import ipdb
    ipdb.set_trace()
    wcsfof.writeMatches("testcat.fits")
    wcsfof.sortMatches(0)
    
    #seqs, extns, objs = getMatchArray(catalog)
    
    ## TODO: don't need next two lines anymore?
    #wcs_list = readExposures_wcsfof(visitSummaryRefs[:2])
    # TODO wcs - reprojectTo(fields[field].projection)
    
    ## WCSFit setup:
    print("reading exposures")
    exposures, extensions = readExposuresExtensions(visitSummaryRefs[:2], fieldNumber=f, instrumentNumber=0)
    skipSet = wcsfit.ExtensionObjectSet("")
    colorExtensions = [wcsfit.ColorExtension() for e in extensions]
    matches = wcsfit.MCat()
    wcsfit.readMatches(seqs, extns, objs, matches, extensions, colorExtensions, skipSet, minMatches, usePM)

    fieldProjections = [150, 2.5]
    fieldEpochs = [2015.5]

    instruments = readInstruments(bands[0])

    fitter = wcsfit_test.WCSFit(exposures, extensions, fieldNames, fieldEpochs, fieldProjections, 
                                instruments, matches, fixMapsList=[1, 2])
    print(fitter.exposures)
