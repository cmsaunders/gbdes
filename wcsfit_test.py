import wcsfit
import numpy as np
print('wcsfit package:', wcsfit.__file__)

ARCSEC = np.pi/180./3600.
MILLIARCSEC = 0.001*ARCSEC
DEGREE = np.pi/180.
RESIDUAL_UNIT = MILLIARCSEC
WCS_UNIT = DEGREE


class WCSFoF:
    ## TODO: probably remove this class
    def __init__(self, fields, instruments, exposures, matchRadius=1.0, useAffinities=True, minMatches=2,
                 allowSelfMatches=False):
        
        if matchRadius < 0.0:
            raise ValueError("matchRadius must be larger than 0.0")
        self.matchRadius = matchRadius
        self.useAffinities = useAffinities
        if minMatches < 2:
            raise ValueError("minMatches must be at least 2")
        self.minMatches = minMatches
        self.allowSelfMatches = allowSelfMatches
        
        self.fields = fields
        self.instruments = instruments
        self.exposures = exposures
    
    def run(self):
        pass

class WCSFit:
    
    def __init__(self, exposures, extensions, fieldNames, fieldEpochs, fieldProjections, instruments, matches,
                 useInstrumentList=['.*'], useColorList=[], skipFile='',
                 skipExposureList=[], sysError=2, referenceSysError=2, maxError=100, fixMapsList=[]):
        # TODO: not added: `inputYaml` line 230 of WCSFit.cpp
        # add Field Info? line 247
        # add Instrument Info and Device Info (X,Y min and max)
        
        self.instruments = instruments  # make list of instrument objects here
        self.exposures = exposures  # process exposures, make list of exposure objects here
        self.extensions = extensions  # list of exposure-device combinations
        # TODO: something here to set name of reference extensions, line 337
        self.parameters = 'params here'
        
        
        # TODO: Will this do anything???
        wcsfit.loadPixelMapParser()
        
        self.fieldProjections = fieldProjections
        self.fieldNames = fieldNames
        self.fieldEpochs = fieldEpochs
        
        # Convert error parameters from I/O units to internal
        referenceSysError *= RESIDUAL_UNIT / WCS_UNIT
        sysError *= RESIDUAL_UNIT / WCS_UNIT
        maxError *= RESIDUAL_UNIT / WCS_UNIT
        if sysError > 0.0:
            astrometricCovariance = np.zeros((2, 2))
            astrometricCovariance[0, 0] = sysError*sysError
            astrometricCovariance[1, 1] = sysError*sysError
            for e in exposures:
                if e.instrument >= 0:  # Note: Instrument = 0+ -> data; Instrument = -1 -> reference catalog
                    e.addToAstrometricCovariance(astrometricCovariance)

        if referenceSysError > 0.0:
            astrometricCovariance = np.zeros((2, 2))
            astrometricCovariance[0, 0] = referenceSysError**2
            astrometricCovariance[1, 1] = referenceSysError**2
            for e in exposures:
                if ((e.instrument == wcsfit.REF_INSTRUMENT) or (e.instrument == wcsfit.PM_INSTRUMENT)):
                    e.addToAstrometricCovariance += astrometricCovariance
                    # Note that code in FitSubroutines::makePMDetection() will
                    # add this systematic error only to the Detection::invCov 2d
                    # covariance, not the full 5d PMDetection::pmInvCov calculation,
                    # since we'll assume these 5d projects (Gaia!) have treated errors well.
                    # For single-epoch reference catalogs, PM is probably the largest sys error.
        """
        # Build a preliminary set of pixel maps from the configured YAML files
        pmcInit = wcsfit.PixelMapCollection() # make this from inputYaml
        pmcInit.fixMapComponents(fixMapsList)
        # Perform various checks of maps
        # check for degeneracies, set some maps to Identity as necessary
        
        # Somehow don't need pmcInit anymore -- was it just to check for degeneracies
        # in the maps? Now make the real map collection
        self.mapCollection = wcsfit.PixelMapCollection.createMapCollection(instruments,
                                                                      exposures, extensions,
                                                                      inputYAML)
        # Add WCS for every extension, and reproject into field coordinates
        wcsfit.setupWCS(self.mapCollection)
        
        degen = self.mapCollection.degen()
        initializeOrder = degen.initializationOrder
        for extensionSet in initializationOrder:
            for extension in extensionSet:
                fitDefaulted(self.mapCollection, defaultedExtensions, instruments, exposures)
                
        for extension in extensions:
            fitDefaulted()
        
        # Check there are no more defaulted maps
        
        # Fix map components
        self.mapCollection.fixMaps(fixMapsList)
        self.mapCollection.rebuildParameterVector()
        
        print(f"# Total number of free map elements: {self.mapCollection.nFreeMaps}"
              f" with {self.mapCollection.nParams()} free parameters.")
       
        # Read in the data
        wcsfit.whoNeedsColor(extensions) 
        
        # Before reading objects, we want to set all starting WCS's to go into field coordinates.
        for extension in extensions:
            if (not extension) or (extension.exposure < 0):
                continue
            ifield = exposures.field.index
            extension.startWcs.reprojectTo(fieldProjections[ifield])
        """
        #for icat in catalogs.size: # Are there actually ever more than one?
        #    matchTable = catalogs[icat]
        #    self.matches = readMatches(matchTable, extensions, colorExtensions, skipSet, minMatches, usePM)
        self.matches = matches
        #print("# Total match count:", self.matches.size)
        
        print("Reading catalogs")
        readObjects(extensionTable, exposures, extensions, fieldProjections)

        # Now loop again over all catalogs being used to supply colors,
        # and insert colors into all the Detections they match
        print("Reading colors")
        readColors(extensionTable, colorExtensions)

        print("Purging defective detections and matches")

        # Get rid of Detections with errors too high
        purgeNoisyDetections(maxError, self.matches, exposures, extensions)
                
        print("Purging sparse matches")
        # Get rid of Matches with too few detections
        purgeSparseMatches(minMatches, self.matches)

        print("Purging out-of-range colors")
        # Get rid of Matches with color out of range (note that default color is 0).
        purgeBadColor(minColor, maxColor, self.matches)
        
        print("Reserving matches")
        if reserveFraction > 0.0:
            reserveMatches(self.matches, reserveFraction, randomNumberSeed)
        
        print("Purging underpopulated exposures")
        # Find exposures whose parameters are free but have too few
        # Detections being fit to the exposure model.
        badExposures = findUnderpopulatedExposures(minFitExposures, self.matches, exposures, extensions,
                                                   self.mapCollection)
        
        print('Purging bad exposure parameters and Detections')
        for i in badExposures:
            print("WARNING")
            self.mapCollection.freezeMap(i.first, self.matches, extensions)
        
        if purgeOutput:
            print("Purging unfittable maps")
            self.mapCollection.purgeInvalid()
        
        matchCensus(self.matches)
        
    def runWCSFit(self, chisqTolerance=0.001, verbose=False, divideInPlace=False, clipThresh=5, clipEntireMatch=False):
        #          reserveFraction, randomNumberSeed, skipFile, clipThresh, maxError, sysError, referenceSysError,
        #          freePM, pmEpoch, parallaxPrior, pmPrior, minMatches, minFitExposures, clipEntireMatch, chisqTolerance,
        #          divideInPlace, purgeOutput, inputMaps, fixMaps, useInstruments, skipExposures, outCatalog, outWcs,
        #          starCatalog, colorExposures, minColor, maxColor, verbose, inputTables, fieldNames, fieldProjections,
        #          fieldEpochs):

        # Now do the refitting
        print("Begin fitting process")
        
        ca = wcsfit.coordAlign(self.mapCollection, self.matches)
        
        nclip = 0
        oldthresh = 0.0
        coarsePasses = True
        ca.setRelTolerance(10 * chisqTolerance)
        minimumImprovement = 0.02
        
        while (coarsePasses or (nclip > 0)):
            # Report number of active Matches / Detections in each iteration:
            
            mcount = 0
            dcount = 0
            ca.count(mcount, dcount, False, 2)
            maxdev = 0.0
            dof = 0
            chi = ca.chisqDOF(dof, maxdev, False)
            if (verbose >= 1):
                print(f"Fitting {mcount} matches with {dcount} detections chisq {chi} / {dof} dof,  maxdev "
                      f"{maxdev} sigma")            

            # Do the fit here!!
            chisq = ca.fitOnce(verbose >= 1, divideInPlace)  # save space if selected
            # Note that fitOnce() remaps *all* the matches, including reserved ones.

            dof, max = ca.chisqDOF(False)  # Exclude reserved Matches
            thresh = np.sqrt(chisq/dof) * clipThresh  # ??? change dof to expectedChisq?
            if (verbose >= 1):
                print(f"After iteration: chisq {chisq} / {dof} dof, max deviation {max} new clip threshold"
                      f"at: {thresh} sigma")

            if (thresh >= max or (oldthresh > 0.0 and (1-thresh/oldthresh) < minimumImprovement)):
                # Sigma clipping is no longer doing much.  Quit if we are at full precision,
                # else require full target precision and initiate another pass.
                if (coarsePasses):
                    coarsePasses = False
                    ca.setRelTolerance(chisqTolerance)
                    print("Starting strict tolerance passes")
                    if (clipEntireMatch & verbose >= 1):
                        print("-->clipping full matches")
                    oldthresh = thresh
                    nclip = ca.sigmaClip(thresh, False, clipEntireMatch & ~coarsePasses, verbose >= 1)
                    if (verbose >= 0):
                        print(f"# Clipped {nclip} matches ")
                    continue
                else:
                    # Done!
                    break

            oldthresh = thresh
            # Clip entire matches on final passes if clipEntireMatch=true
            nclip = ca.sigmaClip(thresh, False, clipEntireMatch & ~coarsePasses, verbose >= 1)
            if ((nclip == 0) & coarsePasses):
                # Nothing being clipped; tighten tolerances and re-fit
                coarsePasses = False
                ca.setRelTolerance(chisqTolerance)
                print("Starting strict tolerance passes")
                if (clipEntireMatch & verbose >= 1):
                    print("-->clipping full matches")
                continue
            
            if verbose >= 0:
                print(f"# Clipped {nclip} matches")
            
        wcsfit.clipReserved(ca, clipThresh, minimumImprovement, False, verbose >= 1)
        
        extensionProjections = np.zeros(len(self.extensions))
        for i, extension in enumerate(self.extensions):
            if not extension:
                continue
            iExposure = extension.exposure
            if (iExposure < 0) or (not self.exposures[iExposure]):
                continue
            iField = self.exposures[iExposure].field
            extensionProjections[i] = self.fieldProjections[iField]
            
        wcsfit.reportStatistics(self.matches, self.exposures, self.extensions)
        
        return self.mapCollection, extensionProjections
            
        
if __name__ == "__main__":
    
    w = WCSFit()
    maps, projections = w.runWCSFit()