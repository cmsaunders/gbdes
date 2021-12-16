#ifndef WCSFIT_FIT_H
#define WCSFIT_FIT_H

#include <map>
#include <algorithm>

#include "Std.h"
#include "Astrometry.h"
#include "FitsTable.h"
#include "StringStuff.h"
#include "Pset.h"
#include "PixelMapCollection.h"
#include "YAMLCollector.h"
#include "Match.h"
#include "Instrument.h"

#include "FitSubroutines.h"
#include "WcsSubs.h"
#include "MapDegeneracies.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace astrometry;
using namespace stringstuff;
using img::FTable;

class FitClass {
  public:
    FitClass();
    //FitClass(string inputMaps);
    
    double reserveFraction;
    int randomNumberSeed;
    
    double clipThresh;
    double maxError;
    //double sysError;
    //double referenceSysError;
    //double parallaxPrior;
    //double pmPrior;
    
    int minMatches;
    int minFitExposures;
    
    bool clipEntireMatch;
    double chisqTolerance;
    bool divideInPlace;
    bool purgeOutput;
    
    double minColor;
    double maxColor;

    bool usePM;
    
    int verbose;


    // This is list of regexes of PixelMap names (or instrument names) that should
    // have their parameters held fixed.
    list<string> fixMapList;
    
    // The list of instruments that we will be matching together in this run:
    // have their parameters held fixed.
    //list<string> useInstrumentList;

    // Positional accuracy demanded of numerical solutions for inversion of 
    // pixel maps: 
    //const double worldTolerance = 0.1*MILLIARCSEC/WCS_UNIT;
    // Fractional reduction in RMS required to continue sigma-clipping:
    double minimumImprovement;
    

    // All we care about fields are names and orientations:
    //NameIndex fieldNames;
    //vector<SphericalCoords*> fieldProjections;
    //vector<double> fieldEpochs;
    Fields fields;
    
    vector<unique_ptr<Instrument>> instruments;
    
    // The table of exposures
    vector<unique_ptr<Exposure>> exposures;
    vector<shared_ptr<Exposure>> SPexposures;
    
    // Extension tables:
    vector<unique_ptr<ColorExtension>> colorExtensions;
    vector<unique_ptr<Extension>> extensions;
    vector<shared_ptr<Extension>> SPextensions;
    
    // Class that will build a starting YAML config for all extensions
    //YAMLCollector inputYAML;// = YAMLCollector("", PixelMapCollection::magicKey);

    PixelMapCollection mapCollection;

    //PixelMapCollection* pmcInit;
    
    // List of all Matches - they will hold pointers to all Detections too.
    MCat matches;

    vector<SphericalCoords*> extensionProjections;//extensions.size(), nullptr);  
    
    //set<string> degenerateTypes; //={"Poly","Linear","Constant"};
    //void addInputYAML(string inputMaps);

    void setExposures(vector<shared_ptr<Exposure>> expos, double sysErr, double refSysErr);
    void setExposures(vector<unique_ptr<Exposure>> expos, double sysErr, double refSysErr);
    
    void setExtensions(vector<shared_ptr<Extension>> extens);

    void setRefWCSNames();

    void addMap(YAMLCollector& inputYAML, string mapName, vector<string> mapParams);
    //void addMap(string mapName, vector<string> mapParams);

    void setupMaps(YAMLCollector& inputYAML);//, PixelMapCollection& mapCollection);
    
    void setMatches(vector<int> sequence, vector<LONGLONG> extensions, vector<LONGLONG> objects,
                    ExtensionObjectSet skipSet);

    void setObjects(int i, img::FTable ff, string xKey, string yKey, string idKey, string pmCovKey,
                    vector<string> xyErrKeys, string magKey, int magKeyElement, string magErrKey,
                    int magErrKeyElement, string pmRaKey, string pmDecKey, string parallaxKey);

    void defaultMaps();

    void reprojectWCSs();

    void fit();

    int getMatchLength();

    void cleanup();
      
};

#endif