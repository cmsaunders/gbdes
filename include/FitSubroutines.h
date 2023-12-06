// Constants and subroutines shared by the astrometric and photometric matching/fitting codes
// WCSFit and PhotoFit

#ifndef FITSUBROUTINES_H
#define FITSUBROUTINES_H

// Number of threads to use for reading catalogs
#define CATALOG_THREADS 4

#include <list>
#include <memory>
#include <set>
#include "StringStuff.h"

// Load the kinds of maps we'll want
#include "PixelMapCollection.h"
#include "PhotoMapCollection.h"
#include "PhotoTemplate.h"
#include "TemplateMap.h"

#include "Pset.h"
#include "FTable.h"
#include "FitsTable.h"
#include "Instrument.h"
#include "YAMLCollector.h"

#include "Match.h"
#include "PhotoMatch.h"
#include "Accum.h"


const int REF_INSTRUMENT = -1;  // Instrument for reference objects (no fitting)
const int PM_INSTRUMENT = -2;   // Reference catalog with 5d proper motion / parallax solutions
const int TAG_INSTRUMENT = -3;  // Exposure number for tag objects (no fitting nor contrib to stats)
const int NO_INSTRUMENT = -4;

const string stellarAffinity = "STELLAR";

const double NO_MAG_DATA = -100.;  // Value entered when there is no valid mag or color

const string colorColumnName = "COLOR";
const string colorErrorColumnName = "COLOR_ERR";

// Here is the default character at which to split lists given in parameters strings
const char DefaultListSeperator = ',';

// These structures are used so that we can write templated code that will work
// for either astrometric or photometric fitting.
struct Astro {
    typedef astrometry::Detection Detection;
    typedef astrometry::Match Match;
    typedef astrometry::SubMap SubMap;
    typedef ExtensionBase<SubMap, Detection> Extension;
    typedef ColorExtensionBase<Match> ColorExtension;
    typedef astrometry::PixelMapCollection Collection;
    typedef astrometry::CoordAlign Align;
    typedef astrometry::MCat MCat;
    static void fillDetection(Detection &d, const Exposure &e, astrometry::SphericalCoords &fieldProjection,
                              const img::FTable &table, long irow, string xKey, string yKey,
                              const std::vector<std::string> &xyErrKeys, std::string magKey, std::string magErrKey,
                              int magKeyElement, int magErrKeyElement, bool xColumnIsDouble,
                              bool yColumnIsDouble, bool errorColumnIsDouble, bool magColumnIsDouble,
                              bool magErrColumnIsDouble, double magshift,
                              const astrometry::PixelMap *startWcs, bool isTag);
    static void setColor(Detection &d, double color) { d.color = color; }
    static double getColor(const Detection &d) { return d.color; }

    // saving results wants an array giving the field projection to use for each catalog.
    static void saveResults(const astrometry::MCat &matches, std::string outCatalog, std::string starCatalog,
                            std::vector<astrometry::SphericalCoords *> catalogProjections);
    // writes out the final residuals catalog as in `saveResults`, but returns results as a map
    struct outputCatalog {
        std::vector<int> matchID;
        std::vector<long> catalogNumber;
        std::vector<long> objectNumber;
        std::vector<bool> clip;
        std::vector<bool> reserve;
        std::vector<bool> hasPM;
        std::vector<double> color;
        std::vector<double> xresw;
        std::vector<double> yresw;
        std::vector<double> xpix;
        std::vector<double> ypix;
        std::vector<double> sigpix;
        std::vector<double> xrespix;
        std::vector<double> yrespix;
        std::vector<double> xworld;
        std::vector<double> yworld;
        std::vector<double> covTotalW_00;
        std::vector<double> covTotalW_11;
        std::vector<double> covTotalW_01;
        std::vector<double> chisq;
        std::vector<double> chisqExpected;
    };
    static outputCatalog getOutputCatalog(const astrometry::MCat &matches);

    struct PMCatalog {
        std::vector<int> pmMatchID;
        std::vector<long> pmCatalogNumber;
        std::vector<long> pmObjectNumber;
        std::vector<bool> pmClip;
        std::vector<bool> pmReserve;
        std::vector<std::vector<double>> pmMean;
        std::vector<std::vector<double>> pmInvCov;
        std::vector<double> pmChisq;
        std::vector<double> pmChisqExpected;
    };
    // writes out the residuals for PMDetections as in `saveResults`, but returns results as a map
    static PMCatalog getPMCatalog(const astrometry::MCat &matches);

    struct StarCatalog {
        std::vector<int> starMatchID;
        std::vector<bool> starReserve;
        std::vector<double> starColor;
        std::vector<int> starPMCount;
        std::vector<int> starDetCount;
        std::vector<int> starClipCount;
        std::vector<int> starDOF;
        std::vector<double> starChisq;
        std::vector<double> starX;
        std::vector<double> starY;
        std::vector<double> starPMx;
        std::vector<double> starPMy;
        std::vector<double> starParallax;
        std::vector<std::vector<double>> starInvCov;
    };
    // writes out the catalog of fit star positions (and proper motions if fit) as in `saveResults`, but
    // returns results as a map
    static StarCatalog getStarCatalog(const astrometry::MCat &matches,
                                      std::vector<astrometry::SphericalCoords *> catalogProjections);

    static void reportStatistics(const astrometry::MCat &matches,
                                 const std::vector<std::unique_ptr<Exposure>> &exposures,
                                 const std::vector<std::unique_ptr<Astro::Extension>> &extensions, std::ostream &os);

    // Specialized function for reading detection from a catalog
    // holding proper motion solution.
    // It will incorporate the basic info in (non-PM) Detection d
    // into a new PMDetection object and return it.
    static std::unique_ptr<astrometry::PMDetection> makePMDetection(
            const astrometry::Detection &d, const Exposure &e, const img::FTable &table, long irow,
            std::string xKey, std::string yKey, std::string pmRaKey, std::string pmDecKey, std::string parallaxKey, std::string pmCovKey,
            bool xColumnIsDouble, bool yColumnIsDouble, bool errorColumnIsDouble,
            const astrometry::PixelMap *startWcs);

    static void handlePMDetection(std::unique_ptr<astrometry::PMDetection> pmd, Detection const &d);

    static std::unique_ptr<Match> makeNewMatch(std::unique_ptr<Detection> d, bool usePM);

    static const int isAstro = 1;
};
struct Photo {
    typedef photometry::Detection Detection;
    typedef photometry::Match Match;
    typedef list<std::unique_ptr<Match>> MCat;
    typedef photometry::SubMap SubMap;
    typedef ExtensionBase<SubMap, Detection> Extension;
    typedef ColorExtensionBase<Match> ColorExtension;
    typedef photometry::PhotoMapCollection Collection;
    typedef photometry::PhotoAlign Align;
    static void fillDetection(Detection &d, const Exposure &e, astrometry::SphericalCoords &fieldProjection,
                              const img::FTable &table, long irow, std::string xKey, std::string yKey,
                              const std::vector<std::string> &xyErrKeys, std::string magKey, std::string magErrKey,
                              int magKeyElement, int magErrKeyElement, bool xColumnIsDouble,
                              bool yColumnIsDouble, bool errorColumnIsDouble, bool magColumnIsDouble,
                              bool magErrColumnIsDouble, double magshift,
                              const astrometry::PixelMap *startWcs, bool isTag);
    static void setColor(Detection &d, double color) { d.args.color = color; }
    static double getColor(const Detection &d) { return d.args.color; }
    static void saveResults(const MCat &matches, std::string outCatalog);

    static void reportStatistics(const MCat &matches, const std::vector<std::unique_ptr<Exposure>> &exposures,
                                 const std::vector<std::unique_ptr<Photo::Extension>> &extensions, std::ostream &os);

    // This is a no-op for photometry:
    static std::unique_ptr<astrometry::PMDetection> makePMDetection(
            photometry::Detection const d, const Exposure &e, const img::FTable &table, long irow,
            std::string xKey, std::string yKey, std::string pmRaKey, std::string pmDecKey, std::string parallaxKey, std::string pmCovKey,
            bool xColumnIsDouble, bool yColumnIsDouble, bool errorColumnIsDouble,
            const astrometry::PixelMap *startWcs) {
        return nullptr;
    }

    static void handlePMDetection(std::unique_ptr<astrometry::PMDetection> pmd,
                                  Detection const &d){};  // This is a no-op for photo

    static std::unique_ptr<Match> makeNewMatch(std::unique_ptr<Detection> d, bool usePM) {
        return std::unique_ptr<Match>(new Match(std::move(d)));
    }
    static const int isAstro = 0;
};

// A helper function that strips white space from front/back of a string and replaces
// internal white space with underscores:
void spaceReplace(std::string &s);

// Another helper function to split up a string into a list of whitespace-trimmed strings.
// Get rid of any null ones.
list<std::string> splitArgument(std::string input, const char listSeperator = DefaultListSeperator);

// Read parameters from files on command line and from std input.
// First nRequiredArgs must be present and are not names of parameter files.
// Print the usage string (to cerr) and the default parameters (to cout) if
// there are not enough parameters.
// Then dumps final parameter values to cout.
void processParameters(Pset &parameters, const std::string &usage, int nRequiredArgs, int argc, char *argv[]);

// Take a string with this format:
// <regex>=<replace> [, <regex>=<replace>] ...
// and parse it into a stringstuff::RegexReplacement object that will execute
// the regular expression replacements on an input string.
stringstuff::RegexReplacements parseTranslator(std::string specString, std::string errorDescription);

// This class reads a set of exposure/object number pairs from a file on construction.
// Then it can be used as a function that returns true if a given pair was on the list.
// First it needs the class holding the pairs (user does not need this)

class EOPair {
public:
    EOPair(int exp, long obj) : extensionNumber(exp), objectNumber(obj) {}
    bool operator<(const EOPair &rhs) const {
        return extensionNumber < rhs.extensionNumber ||
               (extensionNumber == rhs.extensionNumber && objectNumber < rhs.objectNumber);
    }

private:
    int extensionNumber;
    long objectNumber;
};

class ExtensionObjectSet {
public:
    ExtensionObjectSet(std::string filename);
    bool operator()(int extensionNumber, long objectNumber) const;
    void clear() { pairs.clear(); }

private:
    std::set<EOPair> pairs;
};

// These functions fill static tables with all the types of Photo/PixelMaps that we
// may be deserializing:
void loadPixelMapParser();
void loadPhotoMapParser();

// Routines to handle table entries that may be float or double, and might be array elements
// First one: parse a key to see if it has form COLUMN[#].  The integer # is returned and the
// key is return as COLUMN.  If no [#], key is unchanged and elementNumber=-1 is returned.
int elementNumber(std::string &key);
// Next one: see if the table column is double (true) else assume float (false)
bool isDouble(img::FTable f, std::string key, int elementNumber);
// Retrieve a double-valued number from either float or double column, element of array or scalar cell
double getTableDouble(img::FTable f, std::string key, int elementNumber, bool isDouble, long irow);

// Class to produce ordering for vector of points to objects
// that have a public "name" member.
template <class T>
class NameSorter {
public:
    NameSorter(const std::vector<std::unique_ptr<T>> &exposures) : ve(exposures) {}
    bool operator()(int i1, int i2) const {
        LessNoCase lnc;
        return lnc(ve[i1]->name, ve[i2]->name);
    }

private:
    const std::vector<std::unique_ptr<T>> &ve;
};

// This function is used to find degeneracies between exposures and device maps.
// Start with list of free & fixed devices as initial degen/ok, same for exposures.
// Will consider as "ok" any device used in an "ok" exposure and vice-versa.
// The last argument says which exposure/device pairs are used together.
void findDegeneracies(std::set<int> &degenerateDevices, std::set<int> &okDevices, std::set<int> &degenerateExposures,
                      std::set<int> &okExposures, const std::vector<std::set<int>> &exposuresUsingDevice);

// Figure out which extensions of the FITS file inputTables
// are Instrument or MatchCatalog extensions.
void inventoryFitsTables(std::string inputTables, std::vector<int> &instrumentHDUs, std::vector<int> &catalogHDUs);

class Fields {
public:
    // Construct directly from vectors of values.
    //
    // This does _not_ take care of setting default epochs, like `read` does.
    //
    // Vectors are passed by value because we generally expect them to be
    // constructed on-the-fly from Python lists or moved from locals in `read`,
    // so by being careful about using moves we can avoid any copies (not that
    // this really matters).
    //
    // Should be exposed to Python.
    Fields(std::vector<std::string> names, std::vector<double> ra, std::vector<double> dec, std::vector<double> epochs);

    // Default constructor.
    //
    // This is currently used only by FitClass, and should be removed if that
    // class is; it's better for this class to guarantee that it has been
    // initialized with meaningful values whenever an instance exists.
    Fields() = default;

    // Read the Fields table from input, copy to output, extract needed info
    // If a usable default is given, unassigned field epochs are set to default.
    //
    // Need not be exposed to Python but may be if useful for debugging.
    static Fields read(std::string inputTables, std::string outCatalog, double defaultEpoch = 0.);

    // Access to names.
    //
    // Should not be exposed to Python.
    NameIndex const &names() const { return _names; }

    // Access to projections.
    //
    // Should not be exposed to Python.
    std::vector<std::unique_ptr<astrometry::SphericalCoords>> const &projections() const { return _projections; }

    // Access to epochs.
    //
    // Should not be exposed to Python.
    std::vector<double> const &epochs() const { return _epochs; }

private:
    NameIndex _names;
    std::vector<std::unique_ptr<astrometry::SphericalCoords>> _projections;
    std::vector<double> _epochs;
};

struct FieldsHelper {
public:
    FieldsHelper(std::vector<std::string> names_, std::vector<double> ra_, std::vector<double> dec_, std::vector<double> epochs_)
            : names(names_), ra(ra_), dec(dec_), epochs(epochs_){};

    std::vector<std::string> names;
    std::vector<double> ra;
    std::vector<double> dec;
    std::vector<double> epochs;
};

// Read in all the instrument extensions and their device info from input
// FITS file, save useful ones and write to output FITS file.
// The useInstrumentList entries are regexes, empty means use all.
// The final bool argument is set true if we have already created
// the outCatalog FITS file.
std::vector<std::unique_ptr<Instrument>> readInstruments(const std::vector<int> &instrumentHDUs,
                                               const list<std::string> &useInstrumentList, std::string inputTables,
                                               std::string outCatalog, bool &outputCatalogAlreadyOpen);

// Read the Exposure table into an array.
// Uses the instruments table.
// Fills the exposureColorPriorities table using the useColorPriorities
// list of regexes to establish priority order.
// useReference exposures is set if we want to use exposures from REFERENCE instrument.
// last Boolean is as above.
std::vector<std::unique_ptr<Exposure>> readExposures(const std::vector<std::unique_ptr<Instrument>> &instruments,
                                           const std::vector<double> &fieldEpochs,
                                           std::vector<int> &exposureColorPriorities,
                                           const list<std::string> &useColorList, std::string inputTables,
                                           std::string outCatalog, const list<std::string> &skipExposureList,
                                           bool useReferenceExposures, bool &outputCatalogAlreadyOpen);

std::vector<std::shared_ptr<astrometry::Wcs>> readWCSs(const img::FTable &extensionTable);

// Read extensions from the table.
// colorExtensions will get filled with ColorExtension objects for color data
// inputYAML is set to produce YAML for all extensions being fit.
// Systematic error for Detections in this Extension will be assigned from
// the column in the extensionTable named sysErrorColumn, if it is non-Null,
// otherwise taken from the value of sysError or referenceSysError as appropriate
template <class S>
std::vector<std::unique_ptr<typename S::Extension>> readExtensions(
        const img::FTable &extensionTable, const std::vector<std::unique_ptr<Instrument>> &instruments,
        const std::vector<std::unique_ptr<Exposure>> &exposures, const std::vector<int> &exposureColorPriorities,
        std::vector<std::unique_ptr<typename S::ColorExtension>> &colorExtensions, astrometry::YAMLCollector &inputYAML,
        bool logging = true);

// fix all maps in a photo/pixelMapCollection whose names match
// any regex in the fixMapList.  Also any instrument whose name
// matches gets all device maps fixed.
template <class S>
void fixMapComponents(typename S::Collection &pmc, const list<std::string> &fixMapList,
                      const std::vector<std::unique_ptr<Instrument>> &instruments);

template <class S>
int findCanonical(Instrument &instr, int iInst, std::vector<std::unique_ptr<Exposure>> &exposures,
                  std::vector<std::unique_ptr<typename S::Extension>> &extensions, typename S::Collection &pmc);

// Add every extension's map to the YAMLCollector and then emit the
// YAML and read into the MapCollection.
// The names of all maps are already in the extension list.
// Returns time spent in critical regions of the addMap() routine.
template <class S>
void createMapCollection(const std::vector<std::unique_ptr<Instrument>> &instruments,
                         const std::vector<std::unique_ptr<Exposure>> &exposures,
                         const std::vector<std::unique_ptr<typename S::Extension>> &extensions,
                         astrometry::YAMLCollector &inputYAML, typename S::Collection &pmc);

// Set the bool in each Extension indicating whether its map uses color
template <class S>
void whoNeedsColor(const std::vector<std::unique_ptr<typename S::Extension>> &extensions);

// Read a MatchCatalog Extension, recording in each extension the
// objects from it that need to be read from catalog.
// skipSet can give objects to ignore;
// Matches with less than minMatches useful inputs are discarded too.
// Discards Detection if requires color but doesn't have it.
template <class S>
void readMatches(const std::vector<int> &seq, const std::vector<LONGLONG> &extn, const std::vector<LONGLONG> &obj,
                 typename S::MCat &matches, const std::vector<std::unique_ptr<typename S::Extension>> &extensions,
                 const std::vector<std::unique_ptr<typename S::ColorExtension>> &colorExtensions,
                 const ExtensionObjectSet &skipSet, int minMatches,
                 bool usePM = false);  // If true, create PMMatches

template <class S>
void readMatches(const img::FTable &table, typename S::MCat &matches,
                 const std::vector<std::unique_ptr<typename S::Extension>> &extensions,
                 const std::vector<std::unique_ptr<typename S::ColorExtension>> &colorExtensions,
                 const ExtensionObjectSet &skipSet, int minMatches,
                 bool usePM = false);  // If true, create PMMatches

// Read each Extension's objects' data from it FITS catalog
// and place into Detection structures.
template <class S>
void readObjects(const img::FTable &extensionTable, const std::vector<std::unique_ptr<Exposure>> &exposures,
                 const std::vector<std::unique_ptr<typename S::Extension>> &extensions,
                 const std::vector<std::unique_ptr<astrometry::SphericalCoords>> &fieldProjections,
                 bool logging = true);  // Give progress updates?

template <class S>
void readObjects_oneExtension(const std::vector<std::unique_ptr<Exposure>> &exposures, int iext, const img::FTable &ff,
                              const std::string &xKey, const std::string &yKey, const std::string &idkey,
                              const std::string &pmCovKey, const std::vector<std::string> &xyErrKeys, const std::string &magKey,
                              const int &magKeyElement, const std::string &magErrKey, const int &magErrKeyElement,
                              const std::string &pmRaKey, const std::string &pmDecKey, const std::string &parallaxKey,
                              const std::vector<std::unique_ptr<typename S::Extension>> &extensions,
                              const std::vector<std::unique_ptr<astrometry::SphericalCoords>> &fieldProjections,
                              bool logging, bool useRows);

// Read color information from files marked as holding such, insert into
// relevant Matches.
template <class S>
void readColors(const img::FTable &extensionTable,
                const std::vector<std::unique_ptr<typename S::ColorExtension>> &colorExtensions,
                bool logging = true);  // Progress reports

// Find all matched Detections that exceed allowable error, then
// delete them from their Match and delete the Detection.
// Does not apply to reference/tag detections.
template <class S>
void purgeNoisyDetections(double maxError, typename S::MCat &matches,
                          const std::vector<std::unique_ptr<Exposure>> &exposures,
                          const std::vector<std::unique_ptr<typename S::Extension>> &extensions);

// Get rid of Matches with too few Detections being fit: delete
// the Match and all of its Detections.
template <class S>
void purgeSparseMatches(int minMatches, typename S::MCat &matches);

// Get rid of Matches with color outside of specified range.
// Color is always ok if it has NODATA value.
// Note that this even kills Detections that do not need a color for their maps.
template <class S>
void purgeBadColor(double minColor, double maxColor, typename S::MCat &matches);

template <class S>
void reserveMatches(typename S::MCat &matches, double reserveFraction, int randomNumberSeed);

template <class S>
std::map<std::string, long> findUnderpopulatedExposures(long minFitExposure, const typename S::MCat &matches,
                                              const std::vector<std::unique_ptr<Exposure>> &exposures,
                                              const std::vector<std::unique_ptr<typename S::Extension>> &extensions,
                                              const typename S::Collection &pmc);

// Fix the parameters of a map, and mark all Detections making
// use of it as clipped so they will not be used in fitting
template <class S>
void freezeMap(std::string mapName, typename S::MCat &matches,
               std::vector<std::unique_ptr<typename S::Extension>> &extensions, typename S::Collection &pmc);

// Report number of unclipped matches and their chisq
template <class S>
void matchCensus(const typename S::MCat &matches, std::ostream &os);

// Map and clip reserved matches
template <class S>
void clipReserved(typename S::Align &ca, double clipThresh, double minimumImprovement, bool clipEntireMatch,
                  bool reportToCerr);

// Function to produce a list of PhotoPriors from an input file
// (for PhotoFit only - has its own source file)
list<photometry::PhotoPrior *> readPriors(std::string filename, const std::vector<std::unique_ptr<Instrument>> &instruments,
                                          const std::vector<std::unique_ptr<Exposure>> &exposures,
                                          const std::vector<std::unique_ptr<Photo::Extension>> &extensions);

#endif
