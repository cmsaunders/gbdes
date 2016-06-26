// Subroutines for use by WCSFit.cpp
#ifndef WCSSUBS_H
#define WCSSUBS_H

#include <list>
#include <set>
#include "StringStuff.h"
#include "NameIndex.h"
#include "Astrometry.h"
#include "FitsTable.h"

// Load the classes needed
#include "FitSubroutines.h"
#include "PixelMapCollection.h"
#include "Match.h"
#include "Instrument.h"
typedef ExtensionBase<astrometry::SubMap, astrometry::Detection> Extension;
typedef ColorExtensionBase<astrometry::Match> ColorExtension;

// Function that will using starting WCS to fit all of the defaulted
// maps used by the selected extensions.  Then will put the
// initialized parameters back into the PMC and clear the defaulted flag.
void fitDefaulted(astrometry::PixelMapCollection& pmc,
		  set<Extension*> useThese,
		  const vector<Instrument*>& instruments,
		  const vector<Exposure*>& exposures);


// Statistics-accumulator class. 
class Accum {
public:
  double sumx;
  double sumy;
  double sumxx;
  double sumyy;
  double sumxxw;
  double sumyyw;
  double sumxw;
  double sumyw;
  double sumwx;
  double sumwy;
  int n;
  double xpix;
  double ypix;
  double xw;
  double yw;
  double sumdof;
  Accum();
  void add(const astrometry::Detection* d, double xoff=0., double yoff=0., double dof=1.);
  double rms() const;
  double reducedChisq() const;
  string summary() const;
};

// Read the Fields table from input, copy to output, extract needed info
void
readFields(string inputTables,
	   string outCatalog,
	   NameIndex& fieldNames,
	   vector<astrometry::SphericalCoords*>& fieldProjections);

#endif