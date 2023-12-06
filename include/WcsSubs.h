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
typedef Astro::Extension Extension;
typedef Astro::ColorExtension ColorExtension;

// Function that will using starting WCS to fit all of the defaulted
// maps used by the selected extensions.  Then will put the
// initialized parameters back into the PMC and clear the defaulted flag.
void fitDefaulted(astrometry::PixelMapCollection &pmc, std::set<Extension *> useThese,
                  const std::vector<std::unique_ptr<Instrument>> &instruments,
                  const std::vector<std::unique_ptr<Exposure>> &exposures, bool logging = true);

// Define and issue WCS for each extension in use, and set projection to
// field coordinates.
void setupWCS(const std::vector<std::unique_ptr<astrometry::SphericalCoords>> &fieldProjections,
              const std::vector<std::unique_ptr<Instrument>> &instruments,
              const std::vector<std::unique_ptr<Exposure>> &exposures, std::vector<std::unique_ptr<Extension>> &extensions,
              astrometry::PixelMapCollection &pmc);

// Analyze the PixelMap to find list of exposures that we can
// initialize first to set up all defaulted device maps.
list<int> pickExposuresToInitialize(const std::vector<std::unique_ptr<Instrument>> &instruments,
                                    const std::vector<std::unique_ptr<Exposure>> &exposures,
                                    const std::vector<std::unique_ptr<Extension>> &extensions,
                                    astrometry::PixelMapCollection &pmc);

#endif
