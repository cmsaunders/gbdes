// Class used to figure out how to break degeneracies between sequential linear/poly
// maps.

#include "FitSubroutines.h"
#include <set>
#include <list>
#include <string>
#include <vector>


template <class S>
class MapDegeneracies {
public:
    // Find potentially degenerate map elements, build graphs
    // of which extensions are using which maps.
    // Bool determines whether we take any map with free parameters vs
    // examining only those with defaulted parameters.
    // The list specifies which types of maps are potentially degenerate
    // with each other. If the list is empty, all maps are considered.
    MapDegeneracies(const std::vector<std::unique_ptr<typename S::Extension>> &extensions_,
                    const typename S::Collection &mapCollection, const std::set<string> &mapTypes_,
                    bool defaulted = false);
    ~MapDegeneracies() {}

    // Determine which of the possible candidate maps
    // must be converted to Identity maps (or frozen) in order to
    // break degeneracies.  Will exit with error if
    // no path is available.
    list<string> replaceWithIdentity(const std::set<string> &candidates);

    // Return the order in which extensions should be initialized
    // to do so without ever having degenerate maps.  Do all in each
    // set at the same time, e.g. a full exposure.
    std::list<std::set<int>> initializationOrder();

private:
    // Find all extensions using exactly one map
    std::set<int> findNondegenerate() const;
    // Remove all edges for a map and any extns that then have no edges
    void eraseMap(string mapname);

    // All the maps, each with all the extensions that use it
    std::map<string, std::set<int>> maps;
    // All the extensions, each with all the maps that use them
    std::map<int, std::set<string>> extns;
    // The extension table (extn numbers are index into this vector)
    const std::vector<std::unique_ptr<typename S::Extension>> &extensions;
};
