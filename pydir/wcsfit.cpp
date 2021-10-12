#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/eigen.h"
#include <pybind11/stl_bind.h>

//# include "ndarray/pybind11.h"
//#include <tuple>
//#include <string.h>

//#include "WCSFit.h"
//#include "../src/WCSFit.cpp"

//#include "../src/model.hpp"
//#include "../src/test.cpp"
//#include "../src/subs/FitSubroutines.cpp"
//#include "../../gb_packages/astrometry/src/PiecewiseMap.cpp"
//#include "../../yaml-cpp/src/exceptions.cpp"
#include "FitSubroutines.h"
#include "TPVMap.h"
#include "WCSFoF_match.h"
//#include "Units.h"
/*
#include "Std.h"
#include "Astrometry.h"
#include "FitsTable.h"
#include "StringStuff.h"
#include "Pset.h"
#include "PixelMapCollection.h"
#include "YAMLCollector.h"
#include "Match.h"
#include "Instrument.h"
#include "Match.h"

#include "WcsSubs.h"
#include "StringStuff.h"
//#include "MapDegeneracies.h"
*/
namespace py = pybind11;
//using namespace astrometry;

//using astrometry::PMMatch;


//double multiplyByTwo(double num){
//	Operations op;
//	return op.numeric_product(num, 2.0);
//}

template<class T1, class T2> 
void declareExtension(py::module &m) {
    using Class = ExtensionBase<T1, T2>;

    py::class_<Class>(m, "Extension")
        .def(py::init<>())
        .def_readwrite("exposure", &Class::exposure)
        .def_readwrite("device", &Class::device)
        .def_readwrite("map", &Class::map)
        .def_readwrite("airmass", &Class::airmass)
        .def_readwrite("apcorr", &Class::apcorr)
        .def_readwrite("magshift", &Class::magshift)
        .def_readwrite("wcs", &Class::wcs)
        .def_readwrite("startWcs", &Class::startWcs)
        .def_readwrite("keepers", &Class::keepers);
}

template<class T>
void declareColorExtension(py::module &m) {
    using Class = ColorExtensionBase<T>;

    py::class_<Class>(m, "ColorExtension")
        .def(py::init<>());
}

template<class T = float>
void declareBounds(py::module &m) {
    using Class = Bounds<T>;

    py::class_<Class>(m, "Bounds")
        .def(py::init<const T, const T, const T, const T>());

}

template<class T>
void declareVector(py::module &m) {
    using Class = std::vector<T>;

    py::class_<Class>(m, "deviceVector");
}

template<class T>
void declareList(py::module &m, string ListName) {
    using Class = std::list<T>;

    py::class_<Class>(m, ListName.c_str());
}

//PYBIND11_MAKE_OPAQUE(std::list<Point>);

template<class P, int DIM = 2>
void declareMatch(py::module &m) {
    using Class = fof::Match<P, DIM>;
    //using Class2 = std::list<P const*>;
    //declareList<P>(m);

    py::class_<Class>(m, "FoFMatch")
        //.def(py::init<>())
        //.def_static("fromVector", [](Class & self, std::vector<P> v){
        //    for 
        //    self{ std::begin(v), std::end(v)};
        //})
        .def("toVector", [](Class & self){
            //std::vector<const P*> v{ std::begin(self), std::end(self) };
            std::vector<const P*> v{ std::begin(self), std::end(self) };
            return v;
        });
    //py::bind_vector<Class, Class2>(m, "FoFMatch");
}

template<class P, int DIM>
void declareMatchSet(py::module &m) {
    using Class = set<fof::Match<P, DIM>*>;

    py::class_<Class>(m, "MatchSet");
}

template<class P, int DIM>
void declareCatalog(py::module &m, string CatalogType) {
    using Class = fof::Catalog<P, DIM>;

    py::class_<Class>(m, CatalogType.c_str())
        .def("toVector", [](Class & self){
            std::vector<fof::Match<P, DIM>*> v{ std::begin(self), std::end(self) };
            return v;
        });
    //py::class_<Class>(m, CatalogType.c_str());
}

template<class S>
void declareReadMatches(py::module &m) {
    m.def("readMatches", readMatches<S>);
}

PYBIND11_MODULE(wcsfit, m) {
    
    m.doc() = "tmp docstring"; // optional module docstring

    //m.def("multiplyByTwo", &multiplyByTwo,
    //      "obtain the double of a number");
    //m.def("add", &testfn, "Run add");
    //m.def("run", &run, "test main");
    //m.def("testfn", &testfn);
    m.attr("REF_INSTRUMENT") = REF_INSTRUMENT;
    m.attr("PM_INSTRUMENT") = PM_INSTRUMENT;
    m.attr("ARCSEC") = ARCSEC;
    m.attr("MILLIARCSEC") = MILLIARCSEC;
    m.attr("DEGREE") = DEGREE;
    //m.attr("RESIDUAL_UNIT") = RESIDUAL_UNIT;
    //m.attr("WCS_UNIT") = WCS_UNIT;
    m.def("spaceReplace", &spaceReplace);
    m.def("loadPixelMapParser", &loadPixelMapParser);
    //declareReadMatches<Astro>(m);
    m.def("readWrite", py::overload_cast<vector<int>&, vector<LONGLONG>&, vector<LONGLONG>&,
                                         Astro::MCat&, vector<Astro::Extension*>&,
                                         vector<Astro::ColorExtension*>&, ExtensionObjectSet const&,
                                         int, bool>
                                         (&readMatches<Astro>));

    py::class_<astrometry::SphericalCoords>(m, "SphericalCoords");

    // TODO: check which of these are actually needed
    py::class_<Instrument>(m, "Instrument")
        .def(py::init<string>(), py::arg("name_")="")
        .def_readwrite("name", &Instrument::name)
        .def_readwrite("nDevices", &Instrument::nDevices)
        .def_readwrite("band", &Instrument::band)
        .def_readwrite("deviceNames", &Instrument::deviceNames)
        .def("addDevice", &Instrument::addDevice);
    
    py::class_<Exposure>(m, "Exposure")
        .def(py::init<string, astrometry::SphericalCoords const&>())
        .def_readwrite("name", &Exposure::name)
        .def_readwrite("projection", &Exposure::projection)
        .def_readwrite("field", &Exposure::field)
        .def_readwrite("instrument", &Exposure::instrument)
        .def_readwrite("airmass", &Exposure::airmass)
        .def_readwrite("exptime", &Exposure::exptime)
        .def_readwrite("mjd", &Exposure::mjd)
        .def_readwrite("apcorr", &Exposure::apcorr)
        .def_readwrite("epoch", &Exposure::epoch)
        .def_readwrite("projection", &Exposure::projection)
        .def("setAstrometricCovariance", [](Exposure & self, astrometry::Matrix22::Base const& matrix){self.astrometricCovariance = matrix;})
        .def("getAstrometricCovariance", [](Exposure const& self){return static_cast <astrometry::Matrix22::Base const&>(self.astrometricCovariance);})
        .def("addToAstrometricCovariance", [](Exposure & self, astrometry::Matrix22::Base const& matrix){self.astrometricCovariance += matrix;});

    py::class_<astrometry::Detection>(m, "Dectection");
    
    py::class_<astrometry::PixelMap>(m, "PixelMap");

    py::class_<astrometry::SubMap, astrometry::PixelMap>(m, "SubMap")
        .def(py::init<list<astrometry::PixelMap*> const&, string, bool>());

    declareExtension<astrometry::SubMap, astrometry::Detection>(m);

    declareColorExtension<astrometry::Match>(m);
        
    py::class_<astrometry::Gnomonic, astrometry::SphericalCoords>(m, "Gnomonic")
        .def(py::init<astrometry::Orientation const&, bool>(), py::arg("o"), py::arg("shareOrient")=false)
        .def(py::init<astrometry::SphericalCoords const&, astrometry::Orientation&>())
        .def(py::init<>());
        

    py::class_<astrometry::SphericalICRS, astrometry::SphericalCoords>(m, "SphericalICRS")
        .def(py::init<double, double>());
    
    py::class_<astrometry::Orientation>(m, "Orientation")
        .def(py::init<>())
        .def(py::init<astrometry::SphericalCoords&, double>(), py::arg("pole_"), py::arg("pa")=0)
        .def("set", &astrometry::Orientation::set)
        .def("getPole", &astrometry::Orientation::getPole);

    py::class_<astrometry::LinearMap, astrometry::PixelMap>(m, "LinearMap")
        //.def("setVector", [](astrometry::LinearMap & self, astrometry::DVector::Base const& vector){self.v = vector;});
        .def(py::init<astrometry::DVector::Base const&>());

    py::class_<poly2d::Poly2d>(m, "Poly2d")
        .def(py::init<int>())
        .def(py::init<astrometry::DMatrix::Base const&>())
        .def("nCoeffs", &poly2d::Poly2d::nCoeffs)
        .def("getC", &poly2d::Poly2d::getC)
        .def("setC", &poly2d::Poly2d::setC)
        .def("vectorIndex", &poly2d::Poly2d::vectorIndex);

    py::class_<astrometry::PolyMap, astrometry::PixelMap>(m, "PolyMap")
        .def(py::init<poly2d::Poly2d, poly2d::Poly2d, string, Bounds<double>, double>());

    declareBounds<double>(m);

    py::class_<astrometry::Wcs, astrometry::PixelMap>(m, "Wcs")
        .def(py::init<astrometry::PixelMap*, astrometry::SphericalCoords const&, string, double, bool>(),
            py::arg("pm_"), py::arg("nativeCoords_"), py::arg("name")="", py::arg("wScale_")=DEGREE,
            py::arg("shareMap_")=false)
        .def("reprojectTo", &astrometry::Wcs::reprojectTo);

    py::class_<astrometry::Match>(m, "Match");

    py::class_<ExtensionObjectSet>(m, "ExtensionObjectSet")
        .def(py::init<string>());


    // The following are WCSFoF-specific classes:
    py::class_<Point>(m, "Point")
        .def(py::init<double, double, long int, long int, long int>())
        .def_readwrite("extensionNumber", &Point::extensionNumber)
        .def_readwrite("objectNumber", &Point::objectNumber)
        .def_readwrite("exposureNumber", &Point::exposureNumber);

    py::class_<Field>(m, "Field")
        .def(py::init<>())
        .def_readwrite("name", &Field::name)
        .def_readwrite("matchRadius", &Field::matchRadius)
        .def_readwrite("projection", &Field::projection)
        .def_readwrite("extent", &Field::extent)
        .def_readwrite("catalogs", &Field::catalogs);

    declareList<Point const*>(m, "PointList");

    declareMatch<Point, 2>(m);

    declareMatchSet<Point, 2>(m);

    declareCatalog<Point, 2>(m, "PointCat");

    py::bind_map<map<string, PointCat>>(m, "PointCatDict");

    py::class_<Device, Bounds<double>>(m, "Device")
        .def(py::init<>())
        //.def(py::init<double const, double const, double const, double const>())
        .def("setXMin", &Device::setXMin)
        .def("setXMax", &Device::setXMax)
        .def("setYMin", &Device::setYMin)
        .def("setYMax", &Device::setYMax)
        .def_readwrite("name", &Device::name);

    declareVector<Device>(m);

    py::bind_vector<Instr, vector<Device>>(m, "Instr")
        .def(py::init<string>())
        .def_readwrite("name", &Instr::name);
        
    py::class_<Expo>(m, "Expo")
        .def(py::init<>())
        .def_readwrite("name", &Expo::name)
        .def_readwrite("field", &Expo::field)
        .def_readwrite("instrument", &Expo::instrument)
        .def_readwrite("pointing", &Expo::pointing);

    py::class_<FoFClass>(m, "FoFClass")
        .def(py::init<>())
        .def_readwrite("useAffinities", &FoFClass::useAffinities)
        .def_readwrite("minMatches", &FoFClass::minMatches)
        .def_readwrite("allowSelfMatches", &FoFClass::allowSelfMatches)
        .def_readwrite("fields", &FoFClass::fields)
        .def_readwrite("instruments", &FoFClass::instruments)
        .def_readwrite("exposures", &FoFClass::exposures)
        .def_readwrite("allPoints", &FoFClass::allPoints)
        .def_readwrite("sequence", &FoFClass::sequence)
        .def_readwrite("extn", &FoFClass::extn)
        .def_readwrite("obj", &FoFClass::obj)
        .def_readwrite("matchRadius", &FoFClass::matchRadius)
        .def_readwrite("a", &FoFClass::a)
        .def_readwrite("b", &FoFClass::b)
        .def_readwrite("c", &FoFClass::c)
        .def("addTest", &FoFClass::addTest)
        .def("addCatalog", &FoFClass::addCatalog)
        .def("sortMatches", &FoFClass::sortMatches)
        .def("writeMatches", &FoFClass::writeMatches);

}
