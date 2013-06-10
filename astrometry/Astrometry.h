// Routines for the manipulation of astronomical coordinates and times.
#ifndef ASTROMETRY_H
#define ASTROMETRY_H

#include "Std.h"
#include "UseTMV.h"
#include "AstronomicalConstants.h"
#include <iostream>
#include <string>
#include <stdexcept>

#include "TMV_Small.h"
typedef tmv::SmallVector<double,2> Vector2;
typedef tmv::SmallVector<double,3> Vector3;
typedef tmv::SmallMatrix<double,2,2> Matrix22;
typedef tmv::SmallMatrix<double,2,3> Matrix23;
typedef tmv::SmallMatrix<double,3,2> Matrix32;
typedef tmv::SmallMatrix<double,3,3> Matrix33;

namespace astrometry {

  class AstrometryError: public std::runtime_error {
  public:
    AstrometryError(const string &m=""): 
      std::runtime_error("Astrometry Error: " +m) {}
  };

  // ??? Take this out later:
  string  degdms(double degr, int decimalPlaces);
  string  deghms(double degr, int decimalPlaces);
  double  hmsdeg(const string& s);
  double  dmsdeg(const string& s);

  // Specify UT Date + Time
  // Some details here:  it is assumed that the time stored here is UTC.
  // The time used in the JPL ephemeris is TDB (barycentric dynamical), which
  // is within 0.002s of TT (Terrestrial time), which is a coordinate time at
  // Earth's geoid.  TT is a constant (~34s) off from TAI=atomic time.
  // UTC has leap seconds to keep it close to the Earth rotation time, so
  // TT-UTC is a function of time.
  // Subtracting/adding times is done in such a way as to account for leap
  // seconds, though I do not gaurantee success if one of the times is within
  // a minute or so of the leap second itself.  The list of leap seconds is maintained
  // in the file LeapSeconds.h as the function TAIminusUTC().
  // JD's are stored as doubles, which have 52 bits mantissa, so JD should
  // be accurate to 0.1 ms or better, which is <(TDB-TT).
  class UT {
  private:
    // Flags to say whether jd or y/m/d representation is valid
    double jd;
    mutable bool ymdValid;
    mutable int y;
    mutable int mo;
    mutable double d;
    void buildJD();
    void buildYMD() const;
  public:
    // ??? Add validation of date specs to below???
    UT(): jd(0.), ymdValid(false) {}
    UT(double jd_) {set(jd_);}
    UT(int y_, int m_, double d_) {set(y_,m_,d_);}
    UT(int y_, int m_, int d_, int h_, int min_, double s_) {
      set(y_,m_,d_,h_,m_,s_);
    }
    void set(double jd_);
    void set(int y_, int m_, double d_) ;
    void set(int y_, int m_, int d_, int h_, int min_, double s_) ;

    // Note that all the set/get are based on DAY units to match JD's
    void setMJD(double mjd) {set(mjd+MJD0);}
    void setTT(double tt) {set(tt-TTminusUTC(tt)/DAY);}

    double getJD() const {return jd;}
    double getMJD() const {return getJD()-MJD0;}
    double getTT() const {return getJD()+TTminusUTC(getJD())/DAY;}
    void getYMD(int& y_, int& m_, double& d_) const {
      buildYMD(); 
      y_=y; m_=mo; d_=d;
    }
    void getYMDHMS(int& y_, int& m_, int& d, 
		   int& h, int& mn, double& s) const;


    // All of the time intervals use units from AstronomicalConstants.h,
    // which are YEAR.
    static double TTminusUTC(double jd_) ;
    static double TAIminusUTC(double jd_);
    const UT& operator+=(double dt);
    const UT& operator-=(double dt);
    double operator-(const UT& rhs) const;

    void writeYMD(std::ostream& os) const;
    void writeYMDHMS(std::ostream& os) const;
    // ??? writeMPC(std::ostream& os);
    void read(std::istream& is);
    friend std::ostream& operator<<(std::ostream& os, const UT& rhs);
    friend std::istream& operator>>(std::istream& is, UT& rhs);
  };

  ///////////////////////////////////////////////////////////////////////////
  // Spherical coordinate systems
  // The internal format for spherical coordinates is a 3d unit
  // vector, to simplify conversions.  The class's external appearance will be
  // as some 2d representation, either lat/lon or projected coordinates.  
  // The derivs[To|From]3d matrix will give partial derivs of the unit vectors w.r.t.
  // the nominal 2d coordinates or vice-versa
  //
  // Coordinate conversions occur simply by casting one kind of SphericalCoordinate
  // to another.  Internally this is handled by every derived type being able
  // to convert its 3d unit vector to/from the ICRS system.
  ///////////////////////////////////////////////////////////////////////////

  class CartesianCoords;
  class CartesianICRS;
  class CartesianEcliptic;
  class CartesianInvariable;
  class CartesianCustom;

  class SphericalCoords {
  protected:
    Vector3 x;
    // partials in these equations are 3x3 unit vector xforms:
    virtual Vector3 convertToICRS(Matrix33* partials=0) const=0;
    virtual void convertFromICRS(const Vector3& xeq, 
			       Matrix33* partials=0) =0;
    SphericalCoords(const Vector3& x_);
    // Utility for projecting Cartesian vectors to their spherical coords
    void projectIt(const CartesianCoords& rhs,
		   Matrix33* partials=0);
    void projectIt(const CartesianCoords& rhs,
		   Matrix23* partials);
  public:
    SphericalCoords();
    SphericalCoords(double lon, double lat);
    virtual ~SphericalCoords() {};
    virtual SphericalCoords* duplicate() const=0;

    virtual void write(std::ostream& os, int decimalPlaces=3) const;
    virtual void read(std::istream& is);
    friend std::ostream& operator<<(std::ostream& os, 
				    const SphericalCoords& rhs);
    friend std::istream& operator>>(std::istream& is, SphericalCoords& rhs);

    Vector3 getUnitVector() const {return x;}
    void setUnitVector(const Vector3& x_);
    Vector2 getLonLat() const;
    virtual void getLonLat(double& lon, double& lat) const;
    // Derivatives between the 3d unit vectors and the 2d coords of
    // this representation (e.g. lon/lat):
    // Assumes that radius is fixed at unity in 3d.
    virtual Matrix23 derivsTo2d() const;
    virtual Matrix32 derivsFrom2d() const;

    virtual void setLonLat(double lon, double lat);
    virtual void setLonLat(const Vector2& x_);

    // Generic conversion from one coord system to another.
    // The dimension of the partials matrix tells the routine whether
    // the desired derivates are to/from the unit vectors (dim=3) or
    // to/from the nominal angles (dim=2).
    void convertFrom(const SphericalCoords& rhs);
    void convertFrom(const SphericalCoords& rhs, Matrix33& partials);
    void convertFrom(const SphericalCoords& rhs, Matrix32& partials);
    void convertFrom(const SphericalCoords& rhs, Matrix23& partials);
    void convertFrom(const SphericalCoords& rhs, Matrix22& partials);

    // Return great-circle distance (radians) from another point
    double distance(const SphericalCoords& rhs) const;
  };

  std::ostream& operator<<(std::ostream& os, const SphericalCoords& rhs);
  std::istream& operator>>(std::istream& is, SphericalCoords& rhs);

  class SphericalICRS: public SphericalCoords {
  protected:
    virtual Vector3 convertToICRS(Matrix33* partials=0) const;
    virtual void convertFromICRS(const Vector3& xeq, Matrix33* partials=0);
  public:
    SphericalICRS() {};
    SphericalICRS(double ra, double dec): SphericalCoords(ra,dec) {};
    virtual SphericalCoords* duplicate() const {
      return new SphericalICRS(*this);
    }
    // Conversions & projection: 
    explicit SphericalICRS(const SphericalCoords& rhs);
    explicit SphericalICRS(const SphericalCoords& rhs, Matrix33& partials);
    explicit SphericalICRS(const CartesianICRS& rhs);
    explicit SphericalICRS(const CartesianICRS& rhs, Matrix33& partials);
    // This one gives partials of lat/lon w.r.t. input 3d vector:
    explicit SphericalICRS(const CartesianICRS& rhs, Matrix23& partials);
    ~SphericalICRS() {};
    // Here are aliases that give usual ra/dec names
    void getRADec(double& ra, double& dec) const {getLonLat(ra, dec);}
    void setRADec(double ra, double dec) {setLonLat(ra,dec);}
  };

  class SphericalEcliptic: public SphericalCoords {
  protected:
    virtual Vector3 convertToICRS(Matrix33* partials=0) const;
    virtual void convertFromICRS(const Vector3& xeq, Matrix33* partials=0);
  public:
    SphericalEcliptic() {};
    SphericalEcliptic(double lon, double lat): SphericalCoords(lon,lat) {}
    ~SphericalEcliptic() {};
    virtual SphericalCoords* duplicate() const {
      return new SphericalEcliptic(*this);
    }
    // Conversions & projection: 
    explicit SphericalEcliptic(const SphericalCoords& rhs);
    explicit SphericalEcliptic(const SphericalCoords& rhs, Matrix33& partials);
    explicit SphericalEcliptic(const CartesianEcliptic& rhs);
    explicit SphericalEcliptic(const CartesianEcliptic& rhs, Matrix33& partials);
  };

  class SphericalInvariable: public SphericalCoords {
  protected:
    virtual Vector3 convertToICRS(Matrix33* partials=0) const;
    virtual void convertFromICRS(const Vector3& xeq, Matrix33* partials=0);
  public:
    SphericalInvariable() {};
    SphericalInvariable(double lon, double lat): SphericalCoords(lon,lat) {}
    ~SphericalInvariable() {};
    virtual SphericalCoords* duplicate() const {
      return new SphericalInvariable(*this);
    }
    // Conversions & projection: 
    explicit SphericalInvariable(const SphericalCoords& rhs);
    explicit SphericalInvariable(const SphericalCoords& rhs, Matrix33& partials);
    explicit SphericalInvariable(const CartesianInvariable& rhs);
    explicit SphericalInvariable(const CartesianInvariable& rhs, 
				 Matrix33& partials);
  };

  // Define an arbitrary reference pole and orientation
  class Orientation {
  public:
    Orientation(const SphericalCoords& pole_, double pa=0.);
    Orientation(): pole(), zrot(0.) {buildMatrix();}
    const Orientation& operator=(const Orientation& rhs) {
      pole.convertFrom(rhs.pole); 
      zrot=rhs.zrot; rMatrix=rhs.rMatrix; 
      return *this;
    }
    const Matrix33& m() const {return rMatrix;}

    void set(const SphericalCoords& pole_, double pa);
    SphericalICRS getPole() const {return pole;}
    void setPole(const SphericalCoords& pole_);
    double  getZRot() const {return zrot;}
    double  getPA() const {return -getZRot();}
    void setZRot(double zrot_);

    // Set y axis at some N->E PA in the given system:
    // (in the E=x, N=y system, note pos PA is a negative angle
    // in the right-handed sense.)
    void alignToEcliptic(double phi=0.);
    void alignToICRS(double phi=0.);
    void alignToInvariable(double phi=0.);

    Vector3 fromICRS(const Vector3& x) const;
    Vector3 toICRS(const Vector3& x) const;

    void write(std::ostream& os) const;
    void read(std::istream& is);
  private:
    SphericalICRS pole;
    double        zrot;
    Matrix33     rMatrix;
    void	  buildMatrix();
  };

  std::ostream& operator<<(std::ostream& os, const Orientation& o);
  std::istream& operator>>(std::istream& is, Orientation& o);


  // The SphericalCustom and Gnomonic are defined
  // by an Orientation object giving its pole and the position angle of the coord
  // system about this pole.
  //
  // On creation, if shareOrient=false, a private duplicate of the Orientation
  // will be created (and destroyed along with this object).
  // But means copying 2 coords + orient=(2 coords + rot + 3x3 matrix) each
  // time a coordinate is made, which may be slow for processing many coordinates.
  //
  // Use shareOrient=true if you are creating many SphericalCustom objects with 
  // the same orientation and want to reduce copying overhead.  In this case,
  // a pointer to the input Orientation is stored.  You are then responsible for
  // keeping the input Orientation in existence as long as the SphericalCustom is
  // in use, then deleting the Orientation.

  // SphericalCustom and Gnomonic share everything except lat/lon calculation so
  // both will be derived from SphericalCustomBase
  class SphericalCustomBase: public SphericalCoords {
  private:
    const Orientation* orient;
    bool ownOrient;
    virtual Vector3 convertToICRS(Matrix33* partials=0) const;
    virtual void convertFromICRS(const Vector3& xeq, Matrix33* partials=0);
  public:
    SphericalCustomBase(const Orientation& o, bool shareOrient);
    SphericalCustomBase(double xi, double eta, const Orientation& o, bool shareOrient) ;
    virtual ~SphericalCustomBase() {if (ownOrient) delete orient;}
    explicit SphericalCustomBase(const SphericalCustomBase& rhs);
    SphericalCustomBase& operator=(const SphericalCustomBase& rhs);

    // Conversions & projection: 
    explicit SphericalCustomBase(const SphericalCoords& rhs,
				 const Orientation& o,
				 bool shareOrient);
    explicit SphericalCustomBase(const SphericalCoords& rhs, 
				 const Orientation& o,
				 Matrix33& partials,
				 bool shareOrient);
    // For now, conversion from CartesianCustom will share the ReferenceFrame's Orientation:
    explicit SphericalCustomBase(const CartesianCustom& rhs);
    explicit SphericalCustomBase(const CartesianCustom& rhs, Matrix33& partials);

    const Orientation* getOrient() const {return orient;}
  };

  class SphericalCustom: public SphericalCustomBase {
  public:
    SphericalCustom(const Orientation& o, bool shareOrient=false):
      SphericalCustomBase(o, shareOrient) {}
    SphericalCustom(double xi, double eta, const Orientation& o, bool shareOrient=false):
      SphericalCustomBase(xi, eta, o, shareOrient) {}
    virtual ~SphericalCustom() {}
    SphericalCustom(const SphericalCustom& rhs): SphericalCustomBase(rhs) {}
    SphericalCustom& operator=(const SphericalCustom& rhs) {
      SphericalCustomBase::operator=(rhs);
      return *this;
    }
    virtual SphericalCoords* duplicate() const {
      return new SphericalCustom(*this);
    }

    // Conversions & projection: 
    explicit SphericalCustom(const SphericalCoords& rhs,
			     const Orientation& o,
			     bool shareOrient = false): 
      SphericalCustomBase(rhs, o, shareOrient) {}
    explicit SphericalCustom(const SphericalCoords& rhs, 
			     const Orientation& o,
			     Matrix33& partials,
			     bool shareOrient = false): 
      SphericalCustomBase(rhs, o, partials, shareOrient) {}
    // For now, conversion from CartesianCustom will share the ReferenceFrame's Orientation:
    explicit SphericalCustom(const CartesianCustom& rhs):
      SphericalCustomBase(rhs) {}
    explicit SphericalCustom(const CartesianCustom& rhs, Matrix33& partials):
      SphericalCustomBase(rhs) {}
  };

  // The Gnomonic class is gnomonic projection at specified orientation.
  // Differs from SphericalCustom only in projecting from center of sphere
  // onto the tangent plane to get the x/y (=lon/lat) coordinates.
  // So this class overrides only the Lon/Lat-related methods of SphericalCustomBase.
  // Gnomonic project is valid only on the hemisphere nearest the pole.
  class Gnomonic: public SphericalCustomBase {
  public:
    Gnomonic(const Orientation& o, bool shareOrient=false):
      SphericalCustomBase(o, shareOrient) {}
    Gnomonic(double xi, double eta, const Orientation& o, bool shareOrient=false):
      SphericalCustomBase(xi, eta, o, shareOrient) {}
    virtual ~Gnomonic() {}
    Gnomonic(const Gnomonic& rhs): SphericalCustomBase(rhs) {}
    Gnomonic& operator=(const Gnomonic& rhs) {
      SphericalCustomBase::operator=(rhs);
      return *this;
    }
    virtual SphericalCoords* duplicate() const {
      return new Gnomonic(*this);
    }

    // Conversions & projection: 
    explicit Gnomonic(const SphericalCoords& rhs,
		      const Orientation& o,
		      bool shareOrient = false): 
      SphericalCustomBase(rhs, o, shareOrient) {}
    explicit Gnomonic(const SphericalCoords& rhs, 
		      const Orientation& o,
		      Matrix33& partials,
		      bool shareOrient = false): 
      SphericalCustomBase(rhs, o, partials, shareOrient) {}
    // For now, conversion from CartesianCustom will share the ReferenceFrame's Orientation:
    explicit Gnomonic(const CartesianCustom& rhs):
      SphericalCustomBase(rhs) {}
    explicit Gnomonic(const CartesianCustom& rhs, Matrix33& partials):
      SphericalCustomBase(rhs) {}

    // Since TangentPlane is not a lat/lon system, need to override
    // some operations:
    // I/O will give xi/eta in dms
    void write(std::ostream& os, int decimalPlaces=3) const;
    virtual void read(std::istream& is);  
    // and "lonlat" is really xi/eta, so 2d conversions differ:
    void getLonLat(double& lon, double& lat) const;
    virtual void setLonLat(double lon, double lat);
    virtual void setLonLat(const Vector2& x_) {setLonLat(x_[0], x_[1]);}
    virtual Matrix23 derivsTo2d() const;
    virtual Matrix32 derivsFrom2d() const;
  };


  //////////////////////////////////////////////////////////
  // Spatial and angular coordinate systems
  //////////////////////////////////////////////////////////

  // Base class for 3d Cartesian systems
  class CartesianCoords {
  protected:
    Vector3 x;
    void plusEq(const CartesianCoords& rhs) {x+=rhs.x;}
    void minusEq(const CartesianCoords& rhs) {x-=rhs.x;}
    // partials in these equations are 3x3 unit vector xforms:
    virtual Vector3 convertToICRS(Matrix33* partials=0) const=0;
    virtual void convertFromICRS(const Vector3& xeq, 
			       Matrix33* partials=0) =0;
  public:
    CartesianCoords(): x(0.) {}
    CartesianCoords(double x_, double y_, double z_) {
      x[0]=x_; x[1]=y_; x[2]=z_;
    }
    explicit CartesianCoords(const Vector3& x_): x(x_) {}

    // Generic conversion from one coord system to another.
    void convertFrom(const CartesianCoords& rhs);
    void convertFrom(const CartesianCoords& rhs, Matrix33& partials);

    double operator[](int i) const {return x[i];}
    double& operator[](int i) {return x[i];}
    bool operator==(const CartesianCoords& rhs) const {
      return x==rhs.x;
    }

    void setVector(const Vector3 x_) {x=x_;}
    Vector3 getVector() const {return x;}

    const CartesianCoords& operator=(const Vector3& x_) {
      x=x_; return *this;
    }
    const CartesianCoords& operator+=(const Vector3& x_) {
      x+=x_; return *this;
    }
    const CartesianCoords& operator-=(const Vector3& x_) {
      x-=x_; return *this;
    }
    const CartesianCoords& operator*=(double scalar) {
      x*=scalar; return *this;
    }

    double radius() const {return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);}
    void write(std::ostream& os) const;
    void read(std::istream& is);
    friend std::ostream& operator<<(std::ostream& os, 
				    const CartesianCoords& rhs);
    friend std::istream& operator>>(std::istream& is, CartesianCoords& rhs);
  };

  std::ostream& operator<<(std::ostream& os, const CartesianCoords& rhs);
  std::istream& operator>>(std::istream& is, CartesianCoords& rhs);
  
  class CartesianICRS: public CartesianCoords {
  private:
    //  CartesianCoords(rhs) {}
    virtual Vector3 convertToICRS(Matrix33* partials=0) const;
    virtual void convertFromICRS(const Vector3& xeq, Matrix33* partials=0);
  public:
    CartesianICRS() {}
    CartesianICRS(const Vector3& x_): CartesianCoords(x_) {}

    // Constructor with coordinate transformation:
    explicit CartesianICRS(const CartesianCoords& rhs); 
    explicit CartesianICRS(const CartesianCoords& rhs, Matrix33& partials); 
    const CartesianICRS& operator=(const CartesianICRS& rhs) {
      x = rhs.x;
      return *this;
    }

    const CartesianICRS& operator+=(const CartesianICRS& rhs) {
      plusEq(rhs); return *this;
    }
    const CartesianICRS& operator-=(const CartesianICRS& rhs) {
      minusEq(rhs); return *this;
    }
  };

  class CartesianEcliptic: public CartesianCoords {
  private:
    virtual Vector3 convertToICRS(Matrix33* partials=0) const;
    virtual void convertFromICRS(const Vector3& xeq, Matrix33* partials=0);
  public:
    CartesianEcliptic() {}
    CartesianEcliptic(const Vector3& x_): CartesianCoords(x_) {}

    // Constructor with coordinate transformation:
    explicit CartesianEcliptic(const CartesianCoords& rhs); 
    explicit CartesianEcliptic(const CartesianCoords& rhs, Matrix33& partials); 

    const CartesianEcliptic& operator+=(const CartesianEcliptic& rhs) {
      plusEq(rhs); return *this;
    }
    const CartesianEcliptic& operator-=(const CartesianEcliptic& rhs) {
      minusEq(rhs); return *this;
    }
  };

  class CartesianInvariable: public CartesianCoords {
  private:
    virtual Vector3 convertToICRS(Matrix33* partials=0) const;
    virtual void convertFromICRS(const Vector3& xeq, Matrix33* partials=0);
  public:
    CartesianInvariable() {}
    CartesianInvariable(const Vector3& x_): CartesianCoords(x_) {}

    // Constructor with coordinate transformation:
    explicit CartesianInvariable(const CartesianCoords& rhs); 
    explicit CartesianInvariable(const CartesianCoords& rhs, 
				 Matrix33& partials); 

    const CartesianInvariable& operator+=(const CartesianInvariable& rhs) {
      plusEq(rhs); return *this;
    }
    const CartesianInvariable& operator-=(const CartesianInvariable& rhs) {
      minusEq(rhs); return *this;
    }
  };


  // Defines a full Cartesian system and associated Spherical system:
  class ReferenceFrame {
  public:
    CartesianICRS origin;
    Orientation   orient;
    
    ReferenceFrame(const CartesianCoords& origin_,
		   const Orientation& orient_): origin(origin_),
						orient(orient_) {}
    ReferenceFrame(const CartesianCoords& origin_,
		   const SphericalCoords& pole_, 
		   double zrot_=0.): origin(origin_),
				     orient(pole_, zrot_) {}
    ReferenceFrame(): origin(), orient() {}
    void write(std::ostream& os) const;
    void read(std::istream& is);
  };

  std::ostream& operator<<(std::ostream& os, const ReferenceFrame& rhs);
  std::istream& operator>>(std::istream& is, ReferenceFrame& rhs);

  // ??? Reference Frame ownership needs to be clarified
  class CartesianCustom: public CartesianCoords {
  private:
    const ReferenceFrame* ref;
    virtual Vector3 convertToICRS(Matrix33* partials=0) const;
    virtual void convertFromICRS(const Vector3& xeq, Matrix33* partials=0);
  public:
    CartesianCustom(const ReferenceFrame& r): ref(&r) {}
    CartesianCustom(const Vector3& x_,
		    const ReferenceFrame& r): ref(&r), CartesianCoords(x_) {}

    // Constructor with coordinate transformation:
    explicit CartesianCustom(const CartesianCoords& rhs,
			     const ReferenceFrame& r); 
    explicit CartesianCustom(const CartesianCoords& rhs, 
			     const ReferenceFrame& r, 
			     Matrix33& partials); 


    const ReferenceFrame* getFrame() const {return ref;}
    const CartesianCustom& operator+=(const CartesianCustom& rhs) {
      plusEq(rhs); return *this;
    }
    const CartesianCustom& operator-=(const CartesianCustom& rhs) {
      minusEq(rhs); return *this;
    }
  };

} //namespace astrometry

#endif //ASTROMETRY_H
