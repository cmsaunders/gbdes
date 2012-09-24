//  $Id: FITStypes.h,v 1.7 2009/11/02 22:48:53 garyb Exp $
// Definitions, includes, and constants for FITS files and CFITSIO.
#ifndef FITSTYPES_H
#define FITSTYPES_H
#include <typeinfo>

//namespace, but then std libraries get stuck in there too...
namespace FITS {
#include "fitsio.h"	// ??? would be nice to keep this inside FITS


  //const int MAX_FITS_FILES_OPEN=NIOBUF;
  const int MAX_FITS_FILES_OPEN=40;
  const int MAX_COLUMN_NAME_LENGTH=FLEN_VALUE;
  
  // Assumption is that intrinsic int type will be 4 or 8 bits
  // and intrinsic long type could also be either 4 or 8 bits.
  // *** no longer worrying about 2-bit intrinsic ints ***
  // For FITS images, LONGLONG is 8 byte, LONG is 4 byte.
  // Note that CFITSIO has a typedef for LONGLONG
  const bool CLongIsFITSLongLong = (sizeof(long)*8==LONGLONG_IMG);
  //  (else it is assumed to be FITS LONG, 4 bytes)
  const bool CIntIsFITSLong = (sizeof(int)*8==LONG_IMG);
  //  (else assume that int is FITS LONGLONG)

  class FITSError: public std::runtime_error {
  public: 
    FITSError(const string &m=""): std::runtime_error("FITS Error: " + m) {}
  };
  class FITSCantOpen: public FITSError {
  public:
    FITSCantOpen(const string& fname, const string& m=""): 
      FITSError("Cannot open FITS file " + fname + ": " + m) {}
  };

  // Make enumerators out of the CFITSIO constants we'll need
  class Flags {
  private:
    int f;
  public:
    explicit Flags(int i): f(i) {};
    Flags(const Flags& rhs): f(rhs.f) {};
    const Flags& operator=(const Flags& rhs) {f=rhs.f; return *this;}
    Flags operator+(const Flags rhs) const {return Flags(f+rhs.f);}
    Flags operator|(const Flags rhs) const {return Flags(f|rhs.f);}
    Flags operator&(const Flags rhs) const {return Flags(f&rhs.f);}
    bool operator==(const Flags rhs) const {return f == rhs.f;}
    bool operator!=(const Flags rhs) const {return f != rhs.f;}
    Flags operator~() const {return Flags(~f);}
    operator bool() const {return f!=0;}
    int getInt() const {return f;}
  };

  const Flags ReadOnly(0);
  const Flags ReadWrite(1);
  const Flags Overwrite(2);
  const Flags CreateImage(4);

  enum HDUType {HDUImage=IMAGE_HDU, 
		HDUAsciiTable=ASCII_TBL,
		HDUBinTable=BINARY_TBL,
		HDUAny = ANY_HDU };
  // CFITSIO data types: note that these sometimes are used to describe which C intrinsic
  // type is in use (in which case names have obvious relation to intrinsics) and sometimes they
  // describe a storage format within FITS files (in which case they need some machine-dependent
  // translation to intrinsic types).
  enum DataType {Tnull, 
		 Tbit = TBIT, 
		 Tbyte = TBYTE, 
		 Tsbyte = TSBYTE, 
		 Tlogical = TLOGICAL, 
		 Tstring = TSTRING, 
		 Tushort = TUSHORT, 
		 Tshort = TSHORT,
		 Tuint = TUINT,
		 Tint = TINT,
		 Tulong = TULONG,
		 Tlong = TLONG,
		 Tlonglong = TLONGLONG,
		 Tfloat = TFLOAT,
		 Tdouble = TDOUBLE,
		 Tcomplex = TCOMPLEX,
		 Tdblcomplex=TDBLCOMPLEX,
		 Tint32bit=TINT32BIT};

  // Handy function stolen from CCFits, uses RTTI to give back the 
  // enum code of FITS datatype for a given native data type.
  // fitsio.h gives a typedef for LONGLONG
  // Return Tnull if template type has no FITS equivalent.
  template <typename T>
  inline DataType FITSTypeOf() {
    // Note ordering: since integer types can be redundant, prefer int to long to LONGLONG.
    if (typeid(T) == typeid(int))		return Tint;
    if (typeid(T) == typeid(unsigned int))	return Tuint;
    if (typeid(T) == typeid(long))		return Tlong;
    if (typeid(T) == typeid(unsigned long))	return Tulong;
    if (typeid(T) == typeid(LONGLONG))		return Tlonglong;
    return Tnull;
  }
  template <>
  inline DataType FITSTypeOf<bool>() {return Tlogical;}
  template <>
  inline DataType FITSTypeOf<signed char>() {return Tsbyte;}
  template <>
  inline DataType FITSTypeOf<unsigned char>() {return Tbyte;}
  template <>
  inline DataType FITSTypeOf<short>() {return Tshort;}
  template <>
  inline DataType FITSTypeOf<unsigned short>() {return Tushort;}
  template <>
  inline DataType FITSTypeOf<int>() {return Tint;}
  template <>
  inline DataType FITSTypeOf<unsigned int>() {return Tuint;}
  template <>
  inline DataType FITSTypeOf<float>() {return Tfloat;}
  template <>
  inline DataType FITSTypeOf<double>() {return Tdouble;}
  template <>
  inline DataType FITSTypeOf<complex<float> >() {return Tcomplex;}
  template <>
  inline DataType FITSTypeOf<complex<double> >() {return Tdblcomplex;}
  template <>
  inline DataType FITSTypeOf<string>() {return Tstring;}

  // Convert to/from BITPIX keywords to code for intrinsic data type that can hold it
  inline DataType
  Bitpix_to_DataType(const int bitpix) {
    if (bitpix==BYTE_IMG)   return Tbyte;
    if (bitpix==SHORT_IMG)  return Tshort;
    if (bitpix==LONG_IMG)   return CIntIsFITSLong ? Tint : (CLongIsFITSLongLong? Tnull : Tlong);
    if (bitpix==FLOAT_IMG)  return Tfloat;
    if (bitpix==DOUBLE_IMG) return Tdouble;
    if (bitpix==USHORT_IMG) return Tushort;
    if (bitpix==ULONG_IMG)  return CIntIsFITSLong ? Tuint : (CLongIsFITSLongLong? Tnull : Tulong);
    if (bitpix==LONGLONG_IMG) return Tlonglong;
    throw FITSError("Unknown BITPIX value");
  }
  inline int
  DataType_to_Bitpix(const DataType dt) {
    if (dt==Tbyte)	return BYTE_IMG;
    if (dt==Tshort)	return SHORT_IMG;
    if (dt==Tint)	return CIntIsFITSLong ? LONG_IMG : LONGLONG_IMG;
    if (dt==Tlong)	return CLongIsFITSLongLong? LONGLONG_IMG : LONG_IMG;
    if (dt==Tlonglong)	return LONGLONG_IMG;
    if (dt==Tfloat)	return FLOAT_IMG;
    if (dt==Tdouble)	return DOUBLE_IMG;

    if (dt==Tushort)	return USHORT_IMG;
    // Is no unsigned 8-bit format, so just fudge using signed image:
    if (dt==Tuint)	return CIntIsFITSLong ? ULONG_IMG : LONGLONG_IMG;
    if (dt==Tulong)	return CLongIsFITSLongLong? LONGLONG_IMG : ULONG_IMG;
    throw FITSError("Datatype cannot be converted to BITPIX");
  }

}  //namespace FITS

#endif  //FITSTYPES_H
