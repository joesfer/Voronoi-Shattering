#ifndef MPlane_h__
#define MPlane_h__

#if ( MAYA_API_VERSION == 200806 )
	// maya 2008

#include <maya/MPoint.h>

class	MPlane {
public:
	MPlane() : a(0), b(0), c(0), d(0) {}
	void setPlane( double A, double B, double C, double D ) { a = A; b = B; c = C; d = D; }
	double directedDistance( const MPoint& p ) const { return a * p.x + b * p.y + c * p.z - d; }
	const double distance( const MPoint& p ) const { return fabs( directedDistance( p ) ); } 
	MVector normal() const { return MVector( a, b, c ); }
private:
	double a,b,c,d;
};

#else 
	// maya 2011
	#include <maya/MPlane.h>	
#endif

#endif // MPlane_h__
