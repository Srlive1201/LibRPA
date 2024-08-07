#include "vector3_order.h"

template<>
bool operator< ( const Vector3_Order<double> &v1, const Vector3_Order<double> &v2 )
{
	constexpr double threshold = 1E-8;
	if      ( v1.x < v2.x - threshold ) return true;
	else if ( v1.x > v2.x + threshold ) return false;
	if      ( v1.y < v2.y - threshold ) return true;
	else if ( v1.y > v2.y + threshold ) return false;
	if      ( v1.z < v2.z - threshold ) return true;
	else if ( v1.z > v2.z + threshold ) return false;
	return false;
}

template<>
bool operator< ( const Vector3_Order<int> &v1, const Vector3_Order<int> &v2 )
{
	if      ( v1.x < v2.x ) return true;
	else if ( v1.x > v2.x ) return false;
	if      ( v1.y < v2.y ) return true;
	else if ( v1.y > v2.y ) return false;
	if      ( v1.z < v2.z ) return true;
	else if ( v1.z > v2.z ) return false;
	return false;
}

template<>
bool operator== ( const Vector3_Order<double> &v1, const Vector3_Order<double> &v2 )
{
	constexpr double threshold = 1E-8;
	if ( std::abs(v1.x - v2.x) > threshold ) return false;
	if ( std::abs(v1.y - v2.y) > threshold ) return false;
	if ( std::abs(v1.z - v2.z) > threshold ) return false;
	return true;
}

template<>
bool operator== ( const Vector3_Order<int> &v1, const Vector3_Order<int> &v2 )
{
    return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
}
