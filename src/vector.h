/***************************************************************************
 *   Copyright (C) 2009-2012 by Veselin Georgiev, Slavomir Kaslev et al    *
 *   admin@raytracing-bg.net                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef __VECTOR3D_H__
#define __VECTOR3D_H__

#include <stdio.h>
#include <math.h>

struct Vector {
	double x, y, z;
	
	Vector () {}
	Vector(double _x, double _y, double _z) { set(_x, _y, _z); }
	void set(double _x, double _y, double _z)
	{
		x = _x;
		y = _y;
		z = _z;
	}
	void makeZero(void)
	{
		x = y = z = 0.0;
	}
	inline double length(void) const
	{
		return sqrt(x * x + y * y + z * z);
	}
	inline double lengthSqr(void) const
	{
		return (x * x + y * y + z * z);
	}
	inline void scale(double multiplier)
	{
		x *= multiplier;
		y *= multiplier;
		z *= multiplier;
	}
	void operator *= (double multiplier)
	{
		scale(multiplier);
	}
	void operator /= (double divider)
	{
		scale(1.0 / divider);
	}
	inline void normalize(void)
	{
		double lSqr = (x * x + y * y + z * z);
		if (lSqr >= 0.999999 && lSqr < 1.000001) return;
		scale(1.0 / sqrt(lSqr));
	}
	void setLength(double newLength)
	{
		scale(newLength / length());
	}
	void print() const { printf("(%.3lf, %.3lf, %.3lf)", x, y, z); }
	void println() const { printf("(%.3lf, %.3lf, %.3lf)\n", x, y, z); }
	inline double& operator[] (int i)
	{
		return (&x)[i];
	}
	inline const double& operator[] (int i) const
	{
		return (&x)[i];
	}
	int maxDimension() const
	{
		int bi = 0;
		double maxD = fabs(x);
		if (fabs(y) > maxD) { maxD = fabs(y); bi = 1; }
		if (fabs(z) > maxD) { maxD = fabs(z); bi = 2; }
		return bi;
	}
	void operator += (const Vector& rhs) { x += rhs.x; y += rhs.y; z += rhs.z; }
	void operator -= (const Vector& rhs) { x -= rhs.x; y -= rhs.y; z -= rhs.z; }
};

inline Vector operator + (const Vector& a, const Vector& b)
{
	return Vector(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline Vector operator - (const Vector& a, const Vector& b)
{
	return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline Vector operator - (const Vector& a)
{
	return Vector(-a.x, -a.y, -a.z);
}

/// dot product
inline double operator * (const Vector& a, const Vector& b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
/// dot product (functional form, to make it more explicit):
inline double dot(const Vector& a, const Vector& b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
/// cross product
inline Vector operator ^ (const Vector& a, const Vector& b)
{
	return Vector(
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x
	);
}

inline Vector operator * (const Vector& a, double multiplier)
{
	return Vector(a.x * multiplier, a.y * multiplier, a.z * multiplier);
}
inline Vector operator * (double multiplier, const Vector& a)
{
	return Vector(a.x * multiplier, a.y * multiplier, a.z * multiplier);
}
inline Vector operator / (const Vector& a, double divider)
{
	double multiplier = 1.0 / divider;
	return Vector(a.x * multiplier, a.y * multiplier, a.z * multiplier);
}

/// given an unit vector a, create an orhonormed system (a, b, c). Code is deterministic.
inline void orthonormedSystem(const Vector& a, Vector& b, Vector& c)
{
	Vector temp = Vector(1, 0, 0);
	if (fabs(dot(a, temp)) > 0.99) {
		temp = Vector(0, 1, 0);
		if (fabs(dot(a, temp)) > 0.99)
			temp = Vector(0, 0, 1);
	}
	b = a ^ temp;
	b.normalize();
	c = a ^ b;
	c.normalize();
}

inline Vector reflect(const Vector& toBeReflected, const Vector& normal)
{
	return toBeReflected + 2 * (dot(normal, -toBeReflected)) * normal;
}

inline Vector refract(const Vector& i, const Vector& n, float ior)
{
	float NdotI = float(dot(i, n));
	float k = 1 - (ior * ior) * (1 - NdotI * NdotI);
	if (k < 0)
		return Vector(0, 0, 0);
	return ior * i - (ior * NdotI + sqrt(k)) * n;
}

inline Vector faceforward(const Vector& v, const Vector& right)
{
	if (dot(right, v) < 0) return v; else return -v;
}

enum RayFlags {
	RF_SHADOW = 0x0001,
	RF_GLOSSY = 0x0002,
	RF_DIFFUSE = 0x0004,
	RF_DEBUG  = 0x0100,
};

struct Ray {
	Vector start, dir;
	int flags;
	int depth;
	Ray() { flags = false; depth = 0; }
	Ray(const Vector& _start, const Vector& _dir) {
		start = _start;
		dir = _dir;
		depth = 0;
		flags = 0;
	}
};

struct FastRay: public Ray {
	Vector rdir;
	FastRay(const Ray& ray): Ray(ray)
	{
		rdir.x = fabs(dir.x) > 1e-9 ? 1.0 / dir.x : 0.0;
		rdir.y = fabs(dir.y) > 1e-9 ? 1.0 / dir.y : 0.0;
		rdir.z = fabs(dir.z) > 1e-9 ? 1.0 / dir.z : 0.0;
	}
};

#endif // __VECTOR3D_H__
