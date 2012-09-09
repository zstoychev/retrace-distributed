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
#ifndef __COLOR_H__
#define __COLOR_H__

#include "util.h"

inline unsigned convertTo8bit(float x, float gamma)
{
	if (x < 0) x = 0;
	if (x > 1) x = 1;
	if (gamma != 1.0f) x = pow(x, 1.0f / gamma); // output gamma correction
	return nearestInt(x * 255.0f);
}

/// Represents a color, using floatingpoint components in [0..1]
struct Color {
	float r, g, b;
	//
	Color() {}
	Color(float _r, float _g, float _b) //!< Construct a color from floatingpoint values
	{
		setColor(_r, _g, _b);
	}
	explicit Color(unsigned rgbcolor) //!< Construct a color from R8G8B8 value like "0xffce08"
	{
		b = (rgbcolor & 0xff) / 255.0f;
		g = ((rgbcolor >> 8) & 0xff) / 255.0f;
		r = ((rgbcolor >> 16) & 0xff) / 255.0f;
	}
	/// convert to RGB32, with channel shift specifications. The default values are for
	/// the blue channel occupying the least-significant byte
	unsigned toRGB32(int redShift = 16, int greenShift = 8, int blueShift = 0, float gamma = 1.0f) const
	{
		unsigned ir = convertTo8bit(r, gamma);
		unsigned ig = convertTo8bit(g, gamma);
		unsigned ib = convertTo8bit(b, gamma);
		return (ib << blueShift) | (ig << greenShift) | (ir << redShift);
	}
	/// make black
	void makeZero(void)
	{
		r = g = b = 0;
	}
	/// set the color explicitly
	void setColor(float _r, float _g, float _b)
	{
		r = _r;
		g = _g;
		b = _b;
	}
	/// get the intensity of the color (direct)
	inline float intensity(void) const
	{
		return (r + g + b) * 0.3333333f;
	}
	/// get the perceptual intensity of the color
	float intensityPerceptual(void) const
	{
		return (r * 0.299f + g * 0.587f + b * 0.114f);
	}
	/// Accumulates some color to the current
	void operator += (const Color& rhs)
	{
		r += rhs.r;
		g += rhs.g;
		b += rhs.b;
	}
	/// Multiplies some color to the current
	void operator *= (const Color& rhs)
	{
		r *= rhs.r;
		g *= rhs.g;
		b *= rhs.b;
	}
	/// multiplies the color
	void operator *= (float multiplier)
	{
		r *= multiplier;
		g *= multiplier;
		b *= multiplier;
	}
	/// divides the color
	void operator /= (float divider)
	{
		float mul = 1.0f / divider;
		r *= mul;
		g *= mul;
		b *= mul;
	}
	
	float& operator[] (const int x)
	{
		return (&r)[x];
	}
	
	const float& operator[] (const int x) const
	{
		return (&r)[x];
	}
	void clamp(void)
	{
		if (r < 0.0f) r = 0.0f;
		if (r > 1.0f) r = 1.0f;
		
		if (g < 0.0f) g = 0.0f;
		if (g > 1.0f) g = 1.0f;
		
		if (b < 0.0f) b = 0.0f;
		if (b > 1.0f) b = 1.0f;
	}
	void desaturate(float amount)
	{
		float mid = (r + g + b) / 3.0f;
		r = r * (1 - amount) + mid * amount;
		g = g * (1 - amount) + mid * amount;
		b = b * (1 - amount) + mid * amount;
	}
};

/// adds two colors
inline Color operator + (const Color& a, const Color& b)
{
	return Color(a.r + b.r, a.g + b.g, a.b + b.b);
}

/// subtracts two colors
inline Color operator - (const Color& a, const Color& b)
{
	return Color(a.r - b.r, a.g - b.g, a.b - b.b);
}

/// multiplies two colors
inline Color operator * (const Color& a, const Color& b)
{
	return Color(a.r * b.r, a.g * b.g, a.b * b.b);
}

/// multiplies a color by some multiplier
inline Color operator * (const Color& a, float multiplier)
{
	return Color(a.r * multiplier, a.g * multiplier, a.b * multiplier);
}

/// multiplies a color by some multiplier
inline Color operator * (float multiplier, const Color& a)
{
	return Color(a.r * multiplier, a.g * multiplier, a.b * multiplier);
}

/// divides some color
inline Color operator / (const Color& a, float divider)
{
	return Color(a.r / divider, a.g / divider, a.b / divider);
}

inline float colorDifference(const Color& a, const Color& b)
{
	return (fabs(a.r - b.r) + fabs(a.g - b.g) + fabs(a.b - b.b));
}

#endif // __COLOR_H__
