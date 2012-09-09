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
#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "constants.h"
#include <algorithm>
using std::min;
using std::max;

inline double signOf(double x) { return x > 0 ? +1 : -1; }
template <typename T>
inline T sqr(T a) { return a * a; }
inline double toRadians(double angle) { return angle / 180.0 * PI; }
inline double toDegrees(double angle_rad) { return angle_rad / PI * 180.0; }
inline int nearestInt(float x) { return (int) floor(x + 0.5f); }
template <typename T>
inline T abs(T a) { return a < 0 ? -a : a; }

static void findEmptyFN(char fn[], const char* filenamePrefix, const char* filenameSuffix)
{
	for (int i = 0; i < 10000; i++) {
		sprintf(fn, "%s%04d%s", filenamePrefix, i, filenameSuffix);
		FILE* f = fopen(fn, "rb");
		if (!f) return; // file doesn't exist - use that
		fclose(f);
	}
}

class FileRAII {
	FILE* fp;
public:
	FileRAII(FILE *f): fp(f) {}
	~FileRAII() { if (fp) fclose(fp); fp = NULL; }
};


#endif // __UTIL_H__
