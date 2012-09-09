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
 
 /* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#include <math.h>
#include <SDL/SDL.h>
#include "random_generator.h"
#include "constants.h"

// various magic related to MT1
#define M 397
#define MATRIX_A 0x9908b0dfU   /* constant vector a */
#define UPPER_MASK 0x80000000U /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffU /* least significant r bits */

Random::Random(unsigned seed)
{
	this->seed(seed);
	nvi = 0; // no normal variates stored
}

void Random::seed(unsigned s)
{
	mt[0] = s;
	for (mti = 1; mti < MT_LEN; mti++) {
		mt[mti] = (1812433253U * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
	}
}

unsigned Random::_next(void)
{
	unsigned y;
	static const unsigned mag01[2]={ 0x0u, MATRIX_A };
	/* mag01[x] = x * MATRIX_A  for x=0,1 */
	
	if (mti >= MT_LEN) { /* generate N words at one time */
		int kk;
		
		for (kk=0; kk < MT_LEN - M; kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		for (; kk < MT_LEN - 1; kk++) {
			y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
			mt[kk] = mt[kk + (M - MT_LEN)] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		y = (mt[MT_LEN - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
		mt[MT_LEN - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1];
		mti = 0;
	}
	
	y = mt[mti++];
	
	/* Tempering */
	y ^= (y >> 11);
	y ^= (y << 7) & 0x9d2c5680U;
	y ^= (y << 15) & 0xefc60000U;
	y ^= (y >> 18);
	
	return y;
}

int Random::randint(int a, int b)
{
	if (a == b) return a;
	return a + (_next() % (b - a + 1));
}

float Random::randfloat(void)
{
	return float(_next()) * 0.00000000023283064365386962890625f; // 2^-32
}

double Random::randdouble(void)
{
	double z = double(_next()) * 0.00000000023283064365386962890625;
	z += _next();
	return z * 0.00000000023283064365386962890625; // 2^-32
	
}

double Random::gaussian(double mean, double sigma)
{
	// it's not easy to generate a normally-distributed random variable...
	// ... but it's much easier to generate two of them simultaneously
	// so we use buffering; generate two at once, store one in the buffer and return the other
	// on the next call, return the buffered one, and mark the cache as empty
	if (nvi) {
		nvi = 0;
		return nv_temp * sigma + mean;
	}
	nvi = 1;
	// Box-Muller transform:
	double u1 = 1.0 - randdouble();
	double u2 = 1.0 - randdouble();
	double s = sqrt(-2.0 * log(u1));
	double angle = 2 * PI * u2;
	double r1 = s * cos(angle);
	double r2 = s * sin(angle);
	nv_temp = r2;
	return r1 * sigma + mean;
}

void Random::unitDiscSample(double &x, double &y)
{
	double angle = randdouble() * 2 * PI;
	double rad = sqrt(randdouble());
	x = sin(angle) * rad;
	y = cos(angle) * rad;
}

struct HashMapEntry {
	Random r;
	unsigned key;
	char fill[128]; // skip to the next cacheline
};

const int RGENS = 257; // 257 is a prime number
static HashMapEntry rg_table[RGENS];

void initRandom(unsigned seed)
{
	for (int i = 0; i < RGENS; i++)
		rg_table[i].key = 0xffffffff;
	const int MAXWARM = 1223;
	seed ^= 0xbf14ef80; // just in case the user passes '0'...
	// initialize and warm-up the zeroth random generator:
	rg_table[0].r.seed(seed);
	for (int i = 0; i < MAXWARM; i++) rg_table[0].r._next();
	for (int i = 1; i < RGENS; i++) {
		Random& prev = rg_table[i - 1].r;
		Random& next = rg_table[i].r;
		next.seed(prev._next());
		int n = prev.randint(0, MAXWARM - 1);
		for (int i = 0; i < n; i++)
			next._next();
	}
}

Random& getRandomGen(int idx)
{
	unsigned key = idx;
	int i = ((unsigned) idx % (unsigned) RGENS);
	for (int k = 0; k < RGENS; k++) {
		if (rg_table[i].key == key)
			return rg_table[i].r;
		else if (rg_table[i].key == 0xffffffff) {
			rg_table[i].key = key;
			return rg_table[i].r;
		} else {
			i++;
			if (i >= RGENS) i -= RGENS;
		}
	}
	return rg_table[i].r;
}

Random& getRandomGen()
{
	return getRandomGen(SDL_ThreadID());
}

// random generator testing code below (disabled)

#if 0

#include <time.h>
#include <stdio.h>
#include <SDL/SDL.h>
#include "constants.h"
#include "color.h"
#include "sdl.h"

extern Color vfb[VFB_MAX_SIZE][VFB_MAX_SIZE];

static Random* grand;

double genrandpair(void)
{
	double x, y;
	grand->unitDiscSample(x, y);
	return x;
}

double threadid(void)
{
	return SDL_ThreadID();
}

void testSpeed(const char* what, double (*fn)  (void))
{
	unsigned cntPerSec = 1;
	Uint32 c0, c1;
	int remaining = 5;
	bool useful = false;
	while (remaining) {
		double sum = 0;
		c0 = SDL_GetTicks();
		for (unsigned i = 0; i < cntPerSec; i++)
			sum += fn();
		c1 = SDL_GetTicks();
		double t = (c1 - c0)/1000.0;
		sum /= cntPerSec;
		if (useful) printf("%.2lf Million %s per second (avg = %.5lf)\n",
			cntPerSec/t/1000000.0, what, sum);
		if (t > 1) {
			remaining--;
			useful = true;
		} else {
			cntPerSec *= 2;
		}
		
	}
}

static int int_buff[512][512];
static float float_buff[512][512];
static float circle_buff[512][512];
static float norm_buff[512];

static void displayCharts(double npoints)
{
	
	for (int y = 0; y < 511; y++)
		for (int x = 0; x < 511; x++) {
			float f = int_buff[y][x] * 0.2f;
			vfb[y][x] = Color(f, f, f);
	}
	for (int y = 0; y < 511; y++)
		for (int x = 0; x < 512; x++) {
			float f = float_buff[y][x] * 0.2f;
			vfb[y][x+512] = Color(f, f, f);
		}
	for (int y = 0; y < 512; y++)
		for (int x = 0; x < 511; x++) {
			float f = circle_buff[y][x] * 0.2f;
			vfb[y+512][x] = Color(f, f, f);
		}
	const int BORDERS = 16;
	const int NLINES = 8;
	const int SPACING = (512 - 2 * BORDERS) / NLINES;
	double mult = (512 - 2 * BORDERS) * 512 * 0.4 / (npoints);
	for (int x = 0; x < 512; x++) {
		Color base = Color(0, 0, 0);
		if ((int) fabs(x - 256) % ((int) (256.0/3.0)) == 0)
			base = Color(0.2f, 0.2f, 0.2f);
		int ey = nearestInt(512 - BORDERS - norm_buff[x] * mult);
		int sy = 512 - BORDERS;
		for (int y = 0; y < 512; y++) {
			Color f = base;
			if (sy - y >= 0 && (sy - y) % SPACING == 0) f += Color(0.1f, 0.1f, 0.1f);
			if (y == ey) f = Color(0.4f, 0.4f, 1.0f);
			else if (y == sy) f = Color(0.9f, 0.9f, 0.9f);
			else if (y > sy) f.makeZero();
			else if (y > ey) f += Color(0.2f, 0.2f, 0.5f);
			vfb[y+512][x+512] = f;
		}
	}
	for (int i = 0; i < 1024; i++)
		vfb[i][511] = vfb[511][i] = Color(1, 1, 1);
}

void test_random()
{
	initGraphics(1024, 1024);
	Random rnd(time(NULL));
	//
	for (int i = 0; i < 200; i++)
		printf("mt: %u\n", rnd._next());
	
	SDL_WM_SetCaption("Testing MTRandom", NULL);
	
	grand = &rnd;
	
	testSpeed("unitDiscSample()s", genrandpair);
	
	testSpeed("SDL_ThreadID()s", threadid);
	
	double delayFactor = 1000.0;
	// test the random generator graphically:
	for (int i = 0; i < 30000; i++) {
		// generate new 1000 integer points:
		for (int j = 0; j < 1000; j++)
			int_buff[rnd.randint(0, 511)][rnd.randint(0, 511)]++;
		
		// generate new 1000 floatingpoint points:
		for (int j = 0; j < 1000; j++) {
			float x = rnd.randfloat()*512;
			float y = rnd.randfloat()*512;
			int x0 = (int) floor(x);
			int y0 = (int) floor(y);
			int x1 = (x0 + 1) % 512;
			int y1 = (y0 + 1) % 512;
			float p = x - x0;
			float q = y - y0;
			float_buff[y0][x0] += (1 - p) * (1 - q);
			float_buff[y0][x1] += (    p) * (1 - q);
			float_buff[y1][x0] += (1 - p) * (    q);
			float_buff[y1][x1] += (    p) * (    q);
		}
		
		// generate new 1000 circle points:
		for (int j = 0; j < 1000; j++) {
			double cx, cy;
			rnd.unitDiscSample(cx, cy);
			float x = float(cx * 256 + 256);
			float y = float(cy * 256 + 256);
			int x0 = (int) floor(x);
			int y0 = (int) floor(y);
			int x1 = (x0 + 1) % 512;
			int y1 = (y0 + 1) % 512;
			float p = x - x0;
			float q = y - y0;
			circle_buff[y0][x0] += (1 - p) * (1 - q);
			circle_buff[y0][x1] += (    p) * (1 - q);
			circle_buff[y1][x0] += (1 - p) * (    q);
			circle_buff[y1][x1] += (    p) * (    q);
		}
		
		// generate new 5000 normally-distributed points:
		for (int j = 0; j < 5000; j++) {
			double x = rnd.gaussian(256, 256.0/3.0);
			if (x < 0 || x >= 511) continue;
			double p = x - floor(x);
			int x0 = floor(x);
			norm_buff[x0] += (1 - p);
			norm_buff[x0 + 1] += p;
		}
		//
		displayCharts((i + 1) * 5000);
		displayVFB(vfb);
		SDL_Delay((int) delayFactor);
		delayFactor *= 0.9;
	}
	waitForUserExit();
	closeGraphics();
}
#endif
