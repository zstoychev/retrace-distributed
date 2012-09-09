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

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "color.h"
#include "constants.h"
#include "bitmap.h"
#include <algorithm>
using std::swap;

Bitmap::Bitmap()
{
	width = height = -1;
	data = NULL;
	differentiated = false;
}

Bitmap::~Bitmap()
{
	freeMem();
}

void Bitmap::freeMem(void)
{
	if (data) delete [] data;
	data = NULL;
	width = height = -1;
}

int Bitmap::getWidth(void) const { return width; }
int Bitmap::getHeight(void) const { return height; }
bool Bitmap::isOK(void) const { return (data != NULL); }

void Bitmap::generateEmptyImage(int w, int h)
{
	freeMem();
	if (w <= 0 || h <= 0) return;
	width = w;
	height = h;
	data = new Color[w * h];
	memset(data, 0, sizeof(data[0]) * w * h);
}

bool Bitmap::loadImage(const char* filename)
{
	char extension[8];
	int l = (int) strlen(filename);
	if (l > 3) {
		strcpy(extension, &filename[l - 4]);
		for (int i = 0; i < 4; i++) extension[i] = toupper(extension[i]);
		if (!strcmp(extension, ".PFM")) return loadPFM(filename);
	}
	return loadBMP(filename);
}

Color Bitmap::getPixel(int x, int y) const
{
	if (!data || x < 0 || x >= width || y < 0 || y >= height) return Color(0.0f, 0.0f, 0.0f);
	return data[x + y * width];
}

Color Bitmap::getFilteredPixel(float x, float y) const
{
 	if (!data || x < 0 || x >= width || y < 0 || y >= height) return Color(0.0f, 0.0f, 0.0f);
	int x0 = (int) floor(x);
	int y0 = (int) floor(y);
	int x1 = (x0 + 1) % width;
	int y1 = (y0 + 1) % height;
	float p = x - x0;
	float q = y - y0;
    return data[x0 + y0 * width] * (1 - p) * (1 - q)
         + data[x1 + y0 * width] * (    p) * (1 - q)
         + data[x0 + y1 * width] * (1 - p) * (    q)
         + data[x1 + y1 * width] * (    p) * (    q);
}

void Bitmap::setPixel(int x, int y, const Color& color)
{
	if (!data || x < 0 || x >= width || y < 0 || y >= height) return;
	data[x + y * width] = color;
}

class ImageOpenRAII {
	Bitmap *bmp;
public:
	bool imageIsOk;
	FILE* fp;
	ImageOpenRAII(Bitmap *bitmap)
	{
		fp = NULL;
		bmp = bitmap;
		imageIsOk = false;
	}
	~ImageOpenRAII()
	{
		if (!imageIsOk) bmp->freeMem();
		if (fp) fclose(fp); fp = NULL;
	}
	
};

const int BM_MAGIC = 19778;

struct BmpHeader {
	int fs;       // filesize
	int lzero;
	int bfImgOffset;  // basic header size
};
struct BmpInfoHeader {
	int ihdrsize; 	// info header size
	int x,y;      	// image dimensions
	unsigned short channels;// number of planes
	unsigned short bitsperpixel;
	int compression; // 0 = no compression
	int biSizeImage; // used for compression; otherwise 0
	int pixPerMeterX, pixPerMeterY; // dots per meter
	int colors;	 // number of used colors. If all specified by the bitsize are used, then it should be 0
	int colorsImportant; // number of "important" colors (wtf?..)
};

bool Bitmap::loadBMP(const char* filename)
{
	freeMem();
	ImageOpenRAII helper(this);
	
	BmpHeader hd;
	BmpInfoHeader hi;
	Color palette[256];
	int toread = 0;
	unsigned char *xx;
	int rowsz;
	unsigned short sign;
	FILE* fp = fopen(filename, "rb");
	
	if (fp == NULL) return false;
	helper.fp = fp;
	if (!fread(&sign, 2, 1, fp)) return false;
	if (sign != BM_MAGIC) {
		printf("loadBMP: `%s' is not a BMP file.\n", filename);
		return false;
	}
	if (!fread(&hd, sizeof(hd), 1, fp)) return false;
	if (!fread(&hi, sizeof(hi), 1, fp)) return false;
	
	/* header correctness checks */
	if (!(hi.bitsperpixel == 8 || hi.bitsperpixel == 24 ||  hi.bitsperpixel == 32)) {
		printf("loadBMP: Cannot handle file format at %d bpp.\n", hi.bitsperpixel); 
		return false;
	}
	if (hi.channels != 1) {
		printf("loadBMP: cannot load multichannel .bmp!\n");
		return false;
	}
	/* ****** header is OK *******/
	
	// if image is 8 bits per pixel or less (indexed mode), read some pallete data
	if (hi.bitsperpixel <= 8) {
		toread = (1 << hi.bitsperpixel);
		if (hi.colors) toread = hi.colors;
		for (int i = 0; i < toread; i++) {
			unsigned temp;
			if (!fread(&temp, 1, 4, fp)) return false;
			palette[i] = Color(temp);
		}
	}
	toread = hd.bfImgOffset - (54 + toread*4);
	fseek(fp, toread, SEEK_CUR); // skip the rest of the header
	int k = hi.bitsperpixel / 8;
	rowsz = hi.x * k;
	if (rowsz % 4 != 0)
		rowsz = (rowsz / 4 + 1) * 4; // round the row size to the next exact multiple of 4
	xx = new unsigned char[rowsz];
	generateEmptyImage(hi.x, hi.y);
	if (!isOK()) {
		printf("loadBMP: cannot allocate memory for bitmap! Check file integrity!\n");
		delete [] xx;
		return false;
	}
	for (int j = hi.y - 1; j >= 0; j--) {// bitmaps are saved in inverted y
		if (!fread(xx, 1, rowsz, fp)) {
			printf("loadBMP: short read while opening `%s', file is probably incomplete!\n", filename);
			delete [] xx;
			return 0;
		}
		for (int i = 0; i < hi.x; i++){ // actually read the pixels
			if (hi.bitsperpixel > 8)
				setPixel(i, j, Color(xx[i*k+2]/255.0f, xx[i*k+1]/255.0f, xx[i*k]/255.0f));
			else
				setPixel(i, j,  palette[xx[i*k]]);
		}
	}
	delete [] xx;
	
	helper.imageIsOk = true;
	return true;
}

bool Bitmap::saveBMP(const char* filename, float outputGamma)
{
	FILE* fp = fopen(filename, "wb");
	if (!fp) return false;
	BmpHeader hd;
	BmpInfoHeader hi;
	char xx[VFB_MAX_SIZE * 3];


	// fill in the header:
	int rowsz = width * 3;
	if (rowsz % 4)
		rowsz += 4 - (rowsz % 4); // each row in of the image should be filled with zeroes to the next multiple-of-four boundary
	hd.fs = rowsz * height + 54; //std image size
	hd.lzero = 0;
	hd.bfImgOffset = 54;
	hi.ihdrsize = 40;
	hi.x = width; hi.y = height;
	hi.channels = 1;
	hi.bitsperpixel = 24; //RGB format
	// set the rest of the header to default values:
	hi.compression = hi.biSizeImage = 0;
	hi.pixPerMeterX = hi.pixPerMeterY = 0;
	hi.colors = hi.colorsImportant = 0;
	
	fwrite(&BM_MAGIC, 2, 1, fp); // write 'BM'
	fwrite(&hd, sizeof(hd), 1, fp); // write file header
	fwrite(&hi, sizeof(hi), 1, fp); // write image header
	for (int y = height - 1; y >= 0; y--) {
		for (int x = 0; x < width; x++) {
			unsigned t = getPixel(x, y).toRGB32(16, 8, 0, outputGamma);
			xx[x * 3    ] = (0xff     & t);
			xx[x * 3 + 1] = (0xff00   & t) >> 8;
			xx[x * 3 + 2] = (0xff0000 & t) >> 16;
		}
		fwrite(xx, rowsz, 1, fp);
	}
	fclose(fp);
	return true;
}

static bool scanToNextWS(char data[], char res[], int& curIdx, int maxLen)
{
	int j = 0;
	int i = curIdx;
	while (i < maxLen && !isspace(data[i])) {
		res[j++] = data[i++];
	}
	if (i >= maxLen) return false;
	curIdx = i + 1;
	res[j] = 0;
	return true;
}

bool Bitmap::tryLoadPFM(FILE* f)
{
	char data[128], s[128];
	/* 
	 * A description of the file format is available at http://netpbm.sourceforge.net/doc/pfm.html
	 * PFM file start with a text header like
	 *
	 * PF
	 * 1024 768
	 * -1.000
	 *
	 * which describes the image size, endianines (because of the negative third row). After these three
	 * lines, a raster of binary (4-byte float) data follows.
	 */
	if (sizeof(data) != fread(data, 1, sizeof(data), f)) return false;
	int i = 0;
	if (!scanToNextWS(data, s, i, sizeof(data))) return false;
	if (strlen(s) != 2 || s[0] != 'P' || toupper(s[1] != 'F')) return false;
	bool gray = s[1] == 'f';
	int w, h;
	int bigEndian;
	double scaling;
	if (!scanToNextWS(data, s, i, sizeof(data)) || 1 != sscanf(s, "%d", &w)) return false;
	if (!scanToNextWS(data, s, i, sizeof(data)) || 1 != sscanf(s, "%d", &h)) return false;
	if (!scanToNextWS(data, s, i, sizeof(data)) || 1 != sscanf(s, "%lf", &scaling)) return false;
	bigEndian = scaling > 0;
	fseek(f, i, SEEK_SET);
	int lineSize = w * (gray ? 1 : 3) * 4;
	char *line = new char[lineSize];
	float* fline = (float*) line;
	generateEmptyImage(w, h);
	for (int y = 0; y < h; y++) {
		if (lineSize != (int) fread(line, 1, lineSize, f)) {
			delete[] line;
			return false;
		}
		if (bigEndian) {
			for (int i = 0; i < lineSize; i += 4) {
				swap(line[i], line[i + 3]);
				swap(line[i + 1], line[i + 2]);
			}
		}
		if (gray) {
			for (int x = 0; x < w; x++)
				setPixel(x, y, Color(fline[x], fline[x], fline[x]));
		} else {
			for (int x = 0; x < w; x++)
				setPixel(x, y, Color(fline[x*3], fline[x*3+1], fline[x*3+2]));
		}
	}
	delete[] line;
	return true;
}


bool Bitmap::loadPFM(const char* filename)
{
	FILE* f = fopen(filename, "rb");
	if (!f) return false;
	bool loadResult = tryLoadPFM(f);
	fclose(f);
	return loadResult;
}

void Bitmap::differentiate(void)
{
	if (!differentiated) {
		int W = getWidth(), H = getHeight();
		Color* newData = new Color[W*H];
	
		for (int y = 0; y < H; y++) {
			for (int x = 0; x < W; x++) {
				Color us = getPixel(x, y);
				Color right = getPixel((x + 1) % W, y);
				Color down = getPixel(x, (y + 1) % H);
				float dx = us.intensity() - right.intensity();
				float dy = us.intensity() - down.intensity();
				newData[y * W + x] = Color(dx, dy, 0);
			}
		}
		delete[] data;
		data = newData;

		differentiated = true;
	}
}
