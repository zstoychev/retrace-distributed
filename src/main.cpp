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

#include "scene.h"
#include "random_generator.h"
#include "constants.h"
#include <ctime>
#include <cstring>

const char* host = "localhost";
int port = DEFAULT_PORT;
const char* dataDir = "data";

void printUsage() {
	printf("Invalid command arguments!\n");
}

extern int startServer(int argc, char** argv);
extern int startRenderer(int argc, char** argv);

int main(int argc, char** argv)
{
	initRandom((Uint32) time(NULL));

	if (argc < 2) {
		printUsage();
		return 1;
	}
	if (strcmp(argv[1], "-s") == 0 || strcmp(argv[1], "-ss") == 0) {
		return startServer(argc, argv);
	} else if (strcmp(argv[1], "-c") == 0) {
		return startRenderer(argc, argv);
	} else {
		printUsage();
		return 1;
	}

	freeCachedItems();

	return 0;
}
