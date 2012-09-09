bin_PROGRAMS = retrace
retrace_SOURCES = bitmap.cpp \
	camera.cpp \
	cxxptl_sdl.cpp \
	environment.cpp \
	geometry.cpp \
	heightfield.cpp \
  interactivehandler.cpp \
	lights.cpp \
	main.cpp \
	matrix.cpp \
	mesh.cpp \
	random_generator.cpp \
  renderer.cpp \
  renderserver.cpp \
	scene.cpp \
  tcpsocket.cpp \
	sdl.cpp \
	shading.cpp

# set the include path found by configure
AM_CPPFLAGS =  $(LIBSDL_CFLAGS) $(all_includes)

# the library search path.
retrace_LDFLAGS = $(all_libraries)
retrace_LDADD = $(LIBSDL_LIBS)
noinst_HEADERS = bbox.h \
	bitmap.h \
	camera.h \
	color.h \
	constants.h \
	cxxptl_sdl.h \
	environment.h \
	geometry.h \
	heightfield.h \
	lights.h \
	matrix.h \
	mesh.h \
	random_generator.h \
  renderer.h \
  renderserver.h \
	scene.h \
	sdl.h \
	shading.h \
  tcpsocket.h \
	util.h \
	vector.h
