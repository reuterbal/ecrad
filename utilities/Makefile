SOURCES = easy_netcdf.F90 radiation_io.F90 radiation_random_numbers_mix.F90

OBJECTS := $(SOURCES:.F90=.o)
LIBUTILITIES = ../lib/libutilities.a

MAKE_DEPS = ../bin/make_deps.sh
DEPS_FILE = Makefile_deps

all: $(LIBUTILITIES)

$(LIBUTILITIES): $(OBJECTS)
	ar r $(LIBUTILITIES) $(OBJECTS)

%.o: %.F90
	$(FC) $(FCFLAGS) -c $<

$(DEPS_FILE): $(SOURCES)
	$(MAKE_DEPS) $(SOURCES) > $(DEPS_FILE)

clean:
	rm -f *.o $(LIBUTILITIES)

include Makefile_deps
