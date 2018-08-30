CC=g++
CFLAGS=-c -g -Wall -DHAVE_INLINE -O2 -fexceptions -fopenmp  -std=c++0x
IFLAGS=-I include/ -I /usr/include -I Sky/Novas/  
LDFLAGS=-L /usr/lib -l netcdf_c++ -l netcdf -lgsl -lgslcblas -l fftw3  -lcxsparse -lm  -fopenmp 
LDFLAGS_FITS= -lcfitsio -l CCfits




#Common source files
CSOURCES=  Clean/AzElTemplateCalculator.cpp  Clean/Clean.cpp Clean/CleanPCA.cpp Clean/CleanBspline.cpp Clean/CleanHigh.cpp Clean/CleanIterativeMean.cpp   Clean/CleanSelector.cpp Observatory/Detector.cpp Utilities/tinyxml2.cpp Utilities/vector_utilities.cpp Observatory/Array.cpp Observatory/Telescope.cpp Observatory/TimePlace.cpp Sky/Source.cpp Utilities/GslRandom.cpp Analysis/AnalParams.cpp Sky/astron_utilities.cpp Mapmaking/Map.cpp Mapmaking/Observation.cpp Mapmaking/Coaddition.cpp Mapmaking/PointSource.cpp Mapmaking/CompletenessSim.cpp Mapmaking/NoiseRealizations.cpp Mapmaking/WienerFilter.cpp Utilities/gaussFit.cpp Utilities/BinomialStats.cpp Sky/Novas/novas.c Sky/Novas/novascon.c Sky/Novas/nutation.c Sky/Novas/solsys1.c Sky/Novas/eph_manager.c Sky/Novas/readeph0.c Utilities/SBSM.cpp  Utilities/convolution.cpp Utilities/mpfit.cpp Utilities/sparseUtilities.cpp  



#Macana executable definitions
SOURCES=macanap.cpp $(CSOURCES)
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=macana
#fitswriter executable definitions
SOURCES_FITS=Utilities/fitswriter.cpp 
OBJECTS_FITS=$(SOURCES_FITS:.cpp=.o)
EXECUTABLE_FITS=fitswriter

all: $(EXECUTABLE)  $(EXECUTABLE_FITS) 

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)  
	ln -s -f -T macana macanap
	mkdir -p bin 
	cp macana bin/
	cp python/aztec/bin/lmt_reduce bin/
	cp python/aztec/bin/macana_cdirs bin/
		
$(EXECUTABLE_FITS): $(OBJECTS_FITS)
	$(CC) $(OBJECTS_FITS) -o $@ $(LDFLAGS) $(LDFLAGS_FITS)
	cp fitswriter bin/

%.o: %.cpp
	$(CC) $(CFLAGS) $(IFLAGS) $< -o $@

.PHONY: clean

clean:
	rm -rf *.o *~ core */*.o */*~ *.op */*.op $(EXECUTABLE) $(EXECUTABLE_FITS) bin/* */*.pyc
