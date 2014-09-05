OUTPUT_FILE=../c-steg
CXX=g++
CXXFILES= \
main.cpp \
c-steg/JpegUtils.cpp \
c-steg/LossyDataFormatter.cpp \
c-steg/Utils.cpp \
c-steg/StegUtils.cpp \
CStegApi.cpp

CXXFLAGS=-O3 -o $(OUTPUT_FILE) -rdynamic -D_GNU_SOURCE
STATIC_LIBS=./libjpeg/libjpeg.a

all:
	$(CXX) $(CXXFILES) $(STATIC_LIBS) $(CXXFLAGS) -g

release:
	$(CXX) $(CXXFILES) $(STATIC_LIBS) $(CXXFLAGS)
	
clean:
	rm -f $(OUTPUT_FILE) *.o