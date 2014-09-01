CXX=clang
CXXFILES=main.cpp jpeg_tools.cpp lossy_data_formatter.cpp c_steg_api.cpp
CXXFLAGS=-O3 -o c-steg -rdynamic -D_GNU_SOURCE -g
STATIC_LIBS=./jpeg/libjpeg.a

all:
	$(CXX) $(CXXFILES) $(STATIC_LIBS) $(CXXFLAGS)

clean:
	rm -f c-steg *.o