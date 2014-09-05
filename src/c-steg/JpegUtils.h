//
//  JpegUtils.h
//  c-steg
//

#ifndef CSTEG_JPEG_TOOLS_H
#define CSTEG_JPEG_TOOLS_H

#include "Utils.h"

class JpegUtils {
public:
	static int readJpegFile(Image& image, const char* inputFileName);
	static int writeJpegFile(Image& image, const char *outputFileName, int quality);
	static bool writeJpegCoefficients(JpegCoefficients& coefficients, const char* outputFileName, int quality);
	static bool readJpegCoefficients(const char *inputFileName, JpegCoefficients& coefficients);
};

#endif
