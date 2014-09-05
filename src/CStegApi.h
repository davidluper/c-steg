//  c_steg_api.h
//  c-steg
//
//  this header file defines what is meant to be the externally visible api
//  for the c-steg library
//

#ifndef C_STEG_API_H
#define C_STEG_API_H

#include "c-steg/StegUtils.h"

#ifdef __cplusplus
extern "C" {
#endif

class CSteg {
public:
	static bool encodeLossyJpeg(const char* inputFile, const char* outputFile, const char* data,
                              unsigned int dataLength, int quality, int numberOfCompressions);
	static bool decodeLossyJpeg(const char* inputFile, DecodeResponse& decodeResponse);
	static bool encodeLosslessJpeg(const char *inputFileName, const char* outputFileName,
                               unsigned char* data, int dataLength,
                               unsigned short numberOfBitsToSteal, int quality);
	static bool decodeLosslessJpeg(const char *inputFileName, DecodeResponse& decodeResponse);
	static unsigned int getNumberOfCoefficients(const char* inputFileName);
	static int doCoefficientDiff(const char* filePath1, const char* filePath2);
};

#ifdef __cplusplus
}
#endif

#endif
