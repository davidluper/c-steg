//
//  CStegApi.cpp
//  c-steg
//
//  this file contains what is meant to be the external api for
//  the c-steg library
//  

#include "CStegApi.h"
#include "c-steg/Utils.h"
#include "c-steg/LossyDataFormatter.h"

/**
 * this function exposes the ability to encode a lossy byte payload into a jpeg image.  the
 * lossy payload is meant to survive recompressions
 *
 * @param inputFileName
 *     the full path/name to the jpeg to encode information into
 * @param outputFileName
 *     the full path/name to the output image that contains the encoded information
 * @param data
 *     a pointer to the data to encode into the image
 * @param dataLength
 *     the length of the data in the data buffer of information to encode
 * @param quality
 *     the quality number to use when saving the image with the encoded information
 * @param numberOfCompressions
 *     the number of times to uncompress and recompress the image to stabilize the DCT coefficients
 *
 * returns true if successful false otherwise
 */
bool CSteg::encodeLossyJpeg(const char* inputFileName, const char* outputFileName, const char* data,
                              unsigned int dataLength, int quality, int numberOfCompressions) {
    StegUtils::doCompressions(numberOfCompressions, inputFileName, quality);
    StegUtils::purge(inputFileName, numberOfCompressions - 1);

    char file_name_buffer[256];
    memset(file_name_buffer, 0, 256);
    Utils::getFileNameWithoutExtension(inputFileName, file_name_buffer);

    char file_path_buffer[256];
    memset(file_path_buffer, 0, 256);
    int i = numberOfCompressions - 1;
    sprintf(file_path_buffer, "%s-%d.jpeg", file_name_buffer, i);

    struct packet_collection* to_packet_data = LossyDataFormatter::toLossy(data, dataLength);

    bool rtn = StegUtils::encodeLossyJpegFile(file_path_buffer, outputFileName, to_packet_data, quality);

    remove(file_path_buffer);
    free(to_packet_data->data);
    free(to_packet_data);

    return rtn;
}

/*
 * this function exposes the ability to decode a lossy byte payload from a jpeg that was encoded
 * with lossy formatted information
 *
 * @param inputFileName
 *     the full path/name to the jpeg image to decode information from
 * @param dataOut
 *     a pointer that if successful will point to the decoded information
 * @param dataOutLength
 *     a pointer to an unsigned int that if successful will be set to the 
 *     length of the decoded information in data_out
 *
 * returns true if successful false otherwise
 */
bool CSteg::decodeLossyJpeg(const char* inputFileName, DecodeResponse& decodeResponse) {
    struct packet_collection from_packet_data;
    bool rtn = StegUtils::decodeLossyJpegFile(inputFileName, &from_packet_data);

    if (rtn) {
        LossyDataFormatter::fromLossy(&from_packet_data, decodeResponse);
    }
    
    free(from_packet_data.data);
    
    return rtn;
}

/**
 * this function exposes the ability to encode lossless information in a jpeg.  its assumed 
 * when decoding the image it has not been recompressed after the encoding the steg payload.
 *
 * @param inputFileName
 *     the full path/name to the jpeg to encode information into
 * @param outputFileName
 *     the full path/name to the output image that contains the encoded information
 * @param data
 *     a pointer to the data to encode into the image
 * @param dataLength
 *     the length of the data in the data buffer of information to encode
 * @param numberOfBitsToSteal
 *     the number of bits to steal per coefficient when embedding information
 * @param quality
 *     the quality number to use when saving the image with the encoded information
 * 
 * returns true if successful false otherwise
 */
bool CSteg::encodeLosslessJpeg(const char *inputFileName, const char* outputFileName,
                               unsigned char* data, int dataLength,
                               unsigned short numberOfBitsToSteal, int quality) {
    return StegUtils::encodeLosslessJpegFile(inputFileName, outputFileName, data, dataLength, numberOfBitsToSteal, quality);
}

/**
 * this function exposes the ability to decode information that was encoded with the lossless 
 * format.  it is assumed the image has not been recompressed.
 *
 * @param inputFileName
 *     the full path/name to the jpeg image to decode information from
 * @param decodeResponse
 * 	   an object of type DecodeResponse with which to store the decoded
 * 	   output data
 *
 * returns true if successful false otherwise
 */
bool CSteg::decodeLosslessJpeg(const char *inputFileName, DecodeResponse& decodeResponse) {
    return StegUtils::decodeLosslessJpegFile(inputFileName, decodeResponse);
}

/*
 * this function calculates the number of usable coefficients for encoding
 * information into an image.
 *
 * @param inputFileName
 *     the full name/path to the jpeg to calculate the number of usable coefficient for
 *
 * returns the number of usable coefficients for embedding information
 */
unsigned int CSteg::getNumberOfCoefficients(const char *inputFileName) {
    return StegUtils::getNumberOfCoefficients(inputFileName);
}

/**
 * this function does a coefficient diff between to jpeg image files on disk
 *
 * @param filePath1
 *     the full path/name to the first image on disk
 * @param filePath2
 *     the full path/name to the second image on disk
 * returns the number of coefficients that differ between the two sets of coefficients
 */
int CSteg::doCoefficientDiff(const char* filePath1, const char* filePath2) {
	return StegUtils::doCoefficientDiff(filePath1, filePath2);
}
