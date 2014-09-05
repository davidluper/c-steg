//
//  StegUtils.cpp
//  c-steg
//

#include "StegUtils.h"
#include "JpegUtils.h"
#include <time.h>
#include <iostream>
#include <math.h>

#define COEFFICIENT_MAX 32765 //give some room to spare (actual max should be 32767)
#define COEFFICIENT_MIN -32765 //give some room to spare (actual in should be -32767)
#define CSTEG_JPEG_VERSION 1
#define CSTEG_JPEG_VERSION_1_HEADER_LENGTH 51 //16 + 32 + 3
#define VERSION_BIT_LENGTH 16
#define DATA_LENGTH_BIT_LENGTH 32
#define BITS_TO_STEAL_BIT_LENGTH 3

/**
 * this function is used in the lossless encoding function to modify the coefficient
 * to hold a single bit from the data. this is used for encoding the header of the lossless
 * encoder
 *
 * @param coef
 *     the coefficient that need to hold the information from the data
 * @ param coef_bits
 *     the bits from the coefficient to compare against the data_bits to
 *     see if the coefficient needs to be modified
 * @param data_bits
 *     the bits from the data that need to be in the coefficient
 *
 * returns the coefficient modified (if necessary) to hold the information
 * from the data
 */
short StegUtils::getLosslessEncodedShort(short coef, short coefBits, short dataBits) {
    // this will always be a single bit, if the bit from data and coef already
    // equal each other then nothing is necessary
    if (coefBits != dataBits) {
        if (coef >= COEFFICIENT_MAX) {
            //be careful to never over flow
            --coef;
        } else if (coef <= COEFFICIENT_MIN) {
            ++coef;
        } else {
            //if the bits are different then modify the coefficient bit
            //to hold the information from the data
            switch (coef) {
            case 1:
                //never let a coefficient that was not already zero turn into a zero
                ++coef;
                break;
            case -1:
                //never let a coefficient that was not already zero turn into a zero
                --coef;
                break;
            default:
                //randomly either go up one or down one to modifiy the coefficient bit
                if (rand() % 2 == 1) {
                    ++coef;
                } else {
                    --coef;
                }
            }
        }
    }

    return coef;
}

/**
 * this function uses LSB encoding to embed a data payload into the DCT coefficients of
 * a jpeg image.  this functions assumes the decode lossless function will be run on
 * the exact image produced as ouput of this routine.  if the output image from this
 * routine is altered in any way decoding will fail.  if the image is expected to be
 * altered either by cropping or recompression you should use the lossy encode/decode
 * functions
 *
 * @param inputFilName
 *     the complete file path/name to the image to embed information into
 * @param outputFileName
 *     the complete file path/name to the output image that will contain
 *     the image with embedded information
 * @param data
 *     the data (presented as an unsigned char array) to embed into the coefficient
 *     of the specified jpeg
 * @param dataLength
 *     the number of unsigned char bytes to encode from data (the unsigned
 *     char buffer input)
 * @param numberOfBitsToSteal
 *     the number of least significant bits to steal from a coefficient and replace
 *     with data bits
 * @param quality
 *     the quality number to use when saving the jpeg with embeded information
 *
 * returns true if everything went as planned, false if an error occurred
 */
bool StegUtils::encodeLosslessJpegFile(const char *inputFileName, const char* outputFileName,
		const unsigned char* const data, int dataLength, unsigned short numberOfBitsToSteal, int quality) {
	JpegCoefficients inputJpeg;
	if (!JpegUtils::readJpegCoefficients(inputFileName, inputJpeg)) {
		std::cout << "error loading jpeg file: " << inputFileName << std::endl;
		return false;
	}

	// -------------------------------------------------------------------------------------------------
	// begin LSB methodology
	// loop over the coefficients and replace LSBs with data bits
	// -------------------------------------------------------------------------------------------------
	int numberOfCoefficients = 0;
	int dataByteIndex = 0;
	int dataBitIndex = 0;
	int coef = 0;
	int coefBits = 0;
	int dataBits = 0;
	int dataMask = (int) (pow((float) 2, numberOfBitsToSteal) - 1);
	int lengthPayloadBitCount = 0;
	int versionPayloadBitCount = 0;
	int numberOfBitsToStealPayloadCount = 0;
	int lengthPayload = dataLength;
	int dataLengthMinusOne = dataLength - 1;
	unsigned short shortBuf;
	bool headerComplete = false;
	int encodeBitCount = 0;
	bool stop = false;
	time_t t;
	srand((unsigned) time(&t));
	for (JDIMENSION compnum = 0; compnum < inputJpeg.numComponents; ++compnum) {
		if (stop) {
			break;
		}

		for (JDIMENSION rownum = 0; rownum < inputJpeg.heightInBlocks[compnum]; ++rownum) {
			if (stop) {
				break;
			}

			for (JDIMENSION blocknum = 0; blocknum < inputJpeg.widthInBlocks[compnum]; ++blocknum) {
				if (stop) {
					break;
				}

				for (JDIMENSION i = 0; i < DCTSIZE2; ++i) {
					coef = inputJpeg.coefBuffers[compnum][rownum][blocknum][i];
					//only steal coefficients that don't have a zero value
					if (coef != 0) {
						if (headerComplete) {
							//if we have not encoded the entire payload
							if (dataByteIndex < dataLength) {
								// use a short buf here to pack the current
								// byte and the next byte of the data.  this
								// makes it easier for values for bits to steal like
								// 3 which can spill over between bytes
								shortBuf = data[dataByteIndex];
								if (dataByteIndex < dataLengthMinusOne) {
									shortBuf = (shortBuf << 8) | data[dataByteIndex + 1];
								} else {
									shortBuf = shortBuf << 8;
								}

								dataBits = (shortBuf >> (16 - (dataBitIndex + numberOfBitsToSteal))) & dataMask;
								coef = (coef >> numberOfBitsToSteal << numberOfBitsToSteal) | dataBits;
								if (coef == 0) {
									coef = coef | (1 << numberOfBitsToSteal);
								}

								//add the bits we stole to the bit count and byte count
								dataBitIndex += numberOfBitsToSteal;
								if (dataBitIndex >= 8) {
									dataBitIndex %= 8;
									++dataByteIndex;
								}

								inputJpeg.coefBuffers[compnum][rownum][blocknum][i] = coef;
							} else {
								stop = true;
								break;
							}
						} else {
							// the header contains the version number, payload length,
							// and bits to steal encoded by stealing 1 bit the payload
							// will use the number of bits to steal passed into the function
							if (versionPayloadBitCount < VERSION_BIT_LENGTH) {
								coefBits = coef & 1;
								dataBits = (CSTEG_JPEG_VERSION >> versionPayloadBitCount) & 1;

								inputJpeg.coefBuffers[compnum][rownum][blocknum][i] =
										StegUtils::getLosslessEncodedShort(coef, coefBits, dataBits);
								++versionPayloadBitCount;
							} else if (lengthPayloadBitCount < DATA_LENGTH_BIT_LENGTH) {
								coefBits = coef & 1;
								dataBits = (lengthPayload >> lengthPayloadBitCount) & 1;

								inputJpeg.coefBuffers[compnum][rownum][blocknum][i] =
										StegUtils::getLosslessEncodedShort(coef, coefBits, dataBits);
								++lengthPayloadBitCount;
							} else if (numberOfBitsToStealPayloadCount < BITS_TO_STEAL_BIT_LENGTH) {
								//encode the number of bits we stole when encoding this image
								coefBits = coef & 1;
								dataBits = (numberOfBitsToSteal >> numberOfBitsToStealPayloadCount) & 1;

								inputJpeg.coefBuffers[compnum][rownum][blocknum][i] =
										StegUtils::getLosslessEncodedShort(coef, coefBits, dataBits);
								++numberOfBitsToStealPayloadCount;
							} else {
								headerComplete = true;
							}

							++encodeBitCount;
						}

						++numberOfCoefficients;
					}
				}
			}
		}
	}

	return JpegUtils::writeJpegCoefficients(inputJpeg, outputFileName, quality);
}

/**
 * this function decodes the LSB payload from a jpeg image that was embedded through the encode lossless
 * function.  it is assumed the image being decoded is exactly the same (byte for byte) as the output
 * image from the encode lossless function.  if the encoded image was cropped, recompressed, or manipulated
 * in any way, the lossy encode/decode functions should be used.
 *
 * @param inputFileName
 *     the full file path/name to the image to be decoded
 * @param dataOut
 *     a double unsigned char pointer.  this function will figure out what length array the supplied pointer should
 *     point to an reassign the point to point to it.  THIS FUNCTION WILL ALLOCATE MEMORY THAT IS POINTED TO BY
 *     THIS.  the memory allocated by this function should be freed elsewhere!!!!!! responsibly deallocate please!!!
 * @param dataOutLength
 *     this unsigned int will be assigned a value equal to the length of the decoded data_out buffer
 *
 * returns true if everything went as expected, false otherwise.
 */
bool StegUtils::decodeLosslessJpegFile(const char *inputFileName, DecodeResponse& decodeResponse) {
	JpegCoefficients inputJpeg;
	if (!JpegUtils::readJpegCoefficients(inputFileName, inputJpeg)) {
		std::cout << "error loading jpeg file: " << inputFileName << std::endl;
		return false;
	}

	// -------------------------------------------------------------------------------------------
	// begin LSB decoding
	// -------------------------------------------------------------------------------------------
	int dataByteIndex = 0;
	int coef;
	int bufferLength = 0;
	unsigned char* buffer = NULL;
	unsigned int decodeBitCount = 0;
	unsigned char coefBit = 0;
	int coefficient_mask = 0;
	int intBuff = 1;
	int intBuffCount = 0;
	int dataBits = 0;
	bool stop = false;
	int numberOfBitsToSteal = 0;
	bool headerComplete = false;
	int versionPayloadBitCount = 0;
	int numberOfBitsToStealPayloadCount = 0;
	int lengthPayloadBitCount = 0;
	int version = 0;
	for (JDIMENSION compnum = 0; compnum < inputJpeg.numComponents; ++compnum) {
		if (stop) {
			break;
		}

		for (JDIMENSION rownum = 0; rownum < inputJpeg.heightInBlocks[compnum]; ++rownum) {
			if (stop) {
				break;
			}

			for (JDIMENSION blocknum = 0; blocknum < inputJpeg.widthInBlocks[compnum]; ++blocknum) {
				if (stop) {
					break;
				}

				for (JDIMENSION i = 0; i < DCTSIZE2; ++i) {
					coef = inputJpeg.coefBuffers[compnum][rownum][blocknum][i];
					if (coef != 0) {
						if (headerComplete) {
							if (dataByteIndex < bufferLength) {
								dataBits = coef & coefficient_mask;
								intBuff = (intBuff << numberOfBitsToSteal) | dataBits;
								intBuffCount += numberOfBitsToSteal;
								if (intBuffCount >= 8) {
									intBuffCount -= 8;
									buffer[dataByteIndex] = ((intBuff >> intBuffCount) & 255);
									++dataByteIndex;
								}
							} else {
								stop = true;
								break;
							}
						} else {
							if (versionPayloadBitCount < VERSION_BIT_LENGTH) {
								coefBit = coef & 1;

								version = version | (coefBit << versionPayloadBitCount);
								++versionPayloadBitCount;
							} else if (lengthPayloadBitCount < DATA_LENGTH_BIT_LENGTH) {
								coefBit = coef & 1;

								bufferLength = bufferLength | (coefBit << lengthPayloadBitCount);
								++lengthPayloadBitCount;
							} else if (numberOfBitsToStealPayloadCount < BITS_TO_STEAL_BIT_LENGTH) {
								//encode the number of bits we stole when encoding this image
								coefBit = coef & 1;

								numberOfBitsToSteal = numberOfBitsToSteal | (coefBit << numberOfBitsToStealPayloadCount);
								++numberOfBitsToStealPayloadCount;
							} else {
								headerComplete = true;
								buffer = (unsigned char*) malloc(bufferLength);
								memset(buffer, 0, bufferLength);
								coefficient_mask = (int) pow((float) 2, numberOfBitsToSteal) - 1;
							}

							++decodeBitCount;
						}
					}
				}
			}
		}
	}

	decodeResponse.data = buffer;
	decodeResponse.length = bufferLength;

	return true;
}

/*
 * this function encodes shorter amounts of information redundantly into a jpeg image.  this information
 * can potentially survive manipulation such as recompression or cropping.  this function takes a specially
 * formated input that is stored into a packet_collection.  this function puts packeted into the single LSB
 * of non-zero coefficients.  those packets have information beyond just the original data (checksum, rolling sum,
 * packet index, etc.) that allows the decode lossy function to attempt a decode.  the encode and decode methods here
 * are not responsible to formatting the packets of reassembling the original data.  that is done elsewhere
 * this method only encodes supplied data 1 bit at a time into "more reliable/stable" DCT coefficients.  there
 * is a layer above this that does the packet formatting and unpacking.
 *
 * @param inputFileName
 *     the full path/name to a jpeg image on disk to encode
 * @param outputFileName
 *     the full path/name to store the encoded jpeg image on disk
 * @param packetCollection
 *     the packet collection object that holds the formatted packets to encode into the single LSB of
 *     the non-zero coefficients of the input image
 * @param quality
 *     the quality to save the output image with (it would be advisable to make this equal to the
 *     quality of the input image)
 *
 * returns true if everything went as planned, false if an error occurred
 */
bool StegUtils::encodeLossyJpegFile(const char *inputFileName, const char* outputFileName,
		struct packet_collection* packetCollection, int quality) {
	JpegCoefficients inputJpeg;
	if (!JpegUtils::readJpegCoefficients(inputFileName, inputJpeg)) {
		std::cout << "error loading jpeg file: " << inputFileName << std::endl;
		return false;
	}

	int numberOfCoefficients = 0;
	int dataByteIndex = 0;
	int dataBitIndex = 0;
	int coef = 0;
	int coefBit = 0;
	int dataBit = 0;
	int x = 0;
	int y = 0;
	for (JDIMENSION compnum = 0; compnum < inputJpeg.numComponents; ++compnum) {
		for (JDIMENSION rownum = 0; rownum < inputJpeg.heightInBlocks[compnum]; ++rownum) {
			for (JDIMENSION blocknum = 0; blocknum < inputJpeg.widthInBlocks[compnum]; ++blocknum) {
				x = 0;
				y = 0;
				for (JDIMENSION i = 0; i < DCTSIZE2; ++i) {
					coef = inputJpeg.coefBuffers[compnum][rownum][blocknum][i];
					x = (x + 1) % DCTSIZE;
					if (x == 0) {
						++y;
					}

					if (coef != 0 && x > 1 && y > 1) {
						coefBit = coef & 1;
						dataBit = (packetCollection->data[dataByteIndex] >> dataBitIndex) & 1;

						++dataBitIndex;
						dataBitIndex %= 8;
						if (dataBitIndex == 0) {
							++dataByteIndex;
							dataByteIndex %= packetCollection->length;
						}

						if (dataBit != coefBit) {
							switch (coef) {
							case 1:
								++coef;
								break;
							case -1:
								--coef;
								break;
							default:
								if (coef >= COEFFICIENT_MAX) {
									--coef;
								} else if (coef <= COEFFICIENT_MIN) {
									++coef;
								} else {
									if ((coef ^ numberOfCoefficients) == 1) {
										++coef;
									} else {
										--coef;
									}
								}
							}
						}

						inputJpeg.coefBuffers[compnum][rownum][blocknum][i] = coef;
						++numberOfCoefficients;
					}
				}
			}
		}
	}

	return JpegUtils::writeJpegCoefficients(inputJpeg, outputFileName, quality);

	return true;
}

/**
 * this function attempts to decode information encoded into a jpeg image using the lossy data encoding function.
 * this function puts the extracted bit stream into a packet collection. it should be assumed by layers above
 * this one that the packet stream extracted from this image is corrupted.  bits could be flipped or dropped,
 * which means layers above this should not even assume the data is byte aligned!  this decode method is not
 * responsible for reassembling the original data.  that is done elsewhere this method only decodes supplied
 * data 1 bit at a time from non zero DCT coefficients.  there is a layer above this that does the information
 * reassembling from the packets
 *
 * @param input_file_name
 *     the full path/name to the input jpeg file to decode
 * @param packetCollection
 *     the packet collection object to store the extracted bit stream from the single LSB of non-zero DCT
 *     coefficients
 *
 * returns true if the information was extracted, false only if there was an error like file open error
 */
bool StegUtils::decodeLossyJpegFile(const char *inputFileName, struct packet_collection* packetCollection) {
	JpegCoefficients inputJpeg;
	if (!JpegUtils::readJpegCoefficients(inputFileName, inputJpeg)) {
		std::cout << "error loading jpeg file: " << inputFileName << std::endl;
		return false;
	}

	// ----------------------------------------------------------------------------
	// begin LSB decoding 1 bit only for lossy encoding/decoding
	// ----------------------------------------------------------------------------
	int number_of_coefficients = 0;
	int data_byte_index = 0;
	int data_bit_index = 0;
	int coef;
	int buffer_length = (sizeof(unsigned char) * inputJpeg.info->image_width
			* inputJpeg.info->image_height * inputJpeg.info->num_components) / 8;

	unsigned char* buffer = (unsigned char*) malloc(buffer_length);
	memset(buffer, 0, buffer_length);

	int x = 0;
	int y = 0;
	for (JDIMENSION compnum = 0; compnum < inputJpeg.numComponents; ++compnum) {
		for (JDIMENSION rownum = 0; rownum < inputJpeg.heightInBlocks[compnum]; ++rownum) {
			for (JDIMENSION blocknum = 0; blocknum < inputJpeg.widthInBlocks[compnum]; ++blocknum) {
				x = 0;
				y = 0;
				for (JDIMENSION i = 0; i < DCTSIZE2; ++i) {
					coef = inputJpeg.coefBuffers[compnum][rownum][blocknum][i];
					x = (x + 1) % DCTSIZE;
					if (x == 0) {
						++y;
					}

					if (coef != 0 && x > 1 && y > 1) {
						if ((coef & 1) == 1) {
							buffer[data_byte_index] = buffer[data_byte_index] | (1 << data_bit_index);
						}

						++data_bit_index;
						data_bit_index %= 8;
						if (data_bit_index == 0) {
							++data_byte_index;
						}

						++number_of_coefficients;
					}
				}
			}
		}
	}

	packetCollection->data = buffer;
	packetCollection->length = number_of_coefficients / 8;

	return true;
}

/**
 * this function takes an image and runs a series of compressions and then recompressions n times.
 * the series of images is left on disk, use purge to clean up the images.  uncompressing and then
 * recompressing an image makes the coefficients more stable for subsequent recompressions.  this
 * means that when anyone recompresses an image it will destroy less of our steg information.  the
 * images on disk will be named using the input image name.  so /some/path/my.jpeg will be turned
 * into a sequence of /some/path/my-0.jpeg, /some/path/my-1.jpeg, etc.  the final image in the
 * series will be n - 1.
 *
 * @param n
 *     the number of times to perform uncompression/recompression
 * @param filePath
 *     the full path/name to the input image on disk
 * @param quality
 *     the quality to use when compressing the image
 */
void StegUtils::doCompressions(int n, const char* filePath, int quality) {
	char file_path_buffer[256];
	Utils::getFileNameWithoutExtension(filePath, file_path_buffer);

	char file_name_buffer_1[256];
	char file_name_buffer_2[256];
	memset(file_name_buffer_1, 0, 256);

	char* ptr_1 = file_name_buffer_1;
	char* ptr_2 = file_name_buffer_2;

	strcat(ptr_1, filePath);
	Image image;

	for (int i = 0; i < n; ++i) {
		memset(ptr_2, 0, 256);
		sprintf(ptr_2, "%s-%d.jpeg", file_path_buffer, i);

		JpegUtils::readJpegFile(image, ptr_1);
		JpegUtils::writeJpegFile(image, ptr_2, quality);

		char* temp = ptr_1;
		ptr_1 = ptr_2;
		ptr_2 = temp;
	}
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
int StegUtils::doCoefficientDiff(const char* filePath1, const char* filePath2) {
	JpegCoefficients coefficients1;
	JpegUtils::readJpegCoefficients(filePath1, coefficients1);
	JpegCoefficients coefficients2;
	JpegUtils::readJpegCoefficients(filePath2, coefficients2);

	return StegUtils::diffCoefficients(coefficients1, coefficients2);
}

/**
 * this function deletes the intermediate files in a compression series
 *
 * @param jpegFilePath
 *     the original input path/name for the compression series (this file
 *     will not be deleted)
 * @param limit
 *     every file in the compression series (except the original image)
 *     whose n value is strictly less than limit will be deleted
 */
void StegUtils::purge(const char* jpegFilePath, int limit) {
	char jpegPathBuffer[256];
	Utils::getFileNameWithoutExtension(jpegFilePath, jpegPathBuffer);

	char filePathBuffer[256];
	for (int i = 0; i < limit; ++i) {
		memset(filePathBuffer, 0, 256);
		sprintf(filePathBuffer, "%s-%d.jpeg", jpegPathBuffer, i);
		remove(filePathBuffer);
	}
}

/**
 * this file returns the number of usable coefficients a lossless encode can expect when
 * encoding information.  this can be used to check whether payloads will fit into the image
 *
 * @param inputFileName
 *     the path to a jpeg image file to use
 *
 * returns the number of non-zero coefficients data can be embed into.  this number can
 *     be used to see whether data will fit into the image
 */
unsigned int StegUtils::getNumberOfCoefficients(const char *inputFileName) {
	JpegCoefficients inputJpeg;
	if (!JpegUtils::readJpegCoefficients(inputFileName, inputJpeg)) {
		printf("error loading jpeg file %s\n!", inputFileName);
		return -1;
	}

	int number_of_coefficients = 0;
	short coef = 0;

	for (JDIMENSION compnum = 0; compnum < inputJpeg.numComponents; ++compnum) {
		for (JDIMENSION rownum = 0; rownum < inputJpeg.heightInBlocks[compnum]; ++rownum) {
			for (JDIMENSION blocknum = 0; blocknum < inputJpeg.widthInBlocks[compnum]; ++blocknum) {
				for (JDIMENSION i = 0; i < DCTSIZE2; ++i) {
					coef = inputJpeg.coefBuffers[compnum][rownum][blocknum][i];

					// we only use non-zero coefficients
					if (coef != 0) {
						++number_of_coefficients;
					}
				}
			}
		}
	}

	return number_of_coefficients;
}

/*
 * this function loops through two coefficient store structs and returns the number of coefficients
 * between the two that differ.  this is a debug/tool function that probably wont be used in production
 *
 * @param coefficients1
 *     the coefficients from the 1st image in the diff
 * @param coefficients2
 *     the coefficients from the 2nd image in the diff
 *
 * returns the number of coefficients that differ between the two sets of coefficients
 */
int StegUtils::diffCoefficients(JpegCoefficients& coefficients1, JpegCoefficients& coefficients2) {
	int diff = 0;
	int count = 0;
	int num_components = coefficients1.info->num_components;
	size_t block_row_size[num_components];
	int width_in_blocks[num_components];
	int height_in_blocks[num_components];
	int x = 0;
	int y = 0;
	for (JDIMENSION compnum = 0; compnum < coefficients1.info->num_components; ++compnum) {
		for (JDIMENSION rownum = 0; rownum < height_in_blocks[compnum]; ++rownum) {
			for (JDIMENSION blocknum = 0; blocknum < width_in_blocks[compnum]; blocknum++) {
				x = 0;
				y = 0;
				for (JDIMENSION i = 0; i < DCTSIZE2; i++) {
					x = (x + 1) % DCTSIZE;
					if (x == 0) {
						++y;
					}

					int val1 = coefficients1.coefBuffers[compnum][rownum][blocknum][i];
					int val2 = coefficients2.coefBuffers[compnum][rownum][blocknum][i];

					if (val1 != val2) {
						std::cout << val1 << " != " << val2 << ", for " << x << ", " << y << " at rowblock "
								<< rownum << ", " << blocknum << "\n" << std::endl;
						++diff;
					}
				}
			}
		}
	}

	std::cout << "number of diffs = " << diff << std::endl;

	return diff;
}
