//
//  lossy_data_formatter.cpp
//  c-steg
//
//  this file contains functions that format and retrieve data into
//  packets that can survive bits flipping or dissapearing
//

#include "LossyDataFormatter.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#define DATA_FORMATTER_PACKET_LENGTH 4
#define NUMBER_OF_NECESSARY_VALID_PACKETS 6
#define NUMBER_OF_PACKETS_FOR_BYTE_ALIGNMENT 2
#define NUMBER_OF_BITS_TO_REWIND 4
#define TERMINATE -1
#define BAD_PACKET 1
#define FIXED_IT 2
#define BYTE_ALIGNMENT_PROBLEM -2
#define MAX_LOSSY_LENGTH 4096

struct packet {
    unsigned short dataLength;
    unsigned char bytes[DATA_FORMATTER_PACKET_LENGTH];
    unsigned char index;
    unsigned char sum;
    unsigned char checkSum;
};

struct packet_collection_enumerator {
    int currentBit;
    int currentByte;
    unsigned int packetCollectionLength;
};

/**
 * this struct is used in the linked lists for each packet index to
 * keep track of the the majority of what has come in for a particular
 * packet index.  in the event there is some corruption taking the
 * majority vote for a packet helps fight the corruption
 */
struct node {
    unsigned char bytes[DATA_FORMATTER_PACKET_LENGTH];
    struct node* next;
};


class BitEnumerator {
public:
static bool incrementEnu1Bit(struct packet_collection_enumerator* enu);
static bool decrementEnu1Bit(struct packet_collection_enumerator* enu);
static bool incrementEnuNBytes(struct packet_collection_enumerator* enu, unsigned int n);
static bool decrementEnuNBytes(struct packet_collection_enumerator* enu, unsigned int n);
static unsigned char getCurrentByte(struct packet_collection* packetCollection,
		struct packet_collection_enumerator* enu);
static bool hasNextBit(struct packet_collection_enumerator* enu);
static bool hasNextByte(struct packet_collection_enumerator* enu);
static int getRemainingBytes(struct packet_collection_enumerator* enu);
};

/**
 * this function increments the bit stream enumerator for a packet collection 1 bit
 *
 * @param enu
 *     the bit stream enumerator to operate on
 *
 * returns true if the enumerator is left in a valid state (meaning it has more bits) 
 */
bool BitEnumerator::incrementEnu1Bit(struct packet_collection_enumerator* enu) {
    if (enu->currentByte < enu->packetCollectionLength) {
        ++enu->currentBit;
        if (enu->currentBit % 8 == 0) {
            enu->currentBit = 0;
            ++enu->currentByte;
        }
        
        return enu->currentByte < enu->packetCollectionLength;
    } else {
        return false;
    }
}

/*
 * this function decrements the bit stream enumerator for a packet collection 1 bit
 *
 * @param enu
 *     the bit stream enumerator to operate on
 *
 * returns true if the enumerator was left in a valid state (meaning the byte index
 *     is greater than -1). the enumerator can be positioned before the first bit
 *     and that is considered an invalid position.
 */
bool BitEnumerator::decrementEnu1Bit(struct packet_collection_enumerator* enu) {
    if (enu->currentByte > -1) {
        --enu->currentBit;
        if (enu->currentBit < 0) {
            enu->currentBit = 7;
            --enu->currentByte;
        }
        
        return enu->currentByte > -1;
    } else {
        return false;
    }
}

/*
 * this function increments the bit stream enumerator for a packet collection n bytes
 *
 * @param enu
 *     the bit stream enumerator to operate on
 * @param n
 *     the number of bytes to increment the bit stream operator
 *
 * returns true if there are more bytes after the increment operation completed
 */
bool BitEnumerator::incrementEnuNBytes(struct packet_collection_enumerator* enu, unsigned int n) {
    enu->currentByte += n;
    if (enu->currentByte >= enu->packetCollectionLength) {
        enu->currentBit = 0;
        enu->currentByte = enu->packetCollectionLength;
        return false;
    } else {
        return true;
    }
}

/*
 * this function decrements the bit stream enumerator for a pack collection n bytes
 *
 * @param enu
 *     the bit stream enumerator to operate on
 * @param n
 *     the number of bytes to decrement the bit stream operator
 *
 * returns true if the enumerator byte index is greater than or equal to 0 after the 
 *     decrement operation
 */
bool BitEnumerator::decrementEnuNBytes(struct packet_collection_enumerator* enu, unsigned int n) {
    enu->currentByte -= n;
    if (enu->currentByte < 0) {
        enu->currentBit = 7;
        enu->currentByte = -1;
        return false;
    } else {
        return true;
    }
}

/*
 * get the current byte out of the packet collection pointed to by a bit stream enumerator
 *
 * @param packetCollection
 *     the packet collection from which to retrieve a byte
 * @param enu
 *     the bit stream enumerator pointing to a byte to retrieve
 *
 * returns the unsigned char from the packet collection pointed to be the bit stream enumerator
 */
unsigned char BitEnumerator::getCurrentByte(struct packet_collection* packetCollection,
		struct packet_collection_enumerator* enu) {
    return 
        (packetCollection->data[enu->currentByte] >> enu->currentBit) |
        ((packetCollection->data[enu->currentByte + 1] << (8 - enu->currentBit)) & 255);
}

/*
 * this function can be used to see if a bit stream enumerator has more bits
 *
 * @param enu
 *     the bit stream enumerator to check
 *
 * returns true if the bit stream enumerator has more bits
 */
bool BitEnumerator::hasNextBit(struct packet_collection_enumerator* enu) {
    return enu->currentByte < enu->packetCollectionLength;
}

/* 
 * this function can be used to see if a bit stream enumerator has more bytes
 *
 * @param enu
 *     the bit stream enumerator to check
 *
 * returns true if the bit stream enumerator has more bytes
 */
bool BitEnumerator::hasNextByte(struct packet_collection_enumerator* enu) {
    return 
        enu->currentByte < enu->packetCollectionLength - 1 ||
        (enu->currentByte == enu->packetCollectionLength - 1 && enu->currentBit == 0);
}

/**
 * this function calculates and returns the remaining bytes in a bit 
 * stream enumerator
 *
 * @param enu
 *     the bit stream enumerator to use
 *
 * returns the number of bytes remaining in the bit stream enumerator
 */
int BitEnumerator::getRemainingBytes(struct packet_collection_enumerator* enu) {
    int rtn = enu->currentByte;
    if (enu->currentBit != 0) {
        --rtn;
    }

    rtn = enu->packetCollectionLength - rtn;

    if (rtn < 0) { 
        return 0; 
    } else { 
        return rtn;
    }
}

/**
 * this function returns the number of packets needed for transmitting
 * a data payload.
 *
 * @param dataLength
 *     the length of the data being wrapped in packets
 */
unsigned short LossyDataFormatter::getNumberOfPackets(unsigned short dataLength) {
    int numberOfPackets = (dataLength / DATA_FORMATTER_PACKET_LENGTH);
    int numberOfPacketsRemainder = dataLength % DATA_FORMATTER_PACKET_LENGTH;

    if (numberOfPacketsRemainder != 0) {
        ++numberOfPackets;
    }

    return numberOfPackets;
}

/**
 * this function wraps the input data into packets and puts them in a
 * packet collection.  the packet collection is more resilient to the data
 * payload getting corrupted during transmission.
 *
 * @param data
 *     the data to wrap into packets
 * @param dataLength
 *     the length of the data being put into packets
 */
packet_collection* LossyDataFormatter::toLossy(const char* data, unsigned short dataLength) {
    time_t t;
    srand((unsigned) time(&t));

    unsigned short numberOfPackets = getNumberOfPackets(dataLength);

    // each byte in data needs to be included in the length, plus
    // 1 check sum per packet, plus the space for storing the packet
    // length multiple times and the xors for the length (1 per time
    // length is represented)
    unsigned int lossyLength = ((numberOfPackets) * sizeof(struct packet));

    struct packet_collection* packetCollection =
        (struct packet_collection*) malloc(sizeof(struct packet_collection));
    packetCollection->length = lossyLength * sizeof(unsigned char);
    packetCollection->data =
        (unsigned char*) malloc(packetCollection->length * sizeof(unsigned char));

    memset(packetCollection->data, 0, packetCollection->length * sizeof(unsigned char));

    int index = 0;
    struct packet* packets = (struct packet*)&packetCollection->data[0];
    unsigned char packetXor;
    unsigned char packetSum;
    for (int i = 0; i < numberOfPackets; ++i) {
        packetXor = 0;
        packetSum = 0;
        for (int j = 0; j < DATA_FORMATTER_PACKET_LENGTH; ++j) {
            if (index < dataLength) {
                packets[i].bytes[j] = data[index++];
            } else {
                packets[i].bytes[j] = 0;
            }

            packetXor ^= packets[i].bytes[j];
            packetSum += packets[i].bytes[j];
        }

        packets[i].dataLength = dataLength;
        packets[i].index = (i % 256);
        packets[i].sum = (packetSum + packets[i].index) % 256;
        packets[i].checkSum = (packetXor ^ packets[i].index);
    }

    return packetCollection;
}

/**
 * this function inserts a packet into the node buffer.  the node
 * buffer is an array of linked lists.  each packet gets inserted into
 * linked list corresponding to its packet index
 *
 * @param packet
 *     the extracted packet to insert into the node_buffer
 * @param nodeBuffer
 *     the array of node linked lists to insert a packet into
 */
void LossyDataFormatter::insertPacket(struct packet* packet, struct node** nodeBuffer) {
    struct node* newNode = (struct node*) malloc(sizeof(struct node));
    memcpy(newNode->bytes, packet->bytes, DATA_FORMATTER_PACKET_LENGTH);
    //packet->index is an unsigned char and
    //node_buffer is of length 256 so no bounds
    //checking necessary
    newNode->next = nodeBuffer[packet->index];
    nodeBuffer[packet->index] = newNode;
}

/*
 * this function reads a specified number of packets from the packet collection and
 * places them in the buffer.
 *
 * @param packetCollection
 *     the packet collection containing the payload of packets to unpack
 * @param enu
 *     the bit stream enumerator being used for the from lossy operation
 * @param buffer
 *     the buffer to store the the read packets into
 * @param numberOfPackets
 *     the number of packets to read from the packet collection
 *
 * returns true if the packets were able to be read, false if the
 *     requested number of packets exceeds the number of remaining 
 *     bytes in the packet collection
 */
int LossyDataFormatter::read(struct packet_collection* packetCollection, struct packet_collection_enumerator* enu,
          unsigned char* buffer, int numberOfPackets) {
    if (numberOfPackets * sizeof(struct packet) > BitEnumerator::getRemainingBytes(enu)) {
        return 0;
    }

    //this should not happen but if the is some bug 
    //somewhere that lets enu.current_byte be less than
    //0 this will stop from reading invalid memory
    if (enu->currentByte < 0) {
        enu->currentByte = 0;
    }

    for (int h = 0; h < numberOfPackets; ++h) {
        for (int i = 0; i < sizeof(struct packet); ++i) {
            if (!BitEnumerator::hasNextByte(enu)) {
                return 0;
            }
            
            buffer[(h * sizeof(struct packet)) + i] = BitEnumerator::getCurrentByte(packetCollection, enu);
            BitEnumerator::incrementEnuNBytes(enu, 1);
        }
    }

    return 1;
}

/**
 * this function validates a packet. if a bit stream becomes corrupted by dropping
 * bits or flipping bits packets can become invalid.  the checks performed are an
 * xor comparison against a check sum and a rolling sum of the byte payload
 *
 * @param packet
 *     the packet to check for validity
 *
 * returns true if the packet passes validity checks, false otherwise
 */
bool LossyDataFormatter::validatePacket(struct packet* packet) {
    unsigned char packetXor = 0;
    unsigned char packetSum = 0;
    for (int i = 0; i < DATA_FORMATTER_PACKET_LENGTH; ++i) {
        packetXor ^= packet->bytes[i];
        packetSum += packet->bytes[i];
    }

    return 
        (packet->sum == (packetSum + packet->index) % 256) &&
        (packet->checkSum ^ packetXor ^ packet->index == 0);
}

/**
 * this method checks a packet buffer containing multiple packets to see if
 * enough of the packets checkout that it would be prudent to assume no bits
 * were dropped from the packet buffer and it is byte aligned.  a bit being
 * dropped should throw off a large number of packets as opposed a bit being
 * flipped which would only throw off one packet.
 *
 * @param buffer
 *     the buffer containing packets (a byte buffer so is stores the individual
 *     bytes of all the packets sequentially)
 * @param numberOfPackets
 *     the number of packets the buffer contains
 *
 * returns true if the buffer should be considered byte aligned, false otherwise
 */
bool LossyDataFormatter::isByteAligned(unsigned char* buffer, int numberOfPackets) {
    struct packet* packet = (struct packet*)&buffer[0];
    int numberOfValidPackets = 0;
    int numberOfMatchingIndexes = 0;
    int packetIndex = packet[0].index;
    for (int i = 0; i < numberOfPackets; ++i) {
        if (validatePacket(&packet[i])) {
            ++numberOfValidPackets;
        }

        if (i > 0) {
            if (packet[i].index - packet[i - 1].index == 1) {
                ++numberOfMatchingIndexes;
            }
        } else {
            if (packet[i + 1].index - packet[i].index == 1) {
                ++numberOfMatchingIndexes;
            }
        }
    }

    if (numberOfValidPackets >= NUMBER_OF_PACKETS_FOR_BYTE_ALIGNMENT) {
        return true;
    } else if (numberOfPackets >= 2 && numberOfValidPackets >= (numberOfPackets / 2)) {
        return true;
    } else if (numberOfPackets == 1 && numberOfValidPackets == 1) {
        return true;
    } else {
        return false;
    }
}

/**
 * this method handles the situations that arise when a packet is invalid.  a packet can be invalid
 * due to a bit flipping or a bit dropping.  this function tries to figure out which one of those cases
 * happened. if the buffer is not considered to be byte aligned this function will shift backwards 1 bit
 * at a time to try to become byte aligned.
 *
 * @param packetCollection
 *     the packet collection holding the payload being unpacked
 * @param enu
 *     the bit stream enumerator being used for the unpacking of the packet collection
 *
 * returns some constant from the h file if byte alignment can be reached FIXED_IT will be returned,
 *     if there are not enough packets left in the packet collection TERMINATE is returned, or if 
 *     byte alignment can not be reached BYTE_ALIGNMENT_PROBLEM is returned
 */
int LossyDataFormatter::handleFusterCluck(struct packet_collection* packetCollection, struct packet_collection_enumerator* enu) {
    packet_collection_enumerator enuHold;
    packet_collection_enumerator enuStart = *enu;

    int remainingPackets = BitEnumerator::getRemainingBytes(enu) / sizeof(struct packet);
    if (remainingPackets <= 0) {
        return TERMINATE;
    }

    int numberOfPackets =
        remainingPackets < NUMBER_OF_NECESSARY_VALID_PACKETS ?
        remainingPackets :
        NUMBER_OF_NECESSARY_VALID_PACKETS;

    unsigned char buffer[sizeof(struct packet) * numberOfPackets];
     
    // this for loop will read a number of packets, then 
    // check them for byte alignment. if they are not aligned
    // the bit stream enumerator will be decremented one bit
    // NUMBER_OF_BITS_TO_REWIND controls how many bits we 
    // will try to rewind, which correlates to how many
    // bit we think could be dropped in a single packet
    for (int i = 0; i < NUMBER_OF_BITS_TO_REWIND; ++i) {
        // remember where we are before the read, the read is only to
        // determine byte alignment, not to actually advance the enumerator
        enuHold = *enu;
        
        // if read is unsuccessful that means we are at the end of the packet collection
        // alert the calling function to stop
        if (!read(packetCollection, enu, buffer, NUMBER_OF_NECESSARY_VALID_PACKETS)) {
            return TERMINATE;
        }
        
        // put the enumerator back where it was prior to the read
        *enu = enuHold;

        // if we are byte aligned then return, the problem was
        // just a flipped bit
        if (isByteAligned(buffer, numberOfPackets)) {
            return FIXED_IT;
        } 
        
        // if byte alignment could not be determined, move the
        // bit stream enumerator back one to accommodate a dropped
        // bit
        BitEnumerator::decrementEnu1Bit(enu);
    }

    *enu = enuStart;
    return BYTE_ALIGNMENT_PROBLEM;
}

/**
 * the data length is a short and is embedded in every packet.  this way
 * if the length is corrupted in one packet i can make it through in another.  
 * the assumption here is than it will make it through more often than it will 
 * be corrupted.  the when unpacking the length that occurs the most times is
 * considered the actual length.  this function looks for the length that occurred
 * the most.
 *
 * @param dataLengths
 *     the unsigned short array of lengths that came through the packets that
 *     were unpacked
 * @param dataLengthsSize
 *     the number of length values in the data_length array
 *
 * returns the length that occurred the most times in the data_lenghts array
 */
unsigned short LossyDataFormatter::calculateDataLength(unsigned short* dataLengths, int dataLengthsSize) {
    if (dataLengthsSize <= 0) {
        return 0;
    }

    unsigned short uniqueValues[256];
    memset(uniqueValues, 0, 256 * sizeof(unsigned short));
    int uniqueCounts[256];
    memset(uniqueCounts, 0, 256 * sizeof(int));
    int uniqueValuesSize = 0;

    for (int i = 0; i < dataLengthsSize; ++i) {
        int contains = 0;
        for (int j = 0; j < uniqueValuesSize; ++j) {
            if (uniqueValues[j] == dataLengths[i]) {
                contains = 1;
                ++uniqueCounts[j];
                break;
            }
        }

        // we only look at the first 256 uniquely occurring
        // length values.  that is way more than should
        // ever occur but if for some reason more than that
        // occur (or even anything remotely close to that occur)
        // something has gone wrong.  in practice there should be a 
        // relatively small number of lengths that occur (less than 10
        // probably)    
        if (!contains && uniqueValuesSize < 256) {
            uniqueValues[uniqueValuesSize] = dataLengths[i];
            uniqueCounts[++uniqueValuesSize] = 1;
        }
    }

    // look through the unique value histogram and find
    // the value that occurred the most
    unsigned short maxIndex = 0;
    int maxValue = uniqueCounts[0];
    for (int i = 1; i < uniqueValuesSize; ++i) {
        if (uniqueCounts[i] > maxValue) {
            maxIndex = i;
            maxValue = uniqueCounts[i];
        }
    }

    // return the unique value whose count was maximal
    return uniqueValues[maxIndex];
}

/**
 * this function reclaims memory allocated in the linked lists of the node buffer 
 * array of linked lists
 *
 * @param nodeBuffer
 *     a node buffer used when unpacking a lossy data stream
 */
void LossyDataFormatter::freeNodes(struct node** nodeBuffer) {
    for (int i = 0; i < 256; ++i) {
        struct node* node = nodeBuffer[i];
        struct node* hold;
        while (node != NULL) {
            hold = node->next;
            free(node);
            node = hold;
        }
    }
}

/**
 * this function looks over the supplied linked list and logs how many unique
 * values occurred in the data stream for the character at a position of the
 * data stream.  return the character that occurs the most in the linked list
 *
 * @param node
 *     a linked list of packets that occurred while unpacking a packet collection
 * @param bytesOffset
 * 	   the offset of the byte in the node to look at
 * returns
 *     the bytes position to look at in the packet linked list, this corresponds to
 *     the character being sought for the data stream being unpacked.
 */
unsigned char LossyDataFormatter::calculateValue(struct node* node, int bytesOffset) {
    if (node == NULL) {
        return 0;
    }

    unsigned char uniqueValues[256];
    int uniqueValuesSize = 0;
    int uniqueCounts[256];
    
    // loop through the values that occurred for this character position 
    // in the data stream and keep a histogram of how many times a value 
    // occurred.  ultimately we want to return the value that occured the most
    while (node != NULL) {
        int contains = 0;
        for (int j = 0; j < uniqueValuesSize; ++j) {
            if (uniqueValues[j] == node->bytes[bytesOffset]) {
                contains = 1;
                ++uniqueCounts[j];
                break;
            }
        }

        // only look at the 1st 256 unique values for this character in the 
        // data stream.  that should be way more than enough
        if (!contains && uniqueValuesSize < 256) {
            uniqueValues[uniqueValuesSize] = node->bytes[bytesOffset];
            uniqueCounts[++uniqueValuesSize] = 1;
        }

        node = node->next;
    }

    // look through the unique values histogram and find which value occurred
    // the most number of times
    unsigned short maxIndex = 0;
    unsigned char maxValue = uniqueCounts[0];
    for (int i = 1; i < uniqueValuesSize; ++i) {
        if (uniqueCounts[i] > maxValue) {
            maxIndex = i;
            maxValue = uniqueCounts[i];
        }
    }
    
    //return the value that occurred the most number of times
    return uniqueValues[maxIndex];
}

/**
 * this function iterates through the node buffer and gets the character that occurred the most
 * at every position in the output data stream while unpacking the packet collection
 *
 * @param nodeBuffer
 *     the node buffer the packets were stored in
 * @param dataOut
 *     the buffer to store the output data stream to
 * @param dataOutLength
 *     the length of the dataOut buffer
 */
void LossyDataFormatter::reconstructOrinialDataFromPackets(struct node** nodeBuffer, unsigned char* dataOut, unsigned int dataOutLength) {
    unsigned char numberOfPackets = getNumberOfPackets(dataOutLength);
    for (int i = 0; i < numberOfPackets; ++i) {
        for (int j = 0; j < DATA_FORMATTER_PACKET_LENGTH; ++j) {
            unsigned short data_index = (i * DATA_FORMATTER_PACKET_LENGTH) + j;
            if (data_index < dataOutLength && nodeBuffer[i] != NULL) {
                dataOut[data_index] = calculateValue(nodeBuffer[i], j);
            }
        }
    }
}

/**
 * this function unpacks a packet_collection and produces the original input data that was packed into it. 
 * it is assumed the packet collection could have become corrupted and this function will try to not let any
 * corruption stop it from reconstructing the original data.
 *
 * @param packetCollection
 *     a packet collection that was deserialized from a jpeg image
 * @param decodeResponse
 *     an object to store the reassembled data in 
 */
void LossyDataFormatter::fromLossy(struct packet_collection* packetCollection, DecodeResponse& decodeResponse) {
    struct node* nodeBuffer[256];
    memset(nodeBuffer, 0, 256 * sizeof(struct node*));

    unsigned short dataLengths[256];
    memset(dataLengths, 0, 256 * sizeof(unsigned short));    
    int dataLengthsSize = -1;

    unsigned char* payload;
    struct packet_collection_enumerator enu;
    enu.currentBit = 0;
    enu.currentByte = 0;
    enu.packetCollectionLength = packetCollection->length;

    unsigned char buffer[sizeof(struct packet)];
    struct packet* packet = (struct packet*)&buffer[0];

    while (read(packetCollection, &enu, buffer, 1)) {
        int packetIsValid = validatePacket(packet);
        
        // if packet is not valid attempt to figure out why
        // if a bit was dropped handle_fuster_cluck will position
        // the enu to account for the dropped bit
        if (!packetIsValid) {
            int fCluckReturn = handleFusterCluck(packetCollection, &enu);
            if (fCluckReturn == TERMINATE) {
                break;
            }
        } else {
            // if the packet is valid put it in the node buffer at the index
            // of the packet (packet->index). 
            insertPacket(packet, nodeBuffer);
            // keep track of the data lengths we get from the packets
            // data lengths is the size in chars of the original encoded
            // data
            if (dataLengthsSize < 255) {
                dataLengths[++dataLengthsSize] = packet->dataLength;
            }
        }
    }

    // if we did not get any packets that were valid
    if (dataLengthsSize < 1) {
        decodeResponse.data = NULL;
        decodeResponse.length = 0;
        freeNodes(nodeBuffer);
        return;
    }

    unsigned short dataLength = calculateDataLength(dataLengths, dataLengthsSize);    
    // this check prevent allocating ridiculous amounts of memory as a result of an
    // improperly unpacked packet collection.  we should not be encoding anything in 
    // the lossy format that exceeds this amount of memory
    if (dataLength > MAX_LOSSY_LENGTH) {
    	decodeResponse.data = NULL;
    	decodeResponse.length = 0;
        freeNodes(nodeBuffer);
        return;
    }
    
    // allocate the memory for the returned original data
    unsigned char* data = (unsigned char*) malloc(dataLength * sizeof(unsigned char));
    memset(data, 0, dataLength * sizeof(unsigned char));

    // reconstruct the original data from the unpacked packets
    reconstructOrinialDataFromPackets(nodeBuffer, data, dataLength);
    
    //free memory, be responsible
    freeNodes(nodeBuffer);

    decodeResponse.data = data;
    decodeResponse.length = dataLength;
}
