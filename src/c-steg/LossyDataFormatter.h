//
//  lossy_data_formatter.h
//  c-steg
//

#ifndef CSTEG_DATA_FORMATTER_H
#define CSTEG_DATA_FORMATTER_H

#include "Utils.h"

struct packet_collection {
    unsigned int length;
    unsigned char* data;
};

class LossyDataFormatter {
public:
	static struct packet_collection* toLossy(const char* data, unsigned short dataLength);
	static void fromLossy(struct packet_collection* packetCollection, DecodeResponse& decodeResponse);

private:
	static unsigned short getNumberOfPackets(unsigned short dataLength);
	static bool increment_enu_1_bit(struct packet_collection_enumerator* enu);
	static bool decrement_enu_1_bit(struct packet_collection_enumerator* enu);
	static bool increment_enu_n_bytes(struct packet_collection_enumerator* enu, unsigned int n);
	static bool decrement_enu_n_bytes(struct packet_collection_enumerator* enu, unsigned int n);
	static unsigned char get_current_byte(struct packet_collection* packetCollection, struct packet_collection_enumerator* enu);
	static bool has_next_bit(struct packet_collection_enumerator* enu);
	static bool has_next_byte(struct packet_collection_enumerator* enu);
	static void insertPacket(struct packet* packet, struct node** nodeBuffer);
	static int get_remaining_bytes(struct packet_collection_enumerator* enu);
	static int read(struct packet_collection* packet_collection, struct packet_collection_enumerator* enu,
	          unsigned char* buffer, int numberOfPackets);
	static bool validatePacket(struct packet* packet);
	static bool isByteAligned(unsigned char* buffer, int numberOfPackets);
	static int handleFusterCluck(struct packet_collection* packetCollection, struct packet_collection_enumerator* enu);
	static unsigned short calculateDataLength(unsigned short* dataLengths, int dataLengthsSize);
	static void freeNodes(struct node** nodeBuffer);
	static unsigned char calculateValue(struct node* node, int bytesOffset);
	static void reconstructOrinialDataFromPackets(struct node** nodeBuffer, unsigned char* dataOut, unsigned int dataOutLength);
};

#endif
