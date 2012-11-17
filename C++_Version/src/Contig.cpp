/***********************************************************************
 * Contig.cpp
 * Author: jam5pi
 *   Date: Oct 13, 2009
 *
 * <software purpose to be completed>
 ***********************************************************************
 * Revision History
 *
 ***********************************************************************
 * Copyright 2009, Children's Hospital of Cincinnati
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ***********************************************************************/

#include "Contig.h"

Contig::Contig() {
	// TODO Auto-generated constructor stub

}

Contig::Contig(AMOS::Contig_t actg)
	: Sequence(to_string(actg.getEID()), actg.getSeqString(), actg.getQualString())
{
	amosID = actg.getIID();
	velvetID = actg.getEID();
	BOOST_FOREACH(AMOS::Tile_t aTile, actg.getReadTiling()) {
		Tile myTile(aTile);
		tiles.insert(myTile);
	}
}

Contig::~Contig() {
	// TODO Auto-generated destructor stub
}

string Contig::getACE() {
	stringstream ss;
	ss << "CO ctg_" << amosID << " " << length() << " " << readCount() << " 1 U\n";
	string temp = getSequence();
	int i = 0;
	while (i < length()) {
		ss << temp.substr(i,70) << "\n";
		i += 70;
	}
	ss << "\nBQ\n";
	i = 0;
	while (i < length()) {
		int q = (int)getBaseQuality(i) - 32;
		ss << q << " ";
		i++;
		if (i%70 == 0) {
			ss << "\n";
		}
	}
	return ss.str();
}

void Contig::computeQuality(DatastoreHandler *dh) {
	int currOffset = 0;
	deque<string> qual(length(), "");
	BOOST_FOREACH(Tile t, tiles) {
		Read r = dh->getRead(t.getAMOSID());
		int tOffset = t.getOffset();
		if (tOffset - currOffset < 0) {
			for(int i = tOffset; i <= currOffset; i++) {
				qual.push_front(string(""));
			}
			currOffset = tOffset;
		}
		int tLast = tOffset + r.length();
		if (tLast > length()) {
			for(int i = length(); i <= tLast; i++) {
				qual.push_back(string(""));
			}
		}
		deque<string>::iterator currQualRow = qual.begin() + tOffset - currOffset;
		int start = t.getReadStart();
		int stop = t.getReadEnd();
		int dir = 1;
		if (start > stop) {
			// turn the step around and predecrement
			dir = -1;
			start -= 1;
			stop -= 1;
		}
		for (int c = start; c != stop; c += dir) {
			char q = r.getBaseQuality(c);
			currQualRow->append(to_string(q));
			currQualRow++;
		}
	}
	string newQual;
	BOOST_FOREACH(string s, qual) {
		newQual.append(qualityscore(s));
	}

	setQuality(newQual.substr(0 - currOffset, length()));
}

//void Contig::repairEnds(DatastoreHandler *dh, bool extend) {
//
//}

string Contig::qualityscore(string qual) {
	char nullScore = char(Read::fastq_pad);
	int count = 0;
	int total = 0;
	for(int c = 1; c < (int)qual.length(); c++) {
		char q = qual.at(c);
		total += (int(q) - Read::fastq_pad);
		count++;
	}
	char qcode = nullScore;
	if (count > 0) {
		qcode = char(int(total/count) + Read::fastq_pad);
	}
	return to_string(qcode);
}

string Contig::consensus(string bases) {
	int cnt[5];
	for (int i = 0; i < (int)bases.length(); i++) {
		switch (bases.at(i)) {
		case 'A':
		case 'a':
			cnt[0]++;
			break;
		case 'C':
		case 'c':
			cnt[1]++;
			break;
		case 'G':
		case 'g':
			cnt[2]++;
			break;
		case 'T':
		case 't':
			cnt[3]++;
			break;
		default :
			cnt[4]++;
		}
	}
	int maxdex = 5;
	for (int i = maxdex; i >= 0; i--) {
		if (cnt[i] > cnt[maxdex]) {
			maxdex = i;
		}
	}
	switch (maxdex) {
	case 0: return string("A");
	case 1: return string("C");
	case 2: return string("G");
	case 3: return string("T");
	}
	return string("N");
}
