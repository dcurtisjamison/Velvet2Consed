/***********************************************************************
 * DatastoreHandler.cpp
 * Author: jam5pi
 *   Date: Oct 13, 2009
 *
 * An interface to compressed assembly components
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

#include "DatastoreHandler.h"
	const int DatastoreHandler::CREATE;
	const int DatastoreHandler::FETCH;
	const int DatastoreHandler::APPEND;
/*
 * constructor
 */
DatastoreHandler::DatastoreHandler(string dsDir, int direction = FETCH) {
	string ridx = dsDir;
	ridx += "/reads.idx";
	string rdsb = dsDir;
	rdsb += "/reads.dsb";
	string cidx = dsDir;
	cidx += "/contigs.idx";
	string cdsb = dsDir;
	cdsb += "/contigs.dsb";

	//we have a path for all the files
	readIndex = BFS::path(ridx);
	readStore = BFS::path(rdsb);
	contigIndex = BFS::path(cidx);
	contigStore = BFS::path(cdsb);

	//open and check the files
	ios::openmode mode;
	switch (direction) {
	case CREATE:
		mode = ios::out;
		break;
	case APPEND:
		mode = ios::out|ios::app|ios::ate;
		break;
	case FETCH:
	default:
		mode = ios::in;
		break;
	}
	riFile.open (readIndex.string().c_str(), mode);
	rsFile.open (readStore.string().c_str(), mode);
	ciFile.open (contigIndex.string().c_str(), mode);
	csFile.open (contigStore.string().c_str(), mode);

	if (direction == FETCH) {
		makeMap(&rMap, &riFile);
		makeMap(&cMap, &ciFile);
	}
}

DatastoreHandler::~DatastoreHandler() {
	// TODO Auto-generated destructor stub
	riFile.close();
	ciFile.close();
	rsFile.close();
	csFile.close();
}

void DatastoreHandler::storeContig(Contig ctg) {
	int offset = csFile.tellp();
	boost::archive::text_oarchive oa(csFile);
	oa << ctg;
	// write the index
	cerr << "(ciFile) gets " << ctg.getAMOSID() << " " << ctg.getName() << " " << offset << "\n";
	ciFile << ctg.getAMOSID() << " " << ctg.getName() << " " << offset << "\n";
	cout << ctg.getAMOSID() << " " << ctg.getName() << " " << offset << "\n";

}

void DatastoreHandler::storeRead(Read read) {
	int offset = rsFile.tellp();
	boost::archive::text_oarchive oa(rsFile);
	oa << read;
	//write the index
	riFile << read.getAMOSID() << " " << read.getName() << " " << offset << "\n";
}

Contig DatastoreHandler::getContig(int amosID) {
	Contig ctg;
	int offset = cMap[amosID];
	csFile.seekg(offset);
	boost::archive::text_iarchive ia(csFile);
	ia >> ctg;
	return ctg;
}

Read DatastoreHandler::getRead(int amosID) {
	Read read;
	// as a test, let's try fetching a read using the AMOSID as the offset
	int offset = rMap[amosID];
	rsFile.seekg(offset);
	boost::archive::text_iarchive ia(rsFile);
	ia >> read;
	return read;
}

void DatastoreHandler::makeMap (map<int,int>* myMap, fstream* IDX) {
	string line;
	char buff[256];
	IDX->seekg(0);
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	boost::char_separator<char> sep(" ");
	IDX->getline(buff,256);
	while (! IDX->eof()) {
		line = to_string(buff);
		tokenizer tok(line, sep);
		tokenizer::iterator word = tok.begin();
		int id = atoi(word->c_str());
		word++;
		word++;
		int offset = atoi(word->c_str());
		pair<int,int> vt(id,offset);
		myMap->insert(vt);
		IDX->getline(buff,256);
	}
}

vector<int> DatastoreHandler::getContigList () {
	vector<int> idList;
	for (map<int,int>::iterator i = cMap.begin(); i != cMap.end(); i++) {
		idList.push_back(i->first);
	}
	return idList;
}
