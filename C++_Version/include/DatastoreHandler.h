/***********************************************************************
 * DatastoreHandler.h
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

#ifndef DATASTOREHANDLER_H_
#define DATASTOREHANDLER_H_

#include <string.h>
#include <fstream>
using namespace std;

#include "Contig.h"
#include "Read.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>

#include <boost/interprocess/file_mapping.hpp>
#include <boost/filesystem.hpp>
namespace BFS = boost::filesystem;
class Contig;

class DatastoreHandler {
private:
	BFS::path readIndex, readStore, contigIndex, contigStore;
	template <class T> inline std::string to_string (const T& t) {
		std::stringstream ss;
		ss << t;
		return ss.str();
	}

	fstream riFile, rsFile, ciFile, csFile;
	map<int,int> rMap, cMap;
	void makeMap (map<int,int>* myMap, fstream* IDX);

public:
	static const int CREATE = 0;
	static const int APPEND = 1;
	static const int FETCH = 2;
	DatastoreHandler(string dsDir, int direction);
	virtual ~DatastoreHandler();
	void storeContig(Contig ctg);
	void storeRead(Read read);
 	Contig getContig(int amosID);
	Read getRead(int amosID);
	vector<int> getContigList ();
};

#endif /* DATASTOREHANDLER_H_ */
