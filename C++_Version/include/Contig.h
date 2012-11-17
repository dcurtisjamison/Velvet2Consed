/***********************************************************************
 * Contig.h
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

#ifndef CONTIG_H_
#define CONTIG_H_

#include <string>
#include <set>
using namespace std;

#include "Sequence.h"
#include "Tile.h"
#include "Read.h"
#include "DatastoreHandler.h"

#include <foundation_AMOS.hh>
#include <boost/foreach.hpp>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/set.hpp>

struct myComp {
	bool operator() (const int& l, const int& r) const {return l<r;}
};

class DatastoreHandler;

class Contig: public Sequence {

	int amosID;
	string velvetID;
	multiset<Tile, Tile::tileComp> tiles;

	template <class T> inline std::string to_string (const T& t) {
		std::stringstream ss;
		ss << t;
		return ss.str();
	}
	friend class boost::serialization::access;
	template<class Archive> void serialize (Archive & ar, const unsigned int version)  {
		ar & boost::serialization::base_object<Sequence>(*this);
		ar & tiles;
		ar & amosID;
		ar & velvetID;
	}
	string consensus(string bases);
	string qualityscore(string qual);

public:
	Contig();
	Contig(AMOS::Contig_t acon);
	virtual ~Contig();

	void computeQuality(DatastoreHandler *dh);
//	void repairEnds(DatastoreHandler *dh, bool extend);
	int getAMOSID() {return amosID;}
	int readCount() {return tiles.size();}
	multiset<Tile, Tile::tileComp> getTileList() {return tiles;}
	string getACE();
};

#endif /* CONTIG_H_ */
