/***********************************************************************
 * Tile.h
 * Author: jam5pi
 *   Date: Nov 9, 2009
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

#ifndef TILE_H_
#define TILE_H_
#include <sstream>
using namespace std;

#include <foundation_AMOS.hh>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>


class Tile {
private:
	int src, off, clr1, clr2;

public:
	Tile();
	Tile(AMOS::Tile_t aTile);
	virtual ~Tile();
	template<class Archive> void serialize (Archive & ar, const unsigned int version) {
		ar & src;
		ar & off;
		ar & clr1;
		ar & clr2;
	}

	struct tileComp {//: binary_function<Tile*, Tile*, bool> {
		bool operator() (const Tile& p, const Tile& q) {
			return p.off < q.off;
		}
	};// comp;


	int getAMOSID() { return src;}
	int getOffset() {return off;}
	int getReadStart() {return clr1;}
	int getReadEnd() {return clr2;}
	bool isForwardRead() {return clr1 < clr2;}
	string getAF();

};

#endif /* TILE_H_ */
