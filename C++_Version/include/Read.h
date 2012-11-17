/***********************************************************************
 * Read.h
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

#ifndef READ_H_
#define READ_H_

#include <iostream>
#include <sstream>
#include <string>
#include <time.h>
using namespace std;
#include "Sequence.h"

#include <boost/tokenizer.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

class Read: public Sequence {
private:
	int amosID;
	string velID, sep;
	string now;

public:
	static const int fastq_pad = 64;
	Read();
	Read(string nam, string seq, string qual);
	virtual ~Read();

	template<class Archive> void serialize (Archive & ar, const unsigned int version)  {
		ar & boost::serialization::base_object<Sequence>(*this);
		ar & amosID;
		ar & velID;
		ar & now;
		ar & sep;
	}

	void setVelvetID(string id) {velID = id;}
	void setAMOSID(int id) {amosID = id;}

	int getAMOSID() {return amosID;}
	string getVelvetID() {return velID;}

	string toString();
	string getPHD();
	string getAF();
	string getRD();

};

#endif /* READ_H_ */
