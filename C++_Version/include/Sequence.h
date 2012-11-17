/***********************************************************************
 * Sequence.h
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

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <iostream>
#include <sstream>
#include <string>
using namespace std;

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

class Sequence {
private:
	string name, sequence, quality, sep;
public:
	Sequence();
	Sequence(string nam, string seq, string qual);
	virtual ~Sequence();
	template<class Archive> void serialize (Archive & ar, const unsigned int version){
		ar & name;
		ar & sequence;
		ar & quality;
	}

	void setName(string nam) {name = name;}
	void setSequence(string seq) {sequence = seq;}
	void setQuality(string qual) {quality = qual;}
	void setParam(int cnt, string param);

	int length() {return sequence.length();}
	string getName() {return name;}
	string getSequence() {return sequence;}
	string getQuality() {return quality;}
	char getBase(int i) {return sequence.at(i);}
	char getBaseQuality(int i) {return quality.at(i);}

	string toString();
};

#endif /* SEQUENCE_H_ */
