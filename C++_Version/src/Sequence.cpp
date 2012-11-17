/***********************************************************************
 * Sequence.cpp
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

#include "Sequence.h"

Sequence::Sequence() {
	// TODO Auto-generated constructor stub

}
Sequence::Sequence(string nam, string seq, string qual)
	: name(nam), sequence(seq), quality(qual),  sep("|")
{

}

Sequence::~Sequence() {
	// TODO Auto-generated destructor stub
}

string Sequence::toString() {
	stringstream ss;
	ss << name << sep << sequence << sep << quality;
	return ss.str();
}

void Sequence::setParam(int cnt, string param) {
	switch (cnt) {
	case 1:
		name = param;
		break;
	case 2:
		sequence = param;
		break;
	case 3:
		quality = param;
		break;
	default:
		cerr << "don\'t know what to do with " << cnt << " " << param << "\n";
		exit(0);
	}
}

