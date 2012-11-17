/***********************************************************************
 * Read.cpp
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

#include "Read.h"

Read::Read() {
	// TODO Auto-generated constructor stub

}

/*
 * standard constructor
 */
Read::Read(string nam, string seq, string qual)
	: Sequence(nam, seq, qual),  sep("|")
{
	string rdays[7] = {"Sun","Mon","Tues","Wed","Thurs","Fri","Sat"};
	string rmonths[12] = {"Jan","Feb","Mar","Apr","May","Jun",
			"Jul","Aug","Sep","Oct","Nov","Dec"};
	time_t t = time(0);
	tm ts = *localtime(&t);
	stringstream ss;
	ss << rdays[ts.tm_wday] << " " << rmonths[ts.tm_mon] << " " << ts.tm_hour <<
			":" << ts.tm_min << ":" << ts.tm_sec << " " << ts.tm_year;
	now = ss.str();
}

Read::~Read() {
	// TODO Auto-generated destructor stub
}

string Read::toString() {
	stringstream ss;
	ss << amosID << sep << velID << sep << Sequence::toString() ;
	return ss.str();
}

string Read::getPHD() {
	stringstream ss;

	int lind = Sequence::length() - 1;

	ss << "BEGIN_SEQUENCE " << Sequence::getName() << "\n";
	ss << "BEGIN_COMMENT\n";
	ss << "CHROMAT_FILE: noChromatFile\n";
	ss << "CALL_METHOD: Bustard\n";
	ss << "QUALITY_LEVELS: 50\n";
	ss << "TIME: " << now << "\n";
	ss << "TRACE_ARRAY_MIN_INDEX: 0\n";
	ss << "TRACE_ARRAY_MAX_INDEX:" << lind << "\n";
	ss << "CHEM: Solexa\n";
	ss << "END_COMMENT\n";
	ss << "BEGIN_DNA\n";
	for (int i = 0; i < Sequence::length(); i++) {
		char q = Sequence::getBaseQuality(i);
		ss << Sequence::getBase(i) << " " << int(q) - fastq_pad << " "
				<< i + 1 <<  "\n";
	}
	ss << "END_DNA\n";
	ss << "END_SEQUENCE\n";

	return ss.str();
}

string Read::getRD(){
	stringstream ss;
	ss << "RD " << Sequence::getName() << " " << Sequence::length() << " 0 0\n\n";
	ss << "QA 1 " <<  Sequence::length() << " 1 " << Sequence::length() << "\n";
	ss << "CHROMAT_FILE: noChromatFile PHD_FILE: " << Sequence::getName() <<
			".phd TIME: " << now << " CHEM: Solexa" << "\n";
	return ss.str();
}

string Read::getAF() {
	stringstream ss;
	ss << "AF " << velID << " ";
	return ss.str();
}
