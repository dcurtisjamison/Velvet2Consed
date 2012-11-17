/***********************************************************************
 * AMOSHandler.cpp
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

#include "AMOSHandler.h"

AMOSHandler::AMOSHandler() {
	// TODO Auto-generated constructor stub

}

/*
 * standard constructor
 */
AMOSHandler::AMOSHandler(BFS::path inFile) {
	file = inFile;
	input.open(inFile.string().c_str());
	if (!input) {
		cerr <<"Couldn't create an open file for " << inFile.string() << "\n";
		exit(0);
	}
	readHere = input.tellg();
}

AMOSHandler::~AMOSHandler() {
	// TODO Auto-generated destructor stub
}

vector<Contig> AMOSHandler::getContigs(){
	string s, name, sequence, quality;
	vector<Contig> v;

	// read in the sequences
	// TODO add in some error checking
	while (hasMoreSeqs()) {
		v.push_back(getContig());
	}
	return v;
}

Contig AMOSHandler::getContig(){
	AMOS::Message_t amos_mess;//	AMOS::Universal_t obj;
	AMOS::Contig_t aCon;

	input.seekg(readHere); // make sure we are at the current reading site

	try {
		if (amos_mess.read(input)) {
			while (amos_mess.getMessageCode() != AMOS::Contig_t::NCODE) {
				//skip through RED and other non-contig messages
				amos_mess.read(input);
			}
		} else {
			// This error likely means a corrupt AMOS file
			cerr << "FATAL: AMOS Message not read!\n";
			cerr << "\tinput file position = " << input.tellg() <<"\n";
			string line;
			input >> line;
			cerr << "\tinput line = " << line << "\n";
			exit(0);
		}
	}
	catch (AMOS::Exception_t & e) {
		cerr << "AMOS Error: " << endl << e;
		exit(0);
	}
	readHere = input.tellg();
	aCon.readMessage(amos_mess);
	Contig myCon(aCon);
	return myCon;
}

bool AMOSHandler::hasMoreSeqs() {
	input.seekg(readHere);
	input.peek();
	return (input.good());
}
