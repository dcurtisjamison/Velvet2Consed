/***********************************************************************
 * FastaHandler.h
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

#ifndef FASTAHANDLER_H_
#define FASTAHANDLER_H_

#include <string>
#include "Read.h"
#include <fstream>
#include <iostream>
using namespace std;

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>
namespace BFS = boost::filesystem;

class FastaHandler {
private:
	BFS::path file;
	bool fasta;
	fstream input;
	iostream::pos_type readHere;

	void setType();
public:
	FastaHandler(string inFile);
	FastaHandler(BFS::path inFile);
	virtual ~FastaHandler();

	string fileName() {return file.string();}
	vector<Read> getReads();
	Read getRead();
	bool hasMoreSeqs();
};

#endif /* FASTAHANDLER_H_ */
