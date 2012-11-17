/***********************************************************************
 * AMOSHandler.h
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

#ifndef AMOSHANDLER_H_
#define AMOSHANDLER_H_

#include <string>
#include <fstream>
#include <iostream>
using namespace std;

#include "Contig.h"
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>
namespace BFS = boost::filesystem;

#include <foundation_AMOS.hh>

class AMOSHandler {
private:
	BFS::path file;
	ifstream input;
	iostream::pos_type readHere;

public:
	AMOSHandler();
	AMOSHandler(BFS::path inFile);
	~AMOSHandler();

	string fileName() {return file.string();}
	vector<Contig> getContigs();
	Contig getContig();
	bool hasMoreSeqs();
};

#endif /* AMOSHANDLER_H_ */
