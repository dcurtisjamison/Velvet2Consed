/***********************************************************************
 * VelvetHandler.h
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

#ifndef VELVETHANDLER_H_
#define VELVETHANDLER_H_

#include <string>
#include <sstream>
#include <iostream>
using namespace std;

#include "Read.h"
#include "FastaHandler.h"
#include "AMOSHandler.h"
#include "DatastoreHandler.h"

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>
namespace BFS = boost::filesystem;

class VelvetHandler {
private:
	template <class T> inline std::string to_string (const T& t) {
		std::stringstream ss;
		ss << t;
		return ss.str();
	}

	BFS::path vDir;
	vector<string> inFiles;
	string timePoint();

public:
	VelvetHandler(string dir);
	virtual ~VelvetHandler();

	bool getFastaNames();
	void getReads(string faDir, DatastoreHandler *dh);
	void getContigs(DatastoreHandler *dh, DatastoreHandler *dr);
};

#endif /* VELVETHANDLER_H_ */
