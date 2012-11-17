/***********************************************************************
 * FastaHandler.cpp
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

#include "FastaHandler.h"

FastaHandler::FastaHandler(string inFile)
	: file(inFile)
{
	input.open(inFile.c_str(), ios::in);
	readHere = input.tellg();
	setType();
}

FastaHandler::FastaHandler(BFS::path inFile)
	: file (inFile)
{
	input.open(inFile.string().c_str(), ios::in);
	readHere = input.tellg();
	setType();
}

FastaHandler::~FastaHandler() {
	// TODO Auto-generated destructor stub
}

void FastaHandler::setType () {
	char first = input.peek();
	if (first == '>') {
		fasta = true;
	} else {
		if (first != '@') {
			cerr << file.string() << " not a valid fasta or fastq file\n";
			exit(0);
		}
		fasta = false;
	}
}

/*
 * turn a fasta or fastq file into reads
 */
vector<Read> FastaHandler::getReads() {
	string s, name, sequence, quality;
	vector<Read> v;

	// read in the sequences
	// TODO add in some error checking
	while (hasMoreSeqs()) {
		v.push_back(getRead());
	}
	return v;
}

Read FastaHandler::getRead() {
	string s, name, seq, qual;

	input.seekg(readHere);
	getline(input, s);
	name = s.substr(1);
	getline(input, seq);
	if (! fasta) {
		getline(input, s);
		getline(input, qual);
	}
	readHere = input.tellg();
	return Read(name, seq, qual);
}

bool FastaHandler::hasMoreSeqs() {
	input.seekg(readHere);
	input.peek();
	return (! input.eof());
}
