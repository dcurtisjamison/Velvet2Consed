/*
 VelvetToConsed.cpp
 Author: jam5pi
   Date: Oct 13, 2009

Main wrapper for converting Velvet assemblies to Consed
********************************************************************
 Revision History

********************************************************************
 Copyright 2009, Children's Hospital of Cincinnati

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// header stuff
#include <getopt.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

#include <boost/tokenizer.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
namespace bfs = boost::filesystem;

#include "DatastoreHandler.h"
#include "VelvetHandler.h"
#include "AMOSHandler.h"
#include "Tile.h"
#include "Read.h"
#include "Contig.h"

// a place to put our parameters
struct GlobalArgs_t {
	int toDo;
	string vDir;
	string faDir;
	string dsDir;
	string outname;
	vector<int> contig;
} globalArgs;

static const char *optString = "TCEv:f:d:a:o:";

static const struct option longOpts[] = {
		{"Translate_immediately", no_argument, NULL, 'T'},
		{"Create_data_store", no_argument, NULL, 'C'},
		{"Extract_contig", no_argument, NULL, 'E'},
		{"velvet_dir", required_argument, NULL, 'v'},
		{"fasta_dir", required_argument, NULL, 'f'},
		{"dataStore_dir", required_argument, NULL, 'd'},
		{"export_name", required_argument, NULL, 'o'},
		{"contigID", required_argument, NULL, 'a'},
		{NULL, no_argument, NULL, 0}
};

// Quicka and dirty timestamp for logging purposes
string timePoint() {
	time_t t = time(0);
	tm ts = *localtime(&t);
	stringstream ss;
	ss.fill('0');
	ss.width(2);
	ss << ts.tm_hour << ":";
	ss.fill('0');
	ss.width(2);
	ss << ts.tm_min << ":";
	ss.fill('0');
	ss.width(2);
	ss << ts.tm_sec;
	return ss.str();
}


// This subroutine whinges about how to use the system
void printUsage() {
	cout << "USAGE:\n\tVelvetToConsed -T|C|E -v <velvetDir> -f <fastaDir> [-d <dataStoreDir>] [-a contigName]\n";
	cout << "\tSWITCHES\n";
	cout << "\t\tT\t\tTranslate immediately\n";
	cout << "\t\tC\t\tCreate a data store\n";
	cout << "\t\tE\t\tExtract contig\n";
	cout << "\tPARAMETERS\n";
	cout << "\t\t-v <velvetDir>\tPath to velvet directory (required for modes T&C)\n";
	cout << "\t\t-f <fastaDir>\tPath to fastq/fastq files (required for modes T&C)\n";
	cout << "\t\t-o <outputDir>\tOutput files name (optional for modes T&E)\n";
	cout << "\t\t-d <dataStoreDir>\tDatastore path (optional for all modes)\n";
	cout << "\t\t-a contigName[,contigName] \tSpecific contig list to extract (optional for all modes)\n";
	cout << "\tSee manual for more detail\n";
}

// quick wrapper around the boost library
bool fileTest(string dir) {
	if (!bfs::exists(dir)) {
		return false;
	}
	return true;
}

// Sanity check for conflicting and/or missing parameters
bool checkArgs() {
	bool good = 1;
	switch (globalArgs.toDo) {
	case 1 :
	case 3 :
		good = good && fileTest(globalArgs.faDir);
		good = good && fileTest(globalArgs.vDir);
		// fall through on purpose
	case 5 :
		good = good && fileTest(globalArgs.dsDir);
		break;
	}
	return good;
}

// translate without saving a data store
void translate() {

}

// create a new data store
void createDataStore() {

	cout << "[" << timePoint() << "] Creating a new data store\n\tVelvetDir = " <<
			globalArgs.vDir << "\n\tfastaDir = " << globalArgs.faDir <<
			"\n\tdataStoreDir = " << globalArgs.dsDir << "\n";

	if (globalArgs.dsDir == "undef") {
		globalArgs.dsDir = "/tmp";
	}

	if (!checkArgs()) {
		cout << "Sorry, your arguments are invalid.\n";
		printUsage();
		exit(1);
	}

	//define our file handlers
	VelvetHandler *vh = new VelvetHandler(globalArgs.vDir);
	DatastoreHandler *dh = new DatastoreHandler(globalArgs.dsDir, DatastoreHandler::CREATE);
	DatastoreHandler *dr = new DatastoreHandler(globalArgs.dsDir, DatastoreHandler::FETCH);

	// use the velvet handler to put things into the data stores
	if (vh->getFastaNames()) {
		vh->getReads(globalArgs.faDir, dh);
		vh->getContigs(dh, dr);
	}
}

// extract a contig into an ace format file
void extractContig() {
	if (!checkArgs()) {
		cout << "Sorry, your arguments are invalid.\n";
		printUsage();
		exit(1);
	}

	// make a time stamp
	string rdays[7] = {"Sun","Mon","Tues","Wed","Thurs","Fri","Sat"};
	string rmonths[12] = {"Jan","Feb","Mar","Apr","May","Jun",
			"Jul","Aug","Sep","Oct","Nov","Dec"};
	time_t t = time(0);
	tm ts = *localtime(&t);
	stringstream ss;
	ss << rdays[ts.tm_wday] << " " << rmonths[ts.tm_mon] <<
			" " << timePoint() << " " << ts.tm_year;
	string timestamp = ss.str();

	// create the ace directory structure and open output files
	if (globalArgs.outname.length() < 1) {
		globalArgs.outname = string("VelvetExport");
	}
	bfs::path ace(globalArgs.outname);
	bfs::path phd(globalArgs.outname);
	try {
		if (bfs::create_directory(ace)) {
			ace /= "edit_dir";
			phd /= "phdball_dir";
			if (bfs::create_directory(ace) && bfs::create_directory(phd)) {
				ace /= globalArgs.outname;
				phd /= globalArgs.outname;
				ace.replace_extension("ace");
				phd.replace_extension("phdball");
			}
		}
	}
	catch (const std::exception & ex) {
		cerr << "Couldn't create output directories:\n\t" << ex.what() << "\n";
		cerr << "writing files to /tmp\n";
		ace = bfs::path("/tmp/VelvetExport.ace");
		phd = bfs::path("/tmp/VelvetExport.phdball");
	}
	ofstream ACE(ace.string().c_str());
	ofstream PHD(phd.string().c_str());

	// open up the needed files
	cout << "[" << timePoint() << "] Mapping input files into memory\n";
	DatastoreHandler *dh = new DatastoreHandler(globalArgs.dsDir, DatastoreHandler::FETCH);

	// get the contigs
	vector<int> ctgList;
	int rdCnt = 0;
	if (globalArgs.contig.size() == 0) {
		ctgList = dh->getContigList();
	} else {
		ctgList = globalArgs.contig;
	}
	cout <<  "[" << timePoint() <<
			"] Computing read counts for " << ctgList.size() << " contigs\n";
	BOOST_FOREACH(int i, ctgList) {
		Contig ctg = dh->getContig(i);
		rdCnt += ctg.readCount();
	}
	BOOST_FOREACH(int i, ctgList) {
		cout << "[" << timePoint() << "] Extracting contig " << i << "\n";
		Contig ctg = dh->getContig(i);
		ACE << "AS " << ctgList.size() << " " << rdCnt << "\n\n";
		ACE << ctg.getACE() << "\n\n";
		vector<Read> readList; // hang onto the reads, we need to iterate them again
		BOOST_FOREACH(Tile t, ctg.getTileList()) {
			Read curread = dh->getRead(t.getAMOSID());
			readList.push_back(curread);
			ACE << curread.getAF() << t.getAF() << "\n";
		}
		cout << "[" << timePoint() << "] Writing " <<
				readList.size() << " reads for contig " << i << "\n";
		BOOST_FOREACH(Read r, readList) {
			ACE << r.getRD();
			PHD << r.getPHD();
		}
		ACE << "\nWA{\nphdBall velvet2Consed " << timestamp << "\n../phdball_dir/" <<
				phd.filename() << "\n}\n";
	}
	cout << "[" << timePoint() << "] Finished.\n";
}

vector<int> parseContigArg (string list) {
	vector<int> retVal;
	typedef boost::tokenizer<boost::char_separator<char> >
		tokenizer;
	boost::char_separator<char> sep (",");
	tokenizer tok (string(optarg), sep);
	tokenizer::iterator word = tok.begin();
	while (word != tok.end()) {
		retVal.push_back(atoi(word->c_str()));
		word++;
	}
	return retVal;
}

int main(int argc, char* argv[]) {

	int opt = 0;
	int longIndex;
	globalArgs.toDo = 0;
	globalArgs.faDir = "undef";
	globalArgs.dsDir = "/tmp";
	globalArgs.vDir = "undef";
	extern char *optarg;

	string version("vers 0.2b");
	cout << argv[0] << " " << version << "\n©2009 Children's Hospital of Cincinnati \n\n";

//	extern int optind, optopt, opterr;

	// process command line arguments
	opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
	if (opt == -1) {
		cerr << "getop processing fail\n";
		printUsage();
		exit(0);
	}
	while (opt != -1) {
		switch (opt) {
		case 'T':
			globalArgs.toDo += 1;
			break;
		case 'C':
			globalArgs.toDo += 3;
			break;
		case 'E':
			globalArgs.toDo += 5;
			break;
		case 'v':
			globalArgs.vDir = optarg;
			break;
		case 'f':
			globalArgs.faDir = optarg;
			break;
		case 'd':
			globalArgs.dsDir = optarg;
			break;
		case 'o':
			globalArgs.outname = optarg;
			break;
		case 'a':
			globalArgs.contig = parseContigArg(optarg);
			break;
		case '?':
			cerr << "ignored unknown option " << opt << "\n";
//			printUsage();
//			exit(1);
			break;
		default:
			break;
		}
		opt = getopt(argc, argv, optString);
	}

	// decide what we are going to do today
	switch(globalArgs.toDo) {
	case 1:
		createDataStore();
		extractContig();

		break;
	case 3:
		createDataStore();
		break;
	case 5:
		extractContig();
		break;
	default:
		cout << "invalid command choice\n";
	}
}
