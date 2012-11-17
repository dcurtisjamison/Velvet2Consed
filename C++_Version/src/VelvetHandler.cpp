/***********************************************************************
 * VelvetHandler.cpp
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

#include "VelvetHandler.h"

VelvetHandler::VelvetHandler(string dir)
	: vDir(dir)
{

}

VelvetHandler::~VelvetHandler() {
	// TODO Auto-generated destructor stub
}

/*
 * Extract input files from log file
 */
bool VelvetHandler::getFastaNames(){

	string s;
	boost::regex e("\\s\\S*velveth");
	boost::smatch w;
	BFS::path logFile(vDir / "Log");
	ifstream log (logFile.string().c_str());

	cout << "[" << timePoint() << "] Checking " << logFile.string() << " for " << e.str() << "(velveth command)\n";

	while (getline(log,s)) {
		if (boost::regex_search(s, w, e, boost::match_extra)) {
			// line with "velveth" contains our filenames
			typedef boost::tokenizer<boost::char_separator<char> >
			    tokenizer;
			boost::char_separator<char> sep(" ");
			tokenizer tok(s, sep);
			tokenizer::iterator start = tok.begin();
			// first 3 tokens are command related
			for (int i = 0; i < 4; i++) {
				start++;
			}
			while (start != tok.end()) {
				boost::regex m("-"); // ignore switches
				if (!boost::regex_search(*start, w, m, boost::match_extra)) {
					//filenames; order is important
					inFiles.push_back(*start);
				}
				start++;
			}
			if (inFiles.size() > 0) {
				return true;
			}
		}
	}
	return false;
}

/*
 * get reads from sequence & fastq files
 */
void VelvetHandler::getReads(string faDir, DatastoreHandler *dh) {

	vector<Read> velReads, fullReads;
	vector<string> fastaFiles;

	// first we check to make sure the fastq files are accessible
	BFS::path vSeq (vDir / "Sequences");
	if (! BFS::exists(vSeq)) {
		// TODO add a graceful error handler
		cerr << vSeq << " doesn't exist\n";
		exit(0);
	}
	BOOST_FOREACH(string file, inFiles) {
		// check for a faDir, otherwise take the filename at face value

	if (faDir != "null") {
			file.insert(0,"/");
			file.insert(0,faDir);
		}
		if (! BFS::exists(file)) {
			// TODO add a graceful error handler
			cerr << "FATAL!\n\tFile " << file.c_str() << " doesn't exist. Check your directory settings.\n";
			exit(0);
		} else {
			fastaFiles.push_back(file);
		}
	}

	// files are OK, get the reads
	FastaHandler vs(vSeq);
	Read faRead;
	boost::regex a2n("N");
	boost::regex rID("SEQUENCE_(\\d*)_length");
	boost::match_results<std::string::const_iterator> m;
	int readCount = 0;
	BOOST_FOREACH(string file, fastaFiles) {
		FastaHandler fh (file);
		cout << "[" << timePoint() << "] Fetching reads from sequence file " << fh.fileName() << "\n";
		while (fh.hasMoreSeqs()) {
			faRead = fh.getRead();
			if ( ! vs.hasMoreSeqs()) {
				cerr << "Not enough reads in velvet Sequence file!\n";
			} else {
				Read vRead = vs.getRead();
				if (faRead.getSequence() != vRead.getSequence()) {
					string masked = boost::regex_replace(faRead.getSequence(), a2n, "A");
					if (masked != vRead.getSequence()) {
						cerr << "supposed paired sequences unmatchable\n";
						exit(0);
					}
				}
				boost::regex_search(vRead.getName(), m, rID);
				faRead.setVelvetID(string(m[1])); 			//velvet read numbering starts at 0
				faRead.setAMOSID(atoi(faRead.getVelvetID().c_str()) + 1);	//amos read numbering starts at 1
				readCount++;
				dh->storeRead(faRead);
			}
		}
		cout << "[" << timePoint() << "] Processed " << readCount << " reads\n";
	}
}

void VelvetHandler::getContigs(DatastoreHandler *dh, DatastoreHandler *dr) {

	BFS::path aFile (vDir / "velvet_asm.afg");
	if (! BFS::exists(aFile)) {
		//TODO add a graceful error handler
		cerr << aFile << " doesn't exist!\n";
		exit(0);
	}
	AMOSHandler ah(aFile);
	cout << "[" << timePoint() << "] Fetching contigs from " << ah.fileName() << "\n";
	int ctgCount = 0;
	while (ah.hasMoreSeqs()) {
		cout << "[" << timePoint() << "] Reading contig\n";
		Contig ctg = ah.getContig();
		cout << "[" << timePoint() << "] Recomputing quality\n";
		ctg.computeQuality(dr);
//		ctg.repairEnds(dr,1);
		cout << "[" << timePoint() << "] Fetched a contig\n";
		dh->storeContig(ctg);
		ctgCount++;
	}
	cout << "[" << timePoint() <<"] Processed " << ctgCount << " contigs\n";

}

string VelvetHandler::timePoint() {
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
