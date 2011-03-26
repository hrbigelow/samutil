#include "misc.h"
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;


string PrintPlusMinusBool(bool b) { return b ? "+" : "-"; }
string PrintYesNoBool(bool b){ return b ? "yes" : "no"; }
string PrintStrandBool(bool b){ return b ? "SAME" : "OPP"; }

bool ValidOutfile(const string outfile){

	ofstream of;
	ostringstream errmsg;
	of.open(outfile.c_str());
	if (! of) {
		cerr<<"Couldn't open "<<outfile<<" for writing.";
		return false;
	}
	else {
		of.close();
		return true;
	}

}
