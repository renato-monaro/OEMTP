
#include "analysis.h"

using namespace oemtp;
int main(int argc, char* argv[]){
	Circuit Circ1;
	if(argc == 2) { 
		if(Circ1.Interpreter(argv[1]))
			return 1;
		else {
			Circ1.Assembly();
			Circ1.Generate_Output_File();
		}
	}
	else {
		cerr << "Input File not found" << endl;
		return 1;
	}
	return 0;
}
