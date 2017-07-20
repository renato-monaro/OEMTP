
#include "analysis.h"

using namespace oemtp;
int main(int argc, char* argv[]){
	Circuit Circ1;
	bool args = (argc == 2);

#ifdef GNUPLOT
	args = (argc == 2 || argc == 3);
#endif

	if(args) { 
		if(Circ1.Interpreter(argv[1]))
			return 1;
		else {
			Circ1.Assembly();
			#ifdef GNUPLOT
			if (argc == 3)
				Circ1.Define_Folder(argv[2]);
			#endif
			Circ1.Generate_Output_File();
				
		}
	}
	else {
		cerr << "Input File not found" << endl;
		return 1;
	}
	return 0;
}
