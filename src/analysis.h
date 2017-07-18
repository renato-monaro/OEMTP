/**
 *  Copyright 2015 by Renato Monaro and Silvio Giuseppe
 *  Copyright 2016 by Renato Monaro, Silvio Giuseppe and Heitor Kenzo Koga
 *  Copyright 2017 by Renato Monaro, Silvio Giuseppe and Heitor Kenzo Koga
 *
 * This file is part of Open Electromagnetic Transient Program - OEMTP.
 * 
 * OEMTP is free software: you can redistribute 
 * it and/or modify it under the terms of the GNU General Public 
 * License as published by the Free Software Foundation, either 
 * version 3 of the License, or (at your option) any later version.
 * 
 * Some open source application is distributed in the hope that it will 
 * be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @license GPL-3.0+ <http://spdx.org/licenses/GPL-3.0+>
 */

#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "interface_orelay.h"
#include "components.h"
#include "block.h"
#include "ode2.h"
#include "regex"
#include <memory>
#include <algorithm>
#include <fstream>
#include <math.h>
#include <stdlib.h>

namespace oemtp{
	class Circuit{
		public:
    		Circuit();
    		~Circuit();
			bool Interpreter(string file);
			bool Generate_Output_File();
    		bool Join(Component*);
    		bool Join(Block *Bl);
  			bool Assembly();
  			bool Mapping();
  			bool Compute_G();
  			bool Compute_Ih(bool e);
  			bool Compute_V();  			
  			bool Compute_I();
  			bool Reset();
			void Print_Log(string s, int status);
  			int b;
			int Get_Pos_Component(string name);
			double dt;
			double et;			
			vector<string> Separate_Space(string s);
			string Delete_Space(string s);
			string Make_Regex(string name, int nodes, int param);
			string Make_Regex2(string name, int param);
			string Replace_String(string s, string toReplace, string replaceWith);
			long Set_Jump_Step(double t);
    	    
    	protected:
			vector<int> vGet_V; 
			vector<int> vGet_I;
			vector<int> vGet_Torque; 
			vector<int> vGet_Speed;
			vector<int> vOpen;
			vector<double> vOpen_Time;
			vector<int> vClose;
			vector<double> vClose_Time;
			vector<int> vChange;
			vector<double> vChange_Value;
			vector<double> vChange_Time;		
			vector<string> vComponentName;
			vector<unique_ptr<Component>> vComponent;
    	    gsl_matrix *Gpr;
			gsl_matrix *A;
			gsl_matrix *G;
		    gsl_matrix *C;
			gsl_vector *V;
			gsl_vector *V_Pri;
			gsl_vector *I_Hist;
			gsl_vector *I_Hist_Pri;
			gsl_matrix_view Guu,Guk,Gku,Gkk;
			gsl_vector_view Ik,Vk,Iu,Vu;
			gsl_permutation *puu;
			vector<string> Alias;
			vector<unsigned> Nodes;
			unsigned TNodes;
			unsigned Total_Branches;
			unsigned Total_Nodes;
			string OutFile;
			string LogFile;
			int Sampling;
		};
	}
#endif
