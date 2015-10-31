/**
 *  Copyright 2015 by Renato Monaro and Silvio Giuseppe
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

namespace oemtp{
	class Circuit{
		public:
    		Circuit();
    		~Circuit();
    		bool Join(Component*);
    		bool Join(Block *Bl);
  			bool Assembly();
  			bool Mapping();
  			bool Compute_G();
  			bool Compute_Ih(bool e);
  			bool Compute_V();  			
  			bool Compute_I();
  			bool Reset();
  			int b;
    	protected:
    	    vector<Component*> vComponent;
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
		};
	}
#endif
