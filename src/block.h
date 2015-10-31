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

#ifndef BLOCK_H
#define BLOCK_H

#include "components.h"

using namespace std;

namespace oemtp{
	class Block{
		public:
			Block();
			~Block();
			unsigned Get_N_Components();
			Component* Get_Component(unsigned N);
			string Get_Block_Name();
			void Set_Block_Name(string N);
		protected:
			vector<Component*> vComponents;
			string Block_Name;
			};
	class Converter: public Block{
		public:
			Converter(string Name,string Na,string Nb,string Nc,string Np,string Nn, double L, double R,double C, double dT);
			Converter(string Name,string Na,string Nb,string Nc,string Np,string Nn);
			Switch *S1,*S2,*S3,*S4,*S5,*S6;
			Diode *D1,*D2,*D3,*D4,*D5,*D6;
			Inductor *La,*Lb,*Lc;
			Resistor *Ra,*Rb,*Rc;
			Capacitor *Cp,*Cn;
		};	
	}
#endif
