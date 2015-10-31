#include "block.h"

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

using namespace std;
using namespace oemtp;

Block::Block(){
	};

Block::~Block(){
//	for(unsigned k=0;k<vComponents.size();k++)
//		delete vComponents[0];
	};
	
unsigned Block::Get_N_Components(){
	return vComponents.size();
	}
Component* Block::Get_Component(unsigned N){
	return vComponents[N];
	}
string Block::Get_Block_Name(){
	return Block_Name;
	}
void Block::Set_Block_Name(string N){
	Block_Name=N;
	}
	
Converter::Converter(string Name,string Na,string Nb,string Nc,string Np,string Nn, double L, double R, double C, double dt){
	Block_Name=Name;
	S1=new Switch(Np,Na+Block_Name+"002");
	S2=new Switch(Na+Block_Name+"002",Nn);
	S3=new Switch(Np,Nb+Block_Name+"002");
	S4=new Switch(Nb+Block_Name+"002",Nn);
	S5=new Switch(Np,Nc+Block_Name+"002");
	S6=new Switch(Nc+Block_Name+"002",Nn);
	
	D1=new Diode(Na+Block_Name+"002",Np);
	D2=new Diode(Nn,Na+Block_Name+"002");
	D3=new Diode(Nb+Block_Name+"002",Np);
	D4=new Diode(Nn,Nb+Block_Name+"002");
	D5=new Diode(Nc+Block_Name+"002",Np);
	D6=new Diode(Nn,Nc+Block_Name+"002");
	
	La=new Inductor(Na,Na+Block_Name+"001",L,dt);
	Lb=new Inductor(Nb,Nb+Block_Name+"001",L,dt);
	Lc=new Inductor(Nc,Nc+Block_Name+"001",L,dt);
	
	Ra=new Resistor(Na+Block_Name+"001",Na+Block_Name+"002",R);
	Rb=new Resistor(Nb+Block_Name+"001",Nb+Block_Name+"002",R);
	Rc=new Resistor(Nc+Block_Name+"001",Nc+Block_Name+"002",R);
	
	Cp=new Capacitor(Np,"TERRA",C, dt);
	Cn=new Capacitor("TERRA",Nn,C, dt);
	
	vComponents.push_back(S1); vComponents.push_back(S2); vComponents.push_back(S3);
	vComponents.push_back(S4); vComponents.push_back(S5); vComponents.push_back(S6);
	vComponents.push_back(D1); vComponents.push_back(D2); vComponents.push_back(D3);
	vComponents.push_back(D4); vComponents.push_back(D5); vComponents.push_back(D6);
	vComponents.push_back(La); vComponents.push_back(Lb); vComponents.push_back(Lc);
	vComponents.push_back(Ra); vComponents.push_back(Rb); vComponents.push_back(Rc);
	vComponents.push_back(Cp); vComponents.push_back(Cn); 
	}
	
Converter::Converter(string Name,string Na,string Nb,string Nc,string Np,string Nn){
	Block_Name=Name;
	S1=new Switch(Np,Na);
	S2=new Switch(Na,Nn);
	S3=new Switch(Np,Nb);
	S4=new Switch(Nb,Nn);
	S5=new Switch(Np,Nc);
	S6=new Switch(Nc,Nn);
	
	D1=new Diode(Na,Np);
	D2=new Diode(Nn,Na);
	D3=new Diode(Nb,Np);
	D4=new Diode(Nn,Nb);
	D5=new Diode(Nc,Np);
	D6=new Diode(Nn,Nc);
	
	vComponents.push_back(S1); vComponents.push_back(S2); vComponents.push_back(S3);
	vComponents.push_back(S4); vComponents.push_back(S5); vComponents.push_back(S6);
	vComponents.push_back(D1); vComponents.push_back(D2); vComponents.push_back(D3);
	vComponents.push_back(D4); vComponents.push_back(D5); vComponents.push_back(D6);
	}
