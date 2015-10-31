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

#ifdef OPENRELAY
#include "interface_orelay.h"
using namespace orelay;
OEMTPAcquisition::OEMTPAcquisition(double dT){
	T=0;
	dTns=1E9*dT;
	}

void OEMTPAcquisition::Join(Channel<digital>* Ch,oemtp::Switch* Sw,unsigned Pos, bool D){
	SwitchL.push_back(Sw);
	DigitalSw.push_back(Ch);
	AnalogicSw.push_back(NULL);
	TypeSw.push_back(DIGITAL);
	DirectionSw.push_back(D);
	Position.push_back(Pos);
	}
	
void OEMTPAcquisition::Join(Channel<digital>* Ch,oemtp::Switch* Sw,bool D){
	SwitchL.push_back(Sw);
	DigitalSw.push_back(Ch);
	AnalogicSw.push_back(NULL);
	TypeSw.push_back(DIGITAL);
	DirectionSw.push_back(D);
	Position.push_back(0);
	}
	
void OEMTPAcquisition::Join(Channel<digital>* Ch,oemtp::Switch* Sw){
	SwitchL.push_back(Sw);
	DigitalSw.push_back(Ch);
	AnalogicSw.push_back(NULL);
	TypeSw.push_back(DIGITAL);
	DirectionSw.push_back(OUTPUT);
	Position.push_back(0);
	}

void OEMTPAcquisition::Join(Channel<analogic>* Ch,oemtp::Component *Res,unsigned T, unsigned Pos){
	ComponentL.push_back(Res);
	DigitalRes.push_back(NULL);
	AnalogicRes.push_back(Ch);
	TypeRes.push_back(T);
	DirectionRes.push_back(INPUT);
	Position.push_back(Pos);
	}
	
void OEMTPAcquisition::Join(Channel<analogic>* Ch,oemtp::Component *Res,unsigned T){
	ComponentL.push_back(Res);
	DigitalRes.push_back(NULL);
	AnalogicRes.push_back(Ch);
	TypeRes.push_back(T);
	DirectionRes.push_back(INPUT);
	Position.push_back(0);
	}
	
void OEMTPAcquisition::Join(Channel<analogic> *Ch,oemtp::Ammeter *Res){
	ComponentL.push_back(Res);
	DigitalRes.push_back(NULL);
	AnalogicRes.push_back(Ch);
	TypeRes.push_back(TYPE_CURRENT);
	DirectionRes.push_back(INPUT);
	Position.push_back(0);
	}
	
void OEMTPAcquisition::Join(Channel<analogic> *Ch,oemtp::Voltmeter *Res){
	ComponentL.push_back(Res);
	DigitalRes.push_back(NULL);
	AnalogicRes.push_back(Ch);
	TypeRes.push_back(TYPE_VOLTAGE);
	DirectionRes.push_back(INPUT);
	Position.push_back(0);
	}


bool OEMTPAcquisition::Prepare(float){
	return true;
	}

void OEMTPAcquisition::RefreshTime(){
	T+=dTns;
	MasterClock->insert_Value(T);
	}

bool OEMTPAcquisition::Run(){
	for(unsigned k=0;k<ComponentL.size();k++){
		if(DirectionRes[k]==INPUT){
			switch(TypeRes[k]){
				case TYPE_CURRENT:
					Value=ComponentL[k]->Get_I(Position[k]);
					break;	
				case TYPE_VOLTAGE:
					Value=ComponentL[k]->Get_V(Position[k]);
					break;
//				case TYPE_VALUE:
//					Value=ComponentL[k]->Get_Value(Position[k]);
//					break;
				default:
					Value=ComponentL[k]->Get_I(Position[k]);
				}
			AnalogicRes[k]->insert_Value(Value);
			//DigitalRes[k]->insert_Value(true);
			}
		}
	for(unsigned k=0;k<SwitchL.size();k++){
		if(DirectionSw[k]==INPUT){
			switch(TypeSw[k]){
				case TYPE_CURRENT:
					Value=SwitchL[k]->Get_I(Position[k]);
					break;	
				case TYPE_VOLTAGE:
					Value=SwitchL[k]->Get_V(Position[k]);
					break;
//				case TYPE_VALUE:
//					Value=SwitchL[k]->Get_Value(Position[k]);
//					break;
				default:
					Value=SwitchL[k]->Get_I(Position[k]);
				}
			//AnalogicSw[k]->insert_Value(Value);
			DigitalSw[k]->insert_Value(SwitchL[k]->Get_Status());
			}
		else{
			SwitchL[k]->Set_Status(DigitalSw[k]->get_Value());
			}
		}
	return true;
	}
#endif





