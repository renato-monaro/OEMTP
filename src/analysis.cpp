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

#include "analysis.h"

using namespace std;
using namespace oemtp;

Circuit::Circuit(){
	Total_Branches=0;
	Total_Nodes=0;
	Sampling = 1;
	};

Circuit::~Circuit(){
	if(Gpr!=NULL){
		gsl_matrix_free(Gpr);
		gsl_matrix_free(A);
		gsl_matrix_free(G);
		gsl_vector_free(V);
		gsl_vector_free(V_Pri);
		gsl_vector_free(I_Hist);
		gsl_vector_free(I_Hist_Pri);
		gsl_matrix_free(C);
		}
	};
	
bool Circuit::Join(Component *Comp){
	vComponent.push_back(unique_ptr<Component>(Comp));
	for(unsigned k=0;k<Comp->Get_N_Nodes();k++){
		Alias.push_back(Comp->Get_Alias(k));
		Nodes.push_back(0);
		}
	Total_Branches+=Comp->Get_N_Branches();
	return true;
	}
	
bool Circuit::Join(Block *Bl){
	Component *Comp;
	for(unsigned b=0;b<Bl->Get_N_Components();b++){
		Comp=Bl->Get_Component(b);
		vComponent.push_back(unique_ptr<Component>(Comp));
		vComponentName.emplace_back(Bl->Get_Block_Name()+"."+Bl->Get_Component_Block_Name(b));
		for(unsigned k=0;k<Comp->Get_N_Nodes();k++){
			Alias.push_back(Comp->Get_Alias(k));
			Nodes.push_back(0);
			}
		Total_Branches+=Comp->Get_N_Branches();
		}
	return true;
	}
	
bool Circuit::Mapping(){
     for(unsigned k=0;k<Alias.size();k++){
        bool Equal=false;
        for(unsigned j=0;j<k;j++){
            if(Alias[k]==Alias[j]){
                 Nodes[k]=Nodes[j];
                 Equal=true;
                }
            }
            if(Equal==false){
                Nodes[k]=Total_Nodes;
                Total_Nodes++;
            }
        }
	unsigned NodeT,ChNode=Total_Nodes-1;
	TNodes=1;
	for(unsigned k=0;k<Alias.size();k++){
		if(Alias[k]==string("TERRA")){
			//kNodes.push_back(Node[k]);
			NodeT=Nodes[k];
			for(unsigned l=k;l<Nodes.size();l++){
				if(Nodes[l]==NodeT){
					Nodes[l]=ChNode;
					}
				else{
					if(Nodes[l]==Nodes[k]){
						Nodes[l]=NodeT;
						}
					}
				}
			}
		}
	#ifdef DEBUG
	for(unsigned k=0;k<Alias.size();k++){
		cout<<"A:"<<Alias[k]<<" N:"<<Nodes[k]<<endl;
		}
	#endif
	return true;
	}
	
bool Circuit::Assembly(){
	Mapping();
	
	Gpr=gsl_matrix_alloc(Total_Branches,Total_Branches);
	A=gsl_matrix_alloc(Total_Branches,Total_Nodes);
	C=gsl_matrix_alloc(A->size2,A->size1);
	G=gsl_matrix_alloc(Total_Nodes,Total_Nodes);
	V=gsl_vector_alloc(Total_Nodes);
	V_Pri=gsl_vector_alloc(Total_Branches);
	I_Hist=gsl_vector_alloc(Total_Nodes);
	I_Hist_Pri=gsl_vector_alloc(Total_Branches);
	
		
	unsigned T_Br=0,T_Nd=0;
	for(unsigned k=0;k<vComponent.size();k++){
		unsigned Br=vComponent[k]->Get_N_Branches();
		unsigned Nd=vComponent[k]->Get_N_Nodes();
		vComponent[k]->Set_Views(
			gsl_matrix_submatrix(Gpr, T_Br, T_Br, Br,Br),
		   	gsl_vector_subvector(V_Pri,T_Br,Br),
   		   	gsl_vector_subvector(I_Hist_Pri,T_Br,Br));
   		  
   		   	vComponent[k]->Compute_Gpr();
   		   	for(unsigned b=T_Br;b<(Br+T_Br);b++){
	   			gsl_matrix_set(A,b,Nodes[2*b],1);
	   			gsl_matrix_set(A,b,Nodes[2*b+1],-1);
	   			}
	   		T_Br+=Br;
   		   	T_Nd+=Nd;
	   		}
   		   	//Definir 1 e -1 para A  
   		#ifdef DEBUG	
	   	cout<<"----A-----"<<endl;
  		for(unsigned b=0;b<Total_Branches;b++){
  		  	for(unsigned n=0;n<Total_Nodes;n++){
				cout<<gsl_matrix_get(A,b,n)<<" ";
				}
				cout<<endl;
			}
	   	cout<<"---------"<<endl;
	   	#endif
	   	
		Guu=gsl_matrix_submatrix (G, 0, 0, G->size1-TNodes, G->size1-TNodes);
		Gkk=gsl_matrix_submatrix (G, G->size1-TNodes, G->size1-TNodes, TNodes, TNodes);
   	   
		Guk=gsl_matrix_submatrix (G, 0,G->size1-TNodes, G->size1-TNodes, TNodes);
		Gku=gsl_matrix_submatrix (G, G->size2-TNodes,0,TNodes,G->size1-TNodes);
   	   
   	   
		puu=gsl_permutation_alloc(Guu.matrix.size1);
   	   
	   	Vu=gsl_vector_subvector(V,0,G->size1-TNodes);
   	   	Vk=gsl_vector_subvector(V,G->size1-TNodes,TNodes);
   	   	Iu=gsl_vector_subvector(I_Hist,0,G->size1-TNodes);
   	   	Ik=gsl_vector_subvector(I_Hist,G->size1-TNodes,TNodes);
	   	
	   	
		
  	
	   		return true;
		}
		
bool Circuit::Compute_G(){
	bool Modified=false;
	for(unsigned k=0;k<vComponent.size();k++){
		Modified|=vComponent[k]->G_Changed();
		}
	if(!Modified)
		return false; //G not Modified
	for(unsigned k=0;k<vComponent.size();k++){
		vComponent[k]->Compute_Gpr();
		}	
	#ifdef DEBUG
	cout<<"----GPR----"<<endl;
	for(unsigned k=0;k<Gpr->size1;k++){
        for(unsigned j=0;j<Gpr->size2;j++){	
        	cout<<gsl_matrix_get(Gpr,k,j)<<" ";
	    }
	    cout<<endl;
	    }
	    	cout<<"-----------"<<endl;
		#endif

	gsl_matrix_set_zero(C);

	gsl_matrix_set_zero(G);
	
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, A, Gpr,0.0,C); //C=A^T*Gpr
		#ifdef DEBUG
		cout<<"----C----"<<endl;
	for(unsigned k=0;k<C->size1;k++){
        for(unsigned j=0;j<C->size2;j++){	
        	cout<<gsl_matrix_get(C,k,j)<<" ";
	    }
	    cout<<endl;
	    }
	    	cout<<"-----------"<<endl;
	   #endif 	
	gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, C,A,0.0,G); //G=C*A
	
		#ifdef DEBUG
	cout<<"----G----"<<endl;
	for(unsigned k=0;k<G->size1;k++){
        for(unsigned j=0;j<G->size2;j++){	
        	cout<<gsl_matrix_get(G,k,j)<<" ";
	    }
	    cout<<endl;
	    }
	    	cout<<"-----------"<<endl;
	   #endif
	
	int s;
		#ifdef DEBUG
    for(unsigned k=0;k<Guu.matrix.size1;k++){
        for(unsigned j=0;j<Guu.matrix.size2;j++){	
        	cout<<gsl_matrix_get(&Guu.matrix,k,j)<<" ";
	    }
	    cout<<endl;
	   }
	  #endif
	gsl_linalg_LU_decomp (&Guu.matrix, puu, &s);
    /*cout<<"LU"	<<endl;
    for(unsigned k=0;k<Guu.matrix.size1;k++){
        for(unsigned j=0;j<Guu.matrix.size2;j++){	
        	cout<<gsl_matrix_get(&Guu.matrix,k,j)<<" ";
	    }
	    cout<<endl;
	   }*/
	return true; //G Modified
	}
	
bool Circuit::Compute_Ih(bool e){
	for(unsigned k=0;k<vComponent.size();k++){
		vComponent[k]->Compute_Ih(e);
		}	
	gsl_blas_dgemv (CblasTrans, 1.0, A, I_Hist_Pri,0.0,I_Hist); //Ihist=A^T*Ihpri
	return true;
	}
	
bool Circuit::Compute_I(){
	for(unsigned k=0;k<vComponent.size();k++){
		vComponent[k]->Compute_I();
		}
	return true;
	}

bool Circuit::Compute_V(){
/*cout<<"Compute Iu' "<<endl;
cout<<"Guk("<<Guk.matrix.size1<<","<<Guk.matrix.size2<<")"<<endl;
cout<<"Vk("<<Vk.vector.size<<")="<<gsl_vector_get(&Vk.vector,0)<<endl;
cout<<"Iu("<<Iu.vector.size<<")"<<endl;*/
	gsl_blas_dgemv (CblasNoTrans, -1.0,&Guk.matrix, &Vk.vector, 1.0, &Iu.vector);//Iu'=-Guk*Vk+Iu
	
/*	for(unsigned j=0;j<Iu.vector.size;j++){	
		cout<<gsl_vector_get(&Iu.vector,j)<<endl;
		}
	cout<<"Compute LU solve'"<<endl;*/
	
	gsl_linalg_LU_solve (&Guu.matrix, puu, &Iu.vector, &Vu.vector);
	//LU Decomp of G
	//Solve V for (G^-1)I
	/*
	cout<<"Compute V"<<endl;
	
	     for(unsigned j=0;j<Vu.vector.size;j++){	
        	cout<<gsl_vector_get(&Vu.vector,j)<<endl;
	    }
	    cout<<"V"<<endl;
	    	     for(unsigned j=0;j<V->size;j++){	
        	cout<<gsl_vector_get(V,j)<<endl;
	    }
	            		    cout<<"Vk"<<endl;
        		     for(unsigned j=0;j<Vk.vector.size;j++){	
        	cout<<gsl_vector_get(&Vk.vector,j)<<endl;
	    }*/
	    
   	gsl_blas_dgemv (CblasNoTrans, -1.0, A, V,0.0,V_Pri); //V_pri=A*V
	return true;
	}

bool Circuit::Reset(){
	for (unsigned k=0;k<vComponent.size();k++)
		vComponent[k]->Reset();
	return true;
	}

int Circuit::Get_Pos_Component(string name) 
{
	std::vector<string>::iterator it;
	it = find (vComponentName.begin(), vComponentName.end(), name);
	if (it != vComponentName.end())
		return it - vComponentName.begin();
	else
	    return -1;
}

vector<string> Circuit::Separate_Space(string s){
   std::string str(s);
   char split_char = ' ';

   std::istringstream split(str);
   std::vector<std::string> separated;
   for (std::string each; std::getline(split, each, split_char); separated.push_back(each));

	for (int i = 0 ; i < separated.size(); i++){
   		if (separated.at(i) == "") {
			separated.erase(separated.begin() + i);
			i--;
		}
 	}
	return separated;
}

string Circuit::Delete_Space(string s){
   s.erase(std::remove(s.begin(), s.end(), ' '), s.end());
   return s;
}

string Circuit::Make_Regex(string name, int nodes, int param){
	return "(?:^\\s*((?:\\S{1,}\\s*){1})\\s*;\\s*("+name+")\\s*;\\s*((?:\\S{1,}\\s*){"+std::to_string(nodes)+"})\\s*;\\s*((?:\\S{1,}\\s*){"+std::to_string(param)+"}))$";
} 

string Circuit::Make_Regex2(string name, int param){
	return "(?:^\\s*(?:"+name+")\\s* \\s*((?:\\S{1,}\\s*){1}\\s*) \\s*((?:\\S{1,}\\s*){"+std::to_string(param)+"}\\s*))$";
} 

string Circuit::Replace_String(string s, string toReplace, string replaceWith){
    return(s.replace(s.find(toReplace), toReplace.length(), replaceWith));
}

long Circuit::Set_Jump_Step(double t){
	if (t - long(t/dt)*dt > 0.0000000001)
		return long(t/dt) + 1;
	else 
		return long(t/dt);
}

bool Circuit::Generate_Output_File(){
	double eps = 0.000000000001;
	double last_time = -1;
	double min = et + 1;
	double min2 = et + 2;
	double p = 0.1;
	long jump_step = 0;
	bool OFile = 1;
	fstream output;
	std::fstream OEMTPlog;
	if (OutFile != "terminal"){
		output.open(OutFile.c_str(), std::ofstream::out | std::fstream::trunc);
		if (output.fail()){
			cerr << " File error " << endl;
			return 1;
		}
	}
	else 
		OFile = 0;

	OEMTPlog.open (LogFile.c_str(), std::fstream::out | std::fstream::app);
	if (OEMTPlog.fail() ){
		cerr << " File error " << endl;
		return 1;
	}

	else {
		OEMTPlog << endl << "@Simulation log: " << endl << endl;
		if (OFile)
			output << "#Time ";
		else
			cout << "#Time ";

		for (int i = 0; i < vGet_V.size(); i++){
			if (OFile)
				output << "V(" << vComponentName.at(vGet_V.at(i))  <<") ";
			else
				cout << "V(" << vComponentName.at(vGet_V.at(i))  <<") ";
		}
		for (int i = 0; i < vGet_I.size(); i++){
			if (OFile)
				output << "I(" << vComponentName.at(vGet_I.at(i))  <<") ";
			else
				cout << "I(" << vComponentName.at(vGet_I.at(i))  <<") ";
		}
		for (int i = 0; i < vGet_Torque.size(); i++){
			if (OFile)
				output << "Torque(" << vComponentName.at(vGet_Torque.at(i))  <<") ";
			else
				cout << "Torque(" << vComponentName.at(vGet_Torque.at(i))  <<") ";
		}
		for (int i = 0; i < vGet_Speed.size(); i++){
			if (OFile)
				output << "Speed(" << vComponentName.at(vGet_Speed.at(i))  <<") ";
			else
				cout << "Speed(" << vComponentName.at(vGet_Speed.at(i))  <<") ";
		}

		for (long k = 0; k < (long)(et/dt); k++){
			if (jump_step == 0) {
				for (int i = 0; i < vOpen_Time.size(); i++){
					if (fabs(k*dt - vOpen_Time.at(i)) < eps || (last_time+eps < vOpen_Time.at(i) && k*dt-eps > vOpen_Time.at(i) )){
						if (Switch* s = dynamic_cast<Switch*>(vComponent.at(vOpen.at(i)).get())){
							s->Open();		
							min = et + 1;
							OEMTPlog << vComponentName.at(vOpen.at(i)) << " opened at " << k*dt << endl; 
						}
					}
					if (min > vOpen_Time.at(i) && dt*k < vOpen_Time.at(i)  && vOpen_Time.at(i) > 0)
						min = vOpen_Time.at(i);
					if (min2 > vOpen_Time.at(i) &&  dt*k < vOpen_Time.at(i) && vOpen_Time.at(i) > 0 && min2 > min)
						min2 = vOpen_Time.at(i);
				}
				for (int i = 0; i < vClose_Time.size(); i++){
					if (fabs(k*dt - vClose_Time.at(i)) < eps || (last_time+eps < vClose_Time.at(i) && k*dt-eps > vClose_Time.at(i) )){
						if (Switch* s = dynamic_cast<Switch*>(vComponent.at(vClose.at(i)).get())){	
							s->Close();		
							min = et + 1;
							OEMTPlog << vComponentName.at(vClose.at(i)) << " closed at " << k*dt << endl; 
						}
					}
					if (min > vClose_Time.at(i) && dt*k < vClose_Time.at(i) && vClose_Time.at(i) > 0)
						min = vClose_Time.at(i);
					if (min2 > vClose_Time.at(i) && dt*k < vClose_Time.at(i) && vClose_Time.at(i) > 0 && min2 > min)
						min2 = vClose_Time.at(i);
				}
				for (int i = 0; i < vChange_Time.size(); i++){
					if (fabs(k*dt - vChange_Time.at(i)) < eps || (last_time+eps < vChange_Time.at(i) && k*dt-eps > vChange_Time.at(i) )){
						if (Resistor* r = dynamic_cast<Resistor*>(vComponent.at(vChange.at(i)).get())) {
							r->Set_Value(vChange_Value.at(i));
							min = et + 1;
							OEMTPlog << vComponentName.at(vChange.at(i)) << " changed at " << k*dt << endl; 
						}
						else if (Capacitor* c = dynamic_cast<Capacitor*>(vComponent.at(vChange.at(i)).get())){
							c->Set_Value(vChange_Value.at(i));		
							min = et + 1;
							OEMTPlog << vComponentName.at(vChange.at(i)) << " changed at " << k*dt << endl; 
						}					
						else if (Inductor* l = dynamic_cast<Inductor*>(vComponent.at(vChange.at(i)).get())){
							l->Set_Value(vChange_Value.at(i));
							min = et + 1;
							OEMTPlog << vComponentName.at(vChange.at(i)) << " changed at " << k*dt << endl; 
						}
						else if (InductionMachine* im = dynamic_cast<InductionMachine*>(vComponent.at(vChange.at(i)).get())){
							im->Set_Mec_Torque(vChange_Value.at(i));
							min = et + 1;
							OEMTPlog << vComponentName.at(vChange.at(i)) << " mechanical torque changed at " << k*dt << endl; 
						}
						else if (AC_Source* ac = dynamic_cast<AC_Source*>(vComponent.at(vChange.at(i)).get())){
							ac-> Set_Voltage(vChange_Value.at(i));
							ac-> Set_Frequency(vChange_Value.at(i + 1)); 
							ac-> Set_Angle(vChange_Value.at(i + 2));
							min = et + 1;
							OEMTPlog << vComponentName.at(vChange.at(i)) << " changed at " << k*dt << endl; 
						}
						//else if (DC_Source* dc = dynamic_cast<DC_Source*>(vComponent.at(vChange.at(i)).get()))
						//else if (Current_Source* cs = dynamic_cast<Current_Source*>(vComponent.at(vChange.at(i)).get()))	
					}
					if (min > vChange_Time.at(i) &&  dt*k < vChange_Time.at(i) && vChange_Time.at(i) > 0)
						min = vChange_Time.at(i);
					if (min2 > vChange_Time.at(i) && dt*k < vChange_Time.at(i) && vChange_Time.at(i) > 0 && min2 > min)
						min2 = vChange_Time.at(i);
				}
				if (min2 < min)
					min = min2;

				jump_step = Set_Jump_Step(min - k*dt);
			}
			if (Compute_G()){
				Compute_Ih(true);
				Compute_V();
				Compute_I();
				Compute_Ih(true);
				Compute_V();
				Compute_I();	
			}
			else {
				Compute_Ih(false);
				Compute_V();
				Compute_I();
			}
			if (!(k%Sampling)){
				if (OFile) 
					output << endl << k*dt << " ";
				else
					cout << endl << k*dt << " ";

				if (k*dt > et*p) {
				cout << "Loading... " << p*100 << "%" << endl;
				p+= 0.1;
				}

				for (int i = 0; i < vGet_V.size(); i++){
					if (OFile)
						output << vComponent[vGet_V.at(i)]->Get_V(0)  <<" ";
					else
						cout << vComponent[vGet_V.at(i)]->Get_V(0)  <<" ";
				}
					
				for (int i = 0; i < vGet_I.size(); i++){
					if (OFile)
						output << vComponent[vGet_I.at(i)]->Get_I(0)  <<" ";
					else
						cout << vComponent[vGet_I.at(i)]->Get_I(0)  <<" ";
				}

				for (int i = 0; i < vGet_Torque.size(); i++){
					if (InductionMachine* im = dynamic_cast<InductionMachine*>(vComponent.at(vGet_Torque.at(i)).get())) {
						if (OFile)
							output << im->Get_Torque()  <<" ";
						else
							cout << im->Get_Torque()  <<" ";
			
					}
				}

				for (int i = 0; i < vGet_Speed.size(); i++){
					if (InductionMachine* im = dynamic_cast<InductionMachine*>(vComponent.at(vGet_Speed.at(i)).get())) {
						if (OFile)
							output << im->Get_Speed()  <<" ";
						else
							cout << im->Get_Speed()  <<" ";
					}
				}
			}
			last_time = k*dt;
			min2 = et + 2;
			if (jump_step != 0)
				jump_step--;
			if (jump_step < 0)
				jump_step = 0;
			if (min < dt*k)
				min = et + 1;
		}	
	}
	output.close();
	OEMTPlog.close();
	return 0;
}

void Circuit::Print_Log (string s, int status) 
{
	std::fstream OEMTPlog;
	OEMTPlog.open (LogFile.c_str(), std::fstream::out | std::fstream::app);
	if (OEMTPlog.fail() ){
		cerr << LogFile << " not found." << endl;
		return;
	}
	else {
		switch (status){
			case 0: 
				OEMTPlog << "' " << s << " ' --- " << " ERROR: Invalid syntax or Invalid component's parameters/nodes" << endl; 
				break;
			case 1:
				OEMTPlog << "' " << s << " ' --- " << "Success" << endl;  
				break;
			case 2:
				OEMTPlog << " @Simulation:" << endl; 
				break;
			case 3:
				OEMTPlog << " @Circuit:" << endl; 
				break;
			case 4:
				OEMTPlog << " @Measure:" << endl; 
				break;
			case 5: 
				OEMTPlog << " @Action:" << endl; 
				break;
			case 6:
				OEMTPlog << " ERROR: header 'Simulation' not found" << endl; 
				break;
			case 7:
				OEMTPlog << " ERROR: header 'Circuit' not found" << endl; 
				break;
			case 8:
				OEMTPlog << " ERROR: header 'Measure' not found" << endl; 
				break;
			case 9:
				OEMTPlog << "' " << s << " ' --- " << " ERROR: Invalid syntax or Component not found" << endl; 
				break;
			case 10:
				OEMTPlog << "' " << s << " ' --- " << " ERROR: Invalid syntax" << endl; 
				break;
			case 11:
				OEMTPlog << " ERROR: Invalid amount of arguments" << endl; 
				break;
			case 12:
				OEMTPlog << "' " << s << " ' --- " << " ERROR: Invalid number" << endl; 
				break;
			case 13:
				OEMTPlog << "' " << s << " ' --- " << " ERROR: The function is not available to this class" << endl; 
				break;
		}	
	}
}

bool Circuit::Interpreter(string file)
{
	bool CircFound = 0;
	bool MeasFound = 0;
	unsigned int status = 1;
	unsigned int step = 0;
	int pos;
	string line;
	std::smatch matchs;	//std::match_results<string> matchs;
	vector<string> blocks;
	vector<string> nodes;
	vector<string> parameters;
	std::regex HeadSimulation("(?:^\\s*@\\s*(simulation)\\s*)$", std::regex_constants::icase);
	std::regex HeadCircuit("(?:^\\s*@\\s*(circuit)\\s*)$", std::regex_constants::icase);
	std::regex HeadMeasure("(?:^\\s*@\\s*(measure|measures)\\s*)$", std::regex_constants::icase);
	std::regex HeadAction("(?:^\\s*@\\s*(action|actions)\\s*)$", std::regex_constants::icase);
	ifstream ComponentTxt;
	std::fstream OEMTPlog;
	pos = file.find(".txt");
	if (pos != std::string::npos)
		LogFile = Replace_String(file, ".txt", "_log.txt");
	else
		LogFile = "simulation_log.txt";
	OEMTPlog.open(LogFile, std::ofstream::out | std::fstream::trunc); // erase file
	ComponentTxt.open(file.c_str());
	OEMTPlog.close(); 
	OutFile = "terminal";
	if (ComponentTxt.fail() ){
		cerr << file << " not found." << endl;
		return 1;
	}

	else{ 
	    while (ComponentTxt){ 

			while (ComponentTxt && step == 0){
				getline(ComponentTxt, line);
				if(std::regex_match(line, matchs, std::regex("(?:^\\s*#.*)$")) || std::regex_match(line, matchs, std::regex("(?:^\\s*)$")) || line == ""){}
				else if (std::regex_match(line, matchs, HeadSimulation)){
					Print_Log(line, 2);
					step = 1;
				}
				else {
					Print_Log(line, 6);
					return 1;
				}
			}
			while (ComponentTxt && step == 1){
				getline(ComponentTxt, line);
				if (std::regex_match(line, matchs, std::regex("(?:^\\s*#.*)$")) || std::regex_match(line, matchs, std::regex("(?:^\\s*)$")) || line == ""){}
				else if (std::regex_match(line, matchs, std::regex("(?:^\\s*(?:DT|dt)\\s*;\\s*((?:\\S{1,}){1}\\s*))$"))){
					vector <string> aux = Separate_Space(matchs[1]);
					dt = stod(aux.at(0));
					if (dt > 0)
						Print_Log(line, 1);
					else
						Print_Log(line, 12);
				}
				else if (std::regex_match(line, matchs, std::regex("(?:^\\s*(?:ET|et)\\s*;\\s*((?:\\S{1,}){1}\\s*))$"))){
					vector <string> aux = Separate_Space(matchs[1]);
					et = stod(aux.at(0));
					if (et > 0)
						Print_Log(line, 1);
					else
						Print_Log(line, 12);
				}
				else if (std::regex_match(line, matchs, HeadCircuit)){
					Print_Log(line, 3);
					step = 2;
					CircFound = 1;
				}
				else
					Print_Log(line, 10);
	
			}
			while (ComponentTxt && step == 2){

				getline(ComponentTxt, line);

				if (std::regex_match(line, matchs, std::regex("(?:^\\s*#.*)$")) || std::regex_match(line, matchs, std::regex("(?:^\\s*)$")) || line == ""){}
				else if (std::regex_match(line, matchs, std::regex(Make_Regex("resistor", 2, 1), std::regex_constants::icase))){
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join(new Resistor(nodes.at(0), nodes.at(1), stod(parameters.at(0) )));
					Print_Log(line, 1);
				}

				else if (std::regex_match (line, matchs, std::regex(Make_Regex("capacitor", 2, 1), std::regex_constants::icase))){ 
				// Capacitor(string N1,string N2,double C, double dT);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join( new Capacitor(nodes.at(0), nodes.at(1), stod(parameters.at(0)), dt) ); 
					Print_Log(line, 1);
				}

				else if (std::regex_match (line, matchs, std::regex(Make_Regex("inductor", 2, 1), std::regex_constants::icase) )){
				// Inductor(string N1,string N2,double L, double dT);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join( new Inductor(nodes.at(0), nodes.at(1), stod(parameters.at(0)), dt) ); 
					Print_Log(line, 1);
				}

				else if (std::regex_match (line, matchs, std::regex(Make_Regex("ac_source", 2, 4), std::regex_constants::icase) )){
				// AC_Source(string N1,string N2,double Vol,double F, double A,double Res,double dT);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join(new AC_Source(nodes.at(0), nodes.at(1), stod(parameters.at(0)), stod(parameters.at(1)), stod(parameters.at(2)), stod(parameters.at(3)), dt) ); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("switch", 2, 3), std::regex_constants::icase) )){
				// Switch(string N1,string N2,bool IniStatus,double ResON, double ResOFF);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join( new Switch(nodes.at(0), nodes.at(1), stod(parameters.at(0)), stod(parameters.at(1)), stod(parameters.at(2)) )); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("switch", 2, 1), std::regex_constants::icase) )){
				// Switch(string N1,string N2,bool IniStatus);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join( new Switch(nodes.at(0), nodes.at(1), stod(parameters.at(0)) )); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("switch", 2, 0), std::regex_constants::icase) )){
				// Switch(string N1,string N2);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					Join( new Switch(nodes.at(0), nodes.at(1) )); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("switch2", 2, 3), std::regex_constants::icase) )){
				// Switch2(string N1,string N2,bool IniStatus,double ResON, double ResOFF);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join( new Switch2(nodes.at(0), nodes.at(1), stod(parameters.at(0)), stod(parameters.at(1)), stod(parameters.at(2)) )); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("diode", 2, 0), std::regex_constants::icase) )){
				// Diode(string N1,string N2)
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					Join( new Diode(nodes.at(0), nodes.at(1) )); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("diode", 2, 2), std::regex_constants::icase) )){
				// Diode(string N1,string N2,double ResON, double ResOFF)
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join( new Diode(nodes.at(0), nodes.at(1), stod(parameters.at(0)), stod(parameters.at(1)) )); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("voltmeter", 2, 0), std::regex_constants::icase) )){
				// Voltmeter(string N1,string N2);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					Join( new Voltmeter(nodes.at(0), nodes.at(1) )); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("ammeter", 2, 0), std::regex_constants::icase) )){
				// Ammeter(string N1,string N2);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					Join( new Ammeter(nodes.at(0), nodes.at(1) )); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("current_source", 2, 2), std::regex_constants::icase) )){
				// Current_Source(string N1,string N2,double Ic, double R);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join( new Current_Source(nodes.at(0), nodes.at(1), stod(parameters.at(0)), stod(parameters.at(1)) )); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("dc_source", 2, 2), std::regex_constants::icase) )){
				// DC_Source(string N1,string N2,double Vol,double Res);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join( new DC_Source(nodes.at(0), nodes.at(1), stod(parameters.at(0)), stod(parameters.at(1)) )); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("lossless_line", 2, 3), std::regex_constants::icase) )){
				// Lossless_Line(string N1,string N2,double d, double l, double c, double dT);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join( new Lossless_Line(nodes.at(0), nodes.at(1), stod(parameters.at(0)), stod(parameters.at(1)), stod(parameters.at(2)), dt )); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("lossless_line", 2, 4), std::regex_constants::icase) )){
				// Lossless_Line(string N1,string N2,double d, double xl, double yc, double f, double dT);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join(new Lossless_Line(nodes.at(0), nodes.at(1), stod(parameters.at(0)), stod(parameters.at(1)), stod(parameters.at(2)), stod(parameters.at(3)), dt) ); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("lossless_line", 2, 2), std::regex_constants::icase) )){
				// Lossless_Line(string N1,string N2, double Zc, double Tau, double dT);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join( new Lossless_Line(nodes.at(0), nodes.at(1), stod(parameters.at(0)), stod(parameters.at(1)), dt )); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("line", 2, 4), std::regex_constants::icase) )){
				// Line(string N1,string N2,double d, double l, double c, double r, double dT);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join(new Line(nodes.at(0), nodes.at(1), stod(parameters.at(0)), stod(parameters.at(1)), stod(parameters.at(2)), stod(parameters.at(3)), dt) ); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("line", 2, 3), std::regex_constants::icase) )){
				// Line(string N1,string N2, double Zc, double Tau, double r, double dT);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join( new Line(nodes.at(0), nodes.at(1), stod(parameters.at(0)), stod(parameters.at(1)), stod(parameters.at(2)), dt )); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("line", 2, 5), std::regex_constants::icase) )){
				// Line(string N1,string N2,double d, double xl, double yc, double f, double r, double dT);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join(new Line(nodes.at(0), nodes.at(1), stod(parameters.at(0)), stod(parameters.at(1)), stod(parameters.at(2)), stod(parameters.at(3)), stod(parameters.at(4)), dt) ); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("converter", 5, 3), std::regex_constants::icase) )){
				// Converter(string Name,string Na,string Nb,string Nc,string Np,string Nn, double L, double R,double C, double dT);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join(new Converter( Delete_Space(matchs[1]), nodes.at(0), nodes.at(1), nodes.at(2), nodes.at(3), nodes.at(4), stod(parameters.at(0)), stod(parameters.at(1)), stod(parameters.at(2)), dt) ); 
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex("converter", 5, 0), std::regex_constants::icase) )){
				// Converter(string Name,string Na,string Nb,string Nc,string Np,string Nn);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					Join(new Converter( Delete_Space(matchs[1]), nodes.at(0), nodes.at(1), nodes.at(2), nodes.at(3), nodes.at(4) )); 
					Print_Log(line, 1);
				}		

				else if (std::regex_match (line, matchs, std::regex(Make_Regex("inductionmachine", 3, 9), std::regex_constants::icase) )){
				// InductionMachine(string N1, string N2, string N3, double r1, double r2, double l1, double l2, double lh, double jj, double kd, int p, double tr, double Dt);
					vComponentName.emplace_back(Delete_Space(matchs[1]));
					nodes = Separate_Space(matchs[3]);
					parameters = Separate_Space(matchs[4]);
					Join(new InductionMachine(nodes.at(0), nodes.at(1), nodes.at(2), stod(parameters.at(0)), stod(parameters.at(1)), stod(parameters.at(2)), stod(parameters.at(3)), stod(parameters.at(4)), stod(parameters.at(5)), stod(parameters.at(6)), stod(parameters.at(7)), stod(parameters.at(8)), dt) ); 
					Print_Log(line, 1);
				}	

				else if (std::regex_match(line, matchs, HeadMeasure)){
					Print_Log(line, 4);
					step = 3;
					MeasFound = 1;
				}
				else 
					Print_Log(line, 10);
			}
			while (ComponentTxt && step == 3){
				getline(ComponentTxt, line);
				if (std::regex_match(line, matchs, std::regex("(?:^\\s*#.*)$")) || std::regex_match(line, matchs, std::regex("(?:^\\s*)$")) || line == ""){}
				else if (std::regex_match(line, matchs, std::regex("(?:^\\s*(?:output)\\s*;\\s*((?:\\S{1,}){1}\\s*))$",std::regex_constants::icase))){
					vector <string> aux = Separate_Space(matchs[1]);
					OutFile = aux.at(0).c_str();
					Print_Log(line, 1);
				}
				else if (std::regex_match(line, matchs, std::regex("(?:^\\s*(?:sampling)\\s*;\\s*((?:\\S{1,}){1}\\s*))$",std::regex_constants::icase))){
					vector <string> aux = Separate_Space(matchs[1]);
					Sampling = stod(aux.at(0));
					if (Sampling > 0)
						Print_Log(line, 1);
					else
						Print_Log(line, 12);
				}
				else if (std::regex_match (line, matchs, std::regex("(?:^\\s*(?:V|v)\\s*(?:\\(|\\[)\\s*((?:\\S{1,}){1})\\s*(?:\\)|\\])\\s*)$"))){
					vector <string> aux = Separate_Space(matchs[1]);
					string CompName = aux.at(0).c_str();
					int aux3 = Get_Pos_Component(CompName);
					if (aux3 >= 0) {
						vGet_V.push_back(aux3);
						Print_Log(line, 1);
					}
					else
						Print_Log(line, 9);
					
				}
				else if (std::regex_match (line, matchs, std::regex("(?:^\\s*(?:I|i)\\s*(?:\\(|\\[)\\s*((?:\\S{1,}){1})\\s*(?:\\)|\\])\\s*)$"))){
					vector <string> aux = Separate_Space(matchs[1]);
					string CompName = aux.at(0).c_str();
					int aux2 = Get_Pos_Component(CompName);
					if (aux2 >= 0) {
						vGet_I.push_back(aux2);
						Print_Log(line, 1);
					}
					else
						Print_Log(line, 9);
				}
				else if (std::regex_match (line, matchs, std::regex("(?:^\\s*(?:torque)\\s*(?:\\(|\\[)\\s*((?:\\S{1,}){1})\\s*(?:\\)|\\])\\s*)$",std::regex_constants::icase))){
					vector <string> aux = Separate_Space(matchs[1]);
					string CompName = aux.at(0).c_str();
					int aux2 = Get_Pos_Component(CompName);
					if (aux2 >= 0) {
						vGet_Torque.push_back(aux2);
						Print_Log(line, 1);
					}
					else
						Print_Log(line, 9);
				}

				else if (std::regex_match (line, matchs, std::regex("(?:^\\s*(?:speed)\\s*(?:\\(|\\[)\\s*((?:\\S{1,}){1})\\s*(?:\\)|\\])\\s*)$",std::regex_constants::icase))){
					vector <string> aux = Separate_Space(matchs[1]);
					string CompName = aux.at(0).c_str();
					int aux2 = Get_Pos_Component(CompName);
					if (aux2 >= 0) {
						vGet_Speed.push_back(aux2);
						Print_Log(line, 1);
					}
					else
						Print_Log(line, 9);
				}

				else if (std::regex_match(line, matchs, HeadAction)){
					Print_Log(line, 5);
					step = 4;
				}
				else 
					Print_Log(line, 10);
			}
			while (ComponentTxt && step == 4){	
				getline(ComponentTxt, line);
				if (std::regex_match(line, matchs, std::regex("(?:^\\s*#.*)$")) || std::regex_match(line, matchs, std::regex("(?:^\\s*)$")) || line == ""){}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex2("open", 1), std::regex_constants::icase) )){
					vector <string> aux = Separate_Space(matchs[1]);
					string CompName = aux.at(0).c_str();
					int aux2 = Get_Pos_Component(CompName);
					if (aux2 >= 0) {
						Switch* s = dynamic_cast<Switch*>(vComponent.at(aux2).get());
						Switch2* s1 = dynamic_cast<Switch2*>(vComponent.at(aux2).get());
						if (s != NULL || s1 != NULL)
							vOpen.push_back(aux2);
						else {
							Print_Log(line, 13);
							break;
						}
					}
					else {
						Print_Log(line, 9);	
						break;
					}	
					vOpen_Time.push_back(stod(matchs[2]));
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex2("close", 1), std::regex_constants::icase) )){
					vector <string> aux = Separate_Space(matchs[1]);
					string CompName = aux.at(0).c_str();
					int aux2 = Get_Pos_Component(CompName);
					if (aux2 >= 0) {
						Switch* s = dynamic_cast<Switch*>(vComponent.at(aux2).get());
						Switch2* s1 = dynamic_cast<Switch2*>(vComponent.at(aux2).get());
						if (s != NULL || s1 != NULL)
							vClose.push_back(aux2);
						else {
							Print_Log(line, 13);
							break;
						}
					}
					else {
						Print_Log(line, 9);		
						break;
					}
					vClose_Time.push_back(stod(matchs[2]));
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex2("change", 2), std::regex_constants::icase) )){
					vector <string> aux = Separate_Space(matchs[1]);
					vector <string> aux2 = Separate_Space(matchs[2]);
					string CompName = aux.at(0).c_str();
					int aux3 = Get_Pos_Component(CompName);					
					if (aux3 >= 0) {
						Resistor* r = dynamic_cast<Resistor*>(vComponent.at(aux3).get());
						Capacitor* c = dynamic_cast<Capacitor*>(vComponent.at(aux3).get());
						Inductor* l = dynamic_cast<Inductor*>(vComponent.at(aux3).get());
						Current_Source* cs = dynamic_cast<Current_Source*>(vComponent.at(aux3).get());
						DC_Source* dc = dynamic_cast<DC_Source*>(vComponent.at(aux3).get());
						InductionMachine* im = dynamic_cast<InductionMachine*>(vComponent.at(aux3).get());
						if (r != NULL || c != NULL || l != NULL || cs != NULL || dc != NULL || im != NULL )
							vChange.push_back(aux3);
						else {
							Print_Log(line, 13);
							break;
						}
					}
					else {
						Print_Log(line, 9);
						break;
	                }	
					vChange_Time.push_back(stod(aux2.at(0)));
					vChange_Value.push_back(stod(aux2.at(1)));
					
					Print_Log(line, 1);
				}
				else if (std::regex_match (line, matchs, std::regex(Make_Regex2("change", 4), std::regex_constants::icase) )){
					vector <string> aux = Separate_Space(matchs[1]);
					vector <string> aux2 = Separate_Space(matchs[2]);
					string CompName = aux.at(0).c_str();
					int aux3 = Get_Pos_Component(CompName);					
					if (aux3 >= 0) {
						if (AC_Source* ac = dynamic_cast<AC_Source*>(vComponent.at(aux3).get()) ){
								vChange.push_back(aux3);
								vChange.push_back(-1);
								vChange.push_back(-1);
						}
						else {
							Print_Log(line, 13);
							break;
						}
					}
					else {
						Print_Log(line, 9);
						break;
	                }	
					vChange_Time.push_back(stod(aux2.at(0)));
					vChange_Time.push_back(-2);
					vChange_Time.push_back(-2);
					vChange_Value.push_back(stod(aux2.at(1)));
					vChange_Value.push_back(stod(aux2.at(2)));
					vChange_Value.push_back(stod(aux2.at(3)));
					
					Print_Log(line, 1);
				}
				else 
					Print_Log(line, 10);
			}
		}
	if (!CircFound)
		Print_Log(line, 7);
	if (!MeasFound)
		 Print_Log(line, 8);
	
	} 
	ComponentTxt.close();
	return 0;

}
	
	


