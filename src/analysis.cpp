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

#include "analysis.h"

using namespace std;
using namespace oemtp;

Circuit::Circuit(){
	Total_Branches=0;
	Total_Nodes=0;
	};

Circuit::~Circuit(){
	gsl_matrix_free(Gpr);
	gsl_matrix_free(A);
	gsl_matrix_free(G);
	gsl_vector_free(V);
	gsl_vector_free(V_Pri);
	gsl_vector_free(I_Hist);
	gsl_vector_free(I_Hist_Pri);
	gsl_matrix_free(C);
	};
	
bool Circuit::Join(Component *Comp){
	vComponent.push_back(Comp);
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
		vComponent.push_back(Comp);
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
	for(unsigned k=0;k<vComponent.size();k++)
		vComponent[k]->Reset();
	return true;
	}
	
	


