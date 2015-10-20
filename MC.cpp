/*
* MC.cpp
* Simulation of 2-D Rods By GCMC with SUS, No energy
* Author: Yuding Ai
* Date: 2015.10.19
* *************************** MC implementation ********************************
* This simulation follows Successive Umbrella sampling without energy.
* ******************************************************************************
*
*/


#include "MC.h"

MC::MC(long int ST, int LEN,int C, int R, double Z)
{
	VRodlist.clear(); // the list storage the Vertical Rods;
    HRodlist.clear(); // the list storage the Horizantal Rods;
	r = R;
	c = C;
	length = LEN;
	step = ST; // now the step stands for the steps in each SUS window
	z = Z;
	nh=nv=dh=dv=ah=av=0;
}

vector<HR> MC::getVRodlist() 
{
	return VRodlist;
}

vector<HR> MC::getHRodlist() 
{
	return HRodlist;
}

double MC::getTho() const
{
	double tho;	
	tho = double(length*(av+ah-dv-dh))/double(r*c);
	return tho;
}

double MC::getQ() const
{
	double Q;	
	Q = (nv - nh)/(nv + nh);
	return Q;
}

double MC::getAaccp() const
{
	double A;
	A = z*(double(r*c))/(double(av+ah-dv-dh+1.0)*double(length));
	return A;
}

double MC::getDaccp() const
{
	double D;
	D = (double(av+ah-dv-dh)*double(length))/(z*(double(r*c)));
	return D;
}

double MC::getNh() const
{
	return nh;
}
double MC::getNv() const
{
	return nv;
}


void MC::Add(Cells &s,double &prob,double &proba)
{
	int x,y,o; // pick a random position and orientation for the HR to be added;
	x = rand()%c;
	y = rand()%r;
	o = rand()%1;// change it to 1 for  lattice gas case

	if(s.getSquare(x,y).isEmpty()) // if it's open, try to do Addition;
	{
		HR rod(x,y,length,o);

		//======================== Vertical ===============================

		if(o == 0)
		{
			int counter = 0;

			for (int j = 0; j < length-1; j++)
			{					
				// check if the vertical space is wide open
				if(s.getSquare(x,(y+j+1)%r).isOccupied())
				{
					counter++;
				}					
			}
			if(counter == 0)
			{
				if(prob <= proba)
				{					
					// Do addition;
					// push the new rod into the Rodlist;
					VRodlist.push_back(rod);
					av++;
					nv++;// accumulate the # of ver rod;
					// update new N, E and new config;
					for (int i = 0; i < length; i++)
					{	
						s.getSquare(x,(y+i)%r).setStatus(1);
					}
				}
			}								
		}

		else 
		{
        //======================= Horizontal  ============================
			int counter = 0;
			for (int j = 0; j< length-1 ; j++)
			{
				// check if the horizontal space is wide open
				if(s.getSquare((x+1+j)%c,y).isOccupied())
				{
					counter++;
				}							
			}
			if (counter == 0)
			{
				if(prob <= proba)
				{
					//Do addition;
					//push the new rod into the Rodlist;
					HRodlist.push_back(rod);
					ah++;
					nh++;// accumulate the # of hor rod;

					// update new N, E and new config;
					for (int i = 0; i < length; i++)
					{
						s.getSquare((x+i)%c,y).setStatus(1);
					}
				}	
			}							
		}
    }
}

void MC::Del(Cells &s,double &prob,double &probd,double &size)
{
	//Do Del;
	if(nv+nh > 0)// make sure there are rod;
	{
		int indx; // pick a random index from the Rodlist;
		indx = rand()%int(nv+nh);

		//remove Rodlist[indx];
		int x,y;// the position of the target on the cells;

		if(prob <= probd)
		{
			if (indx < nv) // vertical
			{
				x = VRodlist[indx].getX();
				y = VRodlist[indx].getY();

				// --------------------- it's a vertical rod -----------------------			
				for(int i = 0; i<VRodlist[indx].getLength(); i++)
				{
					// update the new config of cells
					s.getSquare(x,(y+i)%r).setStatus(0);
				}
				// remove the target rod from the vector Rodlist;
				VRodlist.erase(VRodlist.begin() + indx);
				nv--;// substract the # of ver rod;
				dv++;
			}

			else
			{
				x = HRodlist[indx - nv].getX();
				y = HRodlist[indx - nv].getY();
				// --------------------- it's a Horizontal rod -----------------------
				for(int i = 0; i<HRodlist[indx-nv].getLength(); i++)
				{
					// update the new config of cells
					s.getSquare((x+i)%c,y).setStatus(0);
				}
				// remove the target rod from the vector Rodlist;
				HRodlist.erase(HRodlist.begin() + indx - nv);
				nh--;// substract the # of hor rod;
				dh++;				
			}
		}										
	}
}


array<double,10000> MC::MCSUS()
{
	Cells s(c,r);  //  setting the lattice;
	//==========================================================  declare instance variables ============================================================= //
	stringstream sh;
	sh.precision(20);
	double addordel;           // the prob to decide either add or del;
	double probd,proba;      // the acceptance prob of addition and deletion; 
	double prob;               // the prob to decide either accept add/del;
	double aaccp,daccp;      // the acceptance probabilities: 
	double V = double(r*c);    // the total lattice size
	double K = double(length); //
    array<double,10000> WF;

	srand(time(NULL));
	long int i = 0; // counter for each window
	double w = 1.0; // a counter that keep track of the index of window
	double fu,fl; // occurrence counter
	

	//================================Start my MC simulation=================================
	while (w <= V/K) // while loop terminate until finish the last window; window[V/K]
	{
		fu = fl = 0; // initially set to 0
		i = 0;
		while(i < step || nv+nh != w)
		// Simulation for each window
		{
			i++;
			// generate a random probability to decide either add or del;
			addordel = rand()%2;
		    double size = nv+nh;			
			
			prob = ((double) rand() / (RAND_MAX)); 

			aaccp = (z*V)/((size+1.0)*K)*(exp(WF[int(size+1)] - WF[int(size)]));
			daccp = (size*K)/(z*V)*(exp(WF[int(size-1)] - WF[int(size)]));	

			probd = min(1.0,daccp);
			proba = min(1.0,aaccp);

	        // ===========================Addition ===================================
			if(addordel == 0) 
			{
				if (size == w - 1.0 ) // only try to add if the size is equal to window index - 1 -> the lower winder side
				{
					//try to Do Addition;
					Add(s,prob,proba);
					if (nv + nh == w)
					{
						fu++; // if successfully added particle, update fu 
						// cout << "!"<< endl;
					}
				}
			}

			// ============================Deletion=============================
			else 
			{
				if (size == w) // only try to delete if the size is equal to the window index -> the upper window side
				{
					//Do deletion;
					Del(s,prob,probd,size);
					if (nv+nh == w - 1.0)
					{
						fl++;//if successfully deleted particle, update fl
					}
				}			
			}			
		}

		// ======================= Update the upper window side ================================
        WF[w] += log(fu/fl);
        cout << fl<<"  "<<fu <<" " <<nv <<endl;
		// ======================= Print out the data into terminal =============================================		
		cout <<"Window: "<< w <<" : "<<"W("<<w<<" : lower) = "<< WF[w-1]<<" "<<"W("<<w<<" : Upper) = "<< WF[w] << endl;
	    w++; // switch into the next window

	}

	for(int i = 0; i< V/K+1; i++)
	{
		sh<<WF[i]<<endl;
	}

	ofstream myfile ("SUSWeight_function.txt");
	myfile.precision(20);
	string data2 = sh.str();
	myfile << data2;
	myfile.close();

	return WF;  
}


void MC::plot(const vector<HR>& VRodlist, const vector<HR>& HRodlist)
{
	stringstream stv,sth;

	for (int i = 0; i < VRodlist.size(); i++)
	{
		stv<< VRodlist[i].getX() << "   "<< VRodlist[i].getY()<<endl;
	}

	ofstream myfilev ("2dplotv.txt");
	string datav = stv.str();
	myfilev << datav;
	myfilev.close();

	for (int j = 0; j < HRodlist.size(); j++)
	{
		sth<< HRodlist[j].getX() << "   "<< HRodlist[j].getY() <<endl;
	}

	ofstream myfileh ("2dploth.txt");
	string datah = sth.str();
	myfileh << datah;
	myfileh.close();
}





