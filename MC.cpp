/*
* MC.cpp
* Simulation of 2-D Rds By GCMC
* Author: Yuding Ai
* Date: 2015.06.05
* *************************** MC implementation ********************************
* This simulation follows grand canonical ensemble.
* In such simulation, Only Energy and N of particles fluctuate, we accomplish EF by
* using both add/del and displacement moves, while NF by particle addition and 
* deletion moves;
* ******************************************************************************
*
* The Acceptance relation for Addition is given by: (Page 12 --M.S.Shell 2012)
*  -- Pad = min [1,((r*c/2)/(N + 1)) * exp{-beta*miu}]
*
* The Acceptance relation for Deletion is given by:
*  -- Pde = min [1,(N/(r*c/2))*exp{beta*miu}]
*
* where miu is the chemical potential
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
	step = ST;
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
	o = rand()%2;// change it to 1 for  lattice gas case

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
	// vector<HR> Rodlist;
	// Rodlist.clear();
	// Rodlist.insert( Rodlist.end(), VRodlist.begin(), VRodlist.end());//merge vertical first
	// Rodlist.insert( Rodlist.end(), HRodlist.begin(), HRodlist.end());//then merge hor
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


void MC::MCRUN()
{
	Cells s(c,r);

	stringstream st;

	double addordel; // the prob to decide either add or del;
	double probd,proba; // the acceptance prob of addition; proba = min(1.0,aaccp);
	double prob; // the prob to decide either accept add/del;
	double aaccp,daccp; 
	double Q; // the fraction of hor and ver particle;
	double tho; // the density 
	double AD;// addition and deletion fraction
	double size;
		
	srand(time(NULL));
	long int i = 0;
	// Histogram his(0, r*c/length, 1); // the histogram of nv

	//================================Start my MC simulation=================================
	while (i<step)
	{
		i++;
		// generate a random probability to decide either add or del;
		addordel = rand()%2;
		size = av+ah-dv-dh;

		// *****************define the probabilities ***********************************// I HAVE TO CHANGE IT FOR LATTICE GAS CASE!!!
		prob = ((double) rand() / (RAND_MAX)); 
		tho = double(length*size)/double(r*c);

		aaccp = z*double(r*c)/(double(nh+nv+1.0)*double(length));
		daccp = (double(nh+nv)*double(length))/(z*double(r*c));

		proba = min(1.0,aaccp);
		probd = min(1.0,daccp);

        // ===========================Addition ===================================
		if(addordel == 0) 
		{
			//Do Addition;
			Add(s,prob,proba);
		}

		// ============================Deletion=============================
		else
		{
			if (size != 0) // make sure there are rods to be del;
			{
				//Do deletion;
				Del(s,prob,probd,size);
			}			
		}

		// ======================= Record the datas =============================================		
        Q = (nv - nh)/(nh + nv);
		AD = (av+ah-dv-dh)/(av+ah+dv+dh);

		if (i%(step/10000) == 0)
		{
			// his.record(nv);
			st << i << "         " << Q <<"        "<< nv << "          "<< nh << "         "<< tho << "         "<< AD<< "         "<< endl;
			cout <<"Process: "<< ((10000*i)/step)/100.00 <<"%"<<"    "<<"SIZE: "<<av+ah-dv-dh<<"    "<<"# of Ver Rod: "<<nv<<"    "<<"# of Hor Rod: "<< nh <<"   "<<"Qis "<<Q <<"   "<<"tho is: "<<tho << endl;
		}
	}
	// Record the data into a txt file
	ofstream myfile3 ("dataplot.dat");
	string data = st.str();
	myfile3 << data;
	myfile3.close();
	// his.plot(0);
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

void MC::Zvs_()
{
	double z = 0;
	double H,V,tho,Q,miubeta,cmiubeta;	
	stringstream st;
	ofstream myfile("dataNvsZ.dat");
	for (int i =0; i < 500; i++)
	{
		MC m(1E9,1,100,100,z);
		z = double (double(10*i)/500.0);
		m.MCRUN();
		H =  m.getNh();
		V = m.getNv();
		tho = m.getTho();
		Q = (H - V)/(H + V);
		miubeta = log(z); // vink, lecture 7,8: page2
		cmiubeta = (log(tho) - log(1-tho));
		cout << i << endl;
		st << z <<"         " << H << "             "<< V<< "             "<<tho<< "             "<<Q<< "             "<< miubeta<< "             "<< cmiubeta << endl;
	}
	string data = st.str();
	myfile << data;
	myfile.close();
}




