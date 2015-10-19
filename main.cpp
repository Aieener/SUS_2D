// My main file for HR_MO+
// Data: July27 2015
// Authur: Yuding Ai
#include "MC.h"
#include <iostream>
using namespace std;

int main()
{
	stringstream t;
	double start = clock();

	// ======================= Plotting the final config ========================
	vector<HR> VR,HR;
	MC m(1E8L,8,64,64,20);
	m.MCRUN();
	VR = m.getVRodlist();
	HR = m.getHRodlist();
	m.plot(VR,HR);
	// ======================= get data for N vs z ========================
	// m.Zvs_();
	// ======================= end of simulation, print out the time =======
	double end = clock();
	cout <<"This simulation takes "<< (double(end-start)/CLOCKS_PER_SEC)<<" sec."<<endl;
	t <<"This simulation takes "<< (double(end-start)/CLOCKS_PER_SEC) <<" sec."<<endl;
	ofstream myfile ("time.txt");
	string datat = t.str();
	myfile << datat;
	myfile.close();	

	return 0;
}