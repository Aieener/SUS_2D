/*
* hardrods.h
* 2-D lattice gas
* Author: Yuding Ai & Eelco Empting
* Date: 2015.06.08
*/
#ifndef HARDRODS_H
#define HARDRODS_H
#include <iostream>
#include <vector>
#include <string>

class HR
{
	private:
		int X; // the x index;
		int Y; // the y index;
		int Orientation; // the orientation of the rod; 2 options: up, right,
		int Length; // the length of the rod

	public:
		//class constants;
		static const int UP = 0;
		static const int RIGHT = 1;


	/*
	* Default constructor
	* @para     none
	* @set X = Y = Orientation = 0; Length = 1
	*/
	HR();

	/*
	* constructor
	* @para     int x, y, L, orientation
	*/
	HR(int x, int y, int l, int o);

	//**** Getters *****//
	int getX() const;
	int getY() const;
	int getLength() const;
	int getOrientation() const;
};

#endif /* HARDRODS_H */
