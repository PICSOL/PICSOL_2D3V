#include "Param.h"

Param::Param()
{
	ifstream in("param.txt", ios::in);
	char buffer[256];
	
	/* data storage path */
	in.getline(buffer, 256);
	PATH = string(buffer).erase(0,1);

	in.getline(buffer,256);
	SHAPE_FACTOR = atoi(buffer);

	/* simulation constant */
	in.getline(buffer, 256);
	MAX_ITERATION = atoi(buffer);

	in.getline(buffer, 256);
	MIN_ERROR  = atof(buffer);

	in.getline(buffer, 256);
	GRID_NUM_X  = atoi(buffer);

	in.getline(buffer, 256);
	GRID_NUM_Y = atoi(buffer);

	in.getline(buffer, 256);
	MAX_LOOP = atoi(buffer);
	
	in.getline(buffer, 256);
	RECORD_CHECK = atoi(buffer);

	in.getline(buffer, 256);
	INITIAL_PARTICLE = atoi(buffer);

	in.getline(buffer, 256);
	DELTA_T = atof(buffer);

	in.getline(buffer, 256);
	DELTA_Te = atof(buffer);

	in.getline(buffer, 256);
	DELTA_Ti = atof(buffer);

	in.getline(buffer, 256);
	DELTA_X = atof(buffer);

	in.getline(buffer, 256);
	DELTA_Y = atof(buffer);

	in.getline(buffer, 256);
	BX = atof(buffer);

	in.getline(buffer, 256);
	BY = atof(buffer);

	/* plasma parameters */
	in.getline(buffer, 256);
	DENSITY = atof(buffer);

	in.getline(buffer, 256);
	ETEMP = atof(buffer);

	in.getline(buffer, 256);
	ITEMP = atof(buffer);

	in.getline(buffer, 256);
	ELECTCHARGE = atof(buffer);

	in.getline(buffer, 256);
	IONCHARGE = atof(buffer);

	in.getline(buffer, 256);
	ELECTMASS = atof(buffer);

	in.getline(buffer, 256);
	IONMASS = atof(buffer);

	in.getline(buffer, 256);
	V_BIAS = atof(buffer);

	in.getline(buffer, 256);
	Y_0 = atof(buffer);

	in.getline(buffer, 256);
	L_0 = atof(buffer);

	in.getline(buffer, 256);
	L_1 = atof(buffer);

	if (in.fail()||in.bad())
	{
		cerr << "CONSTANT FILE LOADING ERROR!";
		system("pasuse");
	}

/*
	ifstream LBValue("leftBound.dat", ios::in);
	ifstream RBValue("rightBound.dat", ios::in);
	double boundValue;

	for (size_t i = 0; i < GRID_NUM_Y; i++)
	{
		LBValue >> boundValue;
		LB.push_back(boundValue);
		RBValue >> boundValue;
		RB.push_back(boundValue);
	}

	if (LBValue.fail() || LBValue.bad()|| RBValue.fail() || RBValue.bad())
	{
		cerr << "BOUNDARY VALUE LOADING ERROR!";
		system("pasuse");
	}
*/
	LBOUND = 0;
	RBOUND = DELTA_X * (GRID_NUM_X - 1);
	XLENGTH = RBOUND - LBOUND;

	DBOUND = 0;
	UBOUND = DELTA_Y * GRID_NUM_Y;
	YLENGTH = UBOUND - DBOUND;

	for (int i = 0; i < GRID_NUM_Y; i++)
	{
		LB.push_back(0.5 * V_BIAS *(1 - tanh(2 * e / L_1 * (abs(i - Y_0) - L_0))));
		RB.push_back(0);
	}

}

void Param::copyFile(string& input,string& output, string& dir)
{
	string outfile = dir + string("\\") + output;
	ifstream in(input.c_str(), ios::in | ios::binary);
	ofstream out(outfile.c_str(), ios::out | ios::binary);

	if (!in || !out)
	{
		cerr << "COPY FILE FAILURE£¡";
		system("pause");
	}
	
	char buffer[256];
	while (!in.eof())
	{
		in.read(buffer, 256);
		out.write(buffer, in.gcount());
	}

	in.close();
	out.close();
}

void Param::saveInput(string& destname)
{
	string parameter("param.txt");
	copyFile(parameter, parameter, destname);
	/*
	string LB("leftBound.dat");
	copyFile(LB, LB, destname);

	string RB("rightBound.dat");
	copyFile(RB, RB, destname);
	*/
}