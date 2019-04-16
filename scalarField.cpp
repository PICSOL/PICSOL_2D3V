#include "scalarField.h"

scalarField::scalarField(Param &param)
{
	if (param.GRID_NUM_X <= 3 || param.GRID_NUM_Y <= 3)
	{
		cerr << " GRID NUMBER MUST BE A POSITIVE INTEGER GREATER THAN THREE!" << endl;
		system("pause");
	}
	this->x_grid_num = param.GRID_NUM_X;
	this->y_grid_num = param.GRID_NUM_Y;
	this->xbound1 = param.LBOUND;
	this->xbound2 = param.RBOUND;
	this->ybound1 = param.DBOUND;
	this->ybound2 = param.UBOUND;
	this->x_inc = param.DELTA_X;
	this->y_inc = param.DELTA_Y;

	for (int i = 0; i < x_grid_num; i++)
	{
		val.push_back(vector<double>(y_grid_num));
	}
}

scalarField::~scalarField()
{
	/* 
	 STL vector doesn't require explicit memory recovery
	 if you wish to use multi-dimension array or your own custom type,
	 beware of the memory leak!
	 */
}

void scalarField::clear()
{
	for (int i = 0; i < x_grid_num; i++)
	{
		val[i].assign(y_grid_num, 0);
	}
}

void scalarField::multiply(double constant)
{
	for (int i = 0; i < x_grid_num; i++)
	{
		for (int j = 0; j < y_grid_num; j++)
		{
			val[i][j] *= constant;
		}
	}
}

void scalarField::equal(double constant)
{
	for (int i = 0; i < x_grid_num; i++)
	{
		for (int j = 0; j < y_grid_num; j++)
		{
			val[i][j] = constant;
		}
	}
}

void scalarField::sum(scalarField& a, scalarField& b)
{
	for (int i = 0; i < x_grid_num; i++)
	{
		for (int j = 0; j < y_grid_num; j++)
		{
			val[i][j] = a.val[i][j] + b.val[i][j];
		}
	}
}

/* save field value in the given directory */
void scalarField::saveInfo(string directory)
{
	ofstream output_file(directory.c_str());
	if (!output_file.is_open())
	{
		cerr << " FIELD OUTPUT FILE OPEN ERROR !" << endl;
		system("pause");
	}
	output_file.setf(ios::fixed);
	output_file.precision(8);

	for (int i = 0; i < x_grid_num; i++)
	{
		for (int j = 0; j < y_grid_num; j++)
		{
            output_file << val[i][j] << '\t';
		}
		output_file << endl;
	}
}