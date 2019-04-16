#include "vectorField.h"

vectorField::vectorField(Param &param)
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
		vector<double> zeros(y_grid_num);
		xval.push_back(zeros);
		yval.push_back(zeros);
		zval.push_back(zeros);
	}
}

vectorField::~vectorField()
{
	/*
	STL vector doesn't require explicit memory recovery
	if you wish to use multi-dimension array or your own custom type,
	beware of the memory leak!
	*/
}

void vectorField::clear()
{
	for (int i = 0; i < x_grid_num; i++)
	{
		xval[i].assign(y_grid_num, 0);
		yval[i].assign(y_grid_num, 0);
		zval[i].assign(y_grid_num, 0);
	}
}

void vectorField::multiply(double constant)
{
	for (int i = 0; i < x_grid_num; i++)
	{
		for (int j = 0; j < y_grid_num; j++)
		{
			xval[i][j] *= constant;
			yval[i][j] *= constant;
			zval[i][j] *= constant;
		}
	}
}

void vectorField::equal(Vector3D<double> vec)
{
	for (int i = 0; i < x_grid_num; i++)
	{
		for (int j = 0; j < y_grid_num; j++)
		{
			xval[i][j] = vec.x;
			yval[i][j] = vec.y;
			zval[i][j] = vec.z;
		}
	}
}

void vectorField::sum(vectorField &a, vectorField &b)
{
	for (int i = 0; i < x_grid_num; i++)
	{
		for (int j = 0; j < y_grid_num; j++)
		{
			xval[i][j] = a.xval[i][j] + b.xval[i][j];
			yval[i][j] = a.yval[i][j] + b.yval[i][j];
			zval[i][j] = a.zval[i][j] + b.zval[i][j];
		}
	}
}

void vectorField::substract(vectorField &a)
{
	for (int i = 0; i < x_grid_num; i++)
	{
		for (int j = 0; j < y_grid_num; j++)
		{
			xval[i][j] -= a.xval[i][j];
			yval[i][j] -= a.yval[i][j];
			zval[i][j] -= a.zval[i][j];
		}
	}
}

void vectorField::add(Vector3D<double> vec)
{
	for (int i = 0; i < x_grid_num; i++)
	{
		for (int j = 0; j < y_grid_num; j++)
		{
			xval[i][j] += vec.x;
			yval[i][j] += vec.y;
			zval[i][j] += vec.z;
		}
	}
}

/* save field value in the given directory */
void vectorField::saveInfo(string directory1, string directory2, string directory3)
{
	ofstream output_file1(directory1.c_str());
	ofstream output_file2(directory2.c_str());
	ofstream output_file3(directory3.c_str());

	if (!output_file1.is_open()|| !output_file2.is_open()|| !output_file3.is_open())
	{
		cerr << " FIELD OUTPUT FILE OPEN ERROR !" << endl;
		system("pause");
	}
	output_file1.setf(ios::fixed);
	output_file1.precision(8);
	output_file2.setf(ios::fixed);
	output_file2.precision(8);
	output_file3.setf(ios::fixed);
	output_file3.precision(8);

	for (int i = 0; i < x_grid_num; i++)
	{
		for (int j = 0; j < y_grid_num; j++)
		{
			output_file1 << xval[i][j] << '\t';
			output_file2 << yval[i][j] << '\t';
			output_file3 << zval[i][j] << '\t';
		}
		output_file1 << endl;
		output_file2 << endl;
		output_file3 << endl;
	}
}