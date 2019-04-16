#pragma once

/*
Vector3D - data structure to store 3D vector
*/

template <class DataType>
class Vector3D
{
public:

	/* variable list */
	DataType x;              // x - component
	DataType y;              // y - component
	DataType z;              // z - component

    /* function list */
	Vector3D() :x(0), y(0), z(0) {}                                           // constructor function 
	Vector3D(DataType x1, DataType y1, DataType z1) :x(x1), y(y1), z(z1){}    // alternative constructor function
	Vector3D<DataType> operator*(const Vector3D<DataType> &t) const;
	Vector3D<DataType> operator*(double t) const;
	Vector3D<DataType> operator/(double t) const;
	Vector3D<DataType> operator+(const Vector3D<DataType> &t) const;
	Vector3D<DataType> operator-(const Vector3D<DataType> &t) const;
	double dot(const Vector3D<DataType> &t) const;
	void multiply(double a);
	double power() const;
	 
};

