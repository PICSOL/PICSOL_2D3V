#include "stdafx.h"
#include "Vector3D.h"

template <class DataType>
Vector3D<DataType> Vector3D<DataType>::operator*(const Vector3D<DataType> &t)const
{
	Vector3D<DataType> cross_product;
	cross_product.x = y * t.z - z * t.y;
	cross_product.y = z * t.x - x * t.z;
	cross_product.z = x * t.y - y * t.x;

	return cross_product;
}

template <class DataType>
Vector3D<DataType> Vector3D<DataType>::operator*(double t)const
{
	Vector3D<DataType> multiply_product;
	multiply_product.x = x * t;
	multiply_product.y = y * t;
	multiply_product.z = z * t;

	return multiply_product;
}

template <class DataType>
Vector3D<DataType> Vector3D<DataType>::operator/(double t)const
{
	Vector3D<DataType> divide_product;
	divide_product.x = x / t;
	divide_product.y = y / t;
	divide_product.z = z / t;

	return divide_product;
}


template <class DataType>
Vector3D<DataType> Vector3D<DataType>::operator+(const Vector3D<DataType> &t)const
{
	Vector3D<DataType> add_product;
	add_product.x = x + t.x;
	add_product.y = y + t.y;
	add_product.z = z + t.z;

	return add_product;
}

template <class DataType>
Vector3D<DataType> Vector3D<DataType>::operator-(const Vector3D<DataType> &t)const
{
	Vector3D<DataType> sub_product;
	sub_product.x = x - t.x;
	sub_product.y = y - t.y;
	sub_product.z = z - t.z;

	return sub_product;
}

template <class DataType>
double Vector3D<DataType>::dot(const Vector3D<DataType> &t) const
{
	return x * t.x + y * t.y + z * t.z;
}

template <class DataType>
void Vector3D<DataType>::multiply(double a)
{
	x *= a;
	y *= a;
	z *= a;
}

template <class DataType>
double Vector3D<DataType>::power() const
{
	return x * x + y * y + z * z;
}