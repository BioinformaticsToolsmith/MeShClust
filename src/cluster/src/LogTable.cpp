#include "LogTable.h"

#include <cmath>
#include <iostream>

LogTable::LogTable() : coeff(1000000 / 2)
{
	uintmax_t size = 1000000;
	double imax = 2;
//	map = new double[size];
	double lsize = log(size);
	for (uintmax_t i = 0; i < size; i++) {
		map[i] = log(imax * (i + 1)) - lsize;
	}
	std::cout << "dmax: " << coeff << std::endl;
}
LogTable::LogTable(uintmax_t size, double imax) : coeff(size / imax)
{
	//map = new double[size];
	double lsize = log(size);
	for (uintmax_t i = 0; i < size; i++) {
		map[i] = log(imax * (i + 1)) - lsize;
	}
	std::cout << "dmax: " << coeff << std::endl;
}

LogTable::~LogTable()
{
	//delete[] map;
}

double LogTable::at(double d) const
{
	size_t idx = d * coeff;
	return map[idx];
}
double LogTable::operator[](double d) const
{
	size_t index = d * coeff;
	return map[index];
}
