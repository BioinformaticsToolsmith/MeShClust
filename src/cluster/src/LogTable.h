#ifndef LOGTABLE_H
#define LOGTABLE_H

#include <stdint.h>
#include <vector>

#define TBLSIZE 1000000
class LogTable {
public:
	LogTable();
	LogTable(uintmax_t _size, double imax=2);
	~LogTable();
	double at(double d) const;
	double operator[](double d) const;
private:
	double map[TBLSIZE];

	const double coeff;
};
#endif
