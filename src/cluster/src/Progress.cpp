#include "Progress.h"
#include <iostream>
Progress::Progress(long num, std::string prefix_)
{
	pmax = num;
	ended = 0;
	pcur = 0;
	prefix = prefix_;
	barWidth = 70 - (prefix.size()+1);
	print();
}

void Progress::print()
{
	double prog = (double)pcur / pmax;
	std::cout << prefix << " [";
	int pos = barWidth * prog;
	for (int i = 0; i < barWidth; i++) {
		if (i < pos) {
			std::cout << "=";
		} else if (i == pos) {
			std::cout << ">";
		} else {
			std::cout << " ";
		}
	}
	std::cout << "] " << int(prog * 100.0) << " %\r";
	std::cout.flush();
}

void Progress::end()
{
	if (!ended) {
		pcur = pmax;
		print();
		std::cout << std::endl;
	}
	ended = true;
}

void Progress::operator++()
{
	pcur++;
	print();
}
void Progress::operator++(int)
{
	print();
	pcur++;
}


void Progress::operator+=(size_t num)
{
	pcur += num;
	print();
}
