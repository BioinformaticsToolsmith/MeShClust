/* -*- C++ -*-
 *
 * main.cpp
 *
 * Author: Benjamin T James
 */
#include "Runner.h"
int main(int argc, char **argv)
{

	// const rlim_t kStackSize = 10024 * 1024L * 1024L;   // min stack size = 1024 Mb
	// struct rlimit rl;
	// int result;

	// result = getrlimit(RLIMIT_STACK, &rl);
	// if (result == 0) {
	// 	if (rl.rlim_cur < kStackSize) {
	// 		rl.rlim_cur = kStackSize;
	// 		result = setrlimit(RLIMIT_STACK, &rl);
	// 		if (result != 0) {
	// 			fprintf(stderr, "setrlimit returned result = %d\n", result);
	// 		}
	// 	}
	// }
	Runner runner(argc, argv);
	return runner.run();
}
