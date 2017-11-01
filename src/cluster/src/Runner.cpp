/* -*- C++ -*-
 *
 * Runner.cpp
 *
 * Author: Benjamin T James
 */
#include <vector>
#include <cmath>
#include <sys/stat.h>
#include <cstdlib>
#include "../../nonltr/ChromListMaker.h"
#include "Runner.h"
#include "Trainer.h"
#include "ClusterFactory.h"
#include "ConcreteMeanShift.h"
#include "MeanShiftContext.h"
Runner::Runner(int argc, char **argv)
{
	get_opts(argc, argv);
	auto pr = find_k();
	longest_length = pr.second;
	if (k == -1) {
		k = pr.first;
	}
	compute_bandwidth();
	srand(time(NULL));
}

int Runner::run() const
{
	if (longest_length <= std::numeric_limits<uint8_t>::max()) {
		return do_run<uint8_t>();
	} else if (longest_length <= std::numeric_limits<uint16_t>::max()) {
		return do_run<uint16_t>();
	} else if (longest_length <= std::numeric_limits<uint32_t>::max()){
		return do_run<uint32_t>();
	} else if (longest_length <= std::numeric_limits<uint64_t>::max()) {
		return do_run<uint64_t>();
	} else {
		throw "Too big sequence";
	}
}
void Runner::get_opts(int argc, char **argv)
{
	for (int i = 1; i < argc; i++) {
		string arg = argv[i];
		if (arg == "--id" && i + 1 < argc) {
			try {
				std::string opt = argv[i+1];
				similarity = std::stod(opt);
				if (similarity <= 0 || similarity >= 1) {
					throw std::invalid_argument("");
				}
			} catch(std::exception e) {
				cerr << "Similarity must be between 0 and 1" << endl;
				exit(EXIT_FAILURE);
			}
			similarity *= 100;
			i++;
		} else if ((arg == "-k" || arg == "--kmer") && i + 1 < argc) {
			k = strtol(argv[i+1], NULL, 10);
			if (errno) {
				perror(argv[i+1]);
				exit(EXIT_FAILURE);
			} else if (k <= 0) {
				fprintf(stderr, "K must be greater than 0.\n");
				exit(EXIT_FAILURE);
			}
			i++;
		} else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
			output = string(argv[i+1]);
			i++;
		} else if ((arg == "-t" || arg == "--threads") && i + 1 < argc) {
			try {
				std::string opt = argv[i+1];
				int threads = std::stoi(opt);
				if (threads <= 0) {
					throw std::invalid_argument("");
				}
			} catch (std::exception e) {
				cerr << "Number of threads must be greater than 0." << endl;
				exit(1);
			}



		} else if ((arg == "-d" || arg == "--delta") && i + 1 < argc) {
			delta = strtol(argv[i+1], NULL, 10);
			if (errno) {
				perror(argv[i+1]);
				exit(EXIT_FAILURE);
			} else if (delta <= 0) {
				fprintf(stderr, "Delta must be greater than 0.\n");
				exit(EXIT_FAILURE);
			}
			i++;
		} else if ((arg == "-t" || arg == "--iter" || arg == "--iterations") && i + 1 < argc) {
			iterations = strtol(argv[i+1], NULL, 10);
			if (errno) {
				perror(argv[i+1]);
				exit(EXIT_FAILURE);
			} else if (iterations <= 0) {
				fprintf(stderr, "Iterations must be greater than 0.\n");
				exit(EXIT_FAILURE);
			}
			i++;
		} else {
			struct stat st;
			stat(argv[i], &st);
			if (S_ISREG(st.st_mode)) {
				files.push_back(argv[i]);
			} else {
				cerr << "Usage: " << *argv << " *.fasta [--id 0.90] [--kmer 3] [--delta 5] [--output output.clstr] [--iterations 20]" << endl;
				exit(EXIT_FAILURE);
			}
		}
	}
	if (files.empty()) {
		cerr << "Usage: " << *argv << " *.fasta [--id 0.90] [--kmer 3] [--delta 5] [--output output.clstr] [--iterations 20]" << endl;
		exit(EXIT_FAILURE);
	}
}

pair<int,uint64_t> Runner::find_k() const
{
	unsigned long long count = 0, length = 0;
        uint64_t longest_seq = 0;
	for (auto f : files) {
		ChromListMaker maker(f);
		auto chromList = maker.makeChromOneDigitList();
		unsigned long long l = 0;
	        for (int i = 0; i < chromList->size(); i++) {
			ChromosomeOneDigit *chrom = dynamic_cast<ChromosomeOneDigit *>(chromList->at(i));
			auto sz = chrom->size();
			l += sz;
			if (sz > longest_seq) {
				longest_seq = sz;
			}

		}
		l /= chromList->size();
		length += l;
	}
	length /= files.size();

	int newk = ceil(log(length) / log(4)) - 1;
	cout << "avg length: " << length << endl;
	cout << "Recommended K: " << newk << endl;
	return make_pair(newk, longest_seq);
}

void Runner::compute_bandwidth()
{
	map<int, int> table;
	similarity /= 100;
	//similarity -= 0.05;
	// table[5] = round(similarity * -905.48 + 912.38);
	// table[4] = round(similarity * -623.78 + 628.58);
	// table[3] = round(similarity * -364.6  + 367.98);
	// table[2] = round(similarity * -221.41 + 222.31);
	// table[1] = round(similarity * -140.58 + 140.48);
	// for (auto kv : table) {
	// 	cout << "K: " << kv.first << " V: " << kv.second << endl;
	// }
	//k = min(5, k);
	cout << "From similarity " << similarity << endl;
//	int ret = table[k];
	int ret = (1.0 - similarity) * 10000;
	cout << "Using bandwidth " << ret << endl;
	bandwidth = ret;
}

double global_mat[4][4] = {{1, -1, -1, -1},
			   {-1, 1, -1, -1},
			   {-1, -1, 1, -1},
			   {-1, -1, -1, 1}};
double global_sigma = -2;
double global_epsilon = -1;
template<class T>
int Runner::do_run() const
{
	using pvec = vector<Point<T> *>;
	using pmap = map<Point<T>*, pvec*>;
	try {
		ClusterFactory<T> factory(k);
		auto points = factory.build_points(files, [&](nonltr::ChromosomeOneDigit *p){ return factory.get_divergence_point(p); });
		Trainer<T> tr(points, 1000, similarity, 40, global_mat, global_sigma, global_epsilon, k);
		tr.train();
		factory.MS(points, bandwidth, similarity, tr, output, iterations, delta);
		return 0;
	} catch (std::exception e) {
		std::cerr << e.what() << endl;
		return 1;
	}
	// factory.MS(points, bandwidth, similarity, output, iterations, delta);
	// return 0;
}


template<class T>
void Runner::print_output(const map<Point<T>*, vector<Point<T>*>*> &partition) const
{
	cout << "Printing output" << endl;
	std::ofstream ofs;
	ofs.open(output, std::ofstream::out);
	int counter = 0;
	for (auto const& kv : partition) {
		if (kv.second->size() == 0) {
			continue;
		}
		ofs << ">Cluster " << counter << endl;
		int pt = 0;
		for (auto p : *kv.second) {
			string s = p->get_header();
			ofs << pt << "\t"  << p->get_length() << "nt, " << s << "... " << endl;
//			string fa = am.get(p->get_id());
//			ofs << writefa(fa) << endl;
			pt++;
		}
		counter++;
	}
	ofs.close();
}
