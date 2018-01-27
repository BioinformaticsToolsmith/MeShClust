/* -*- C++ -*-
 *
 * Runner.cpp
 *
 * Author: Benjamin T James
 */
#include <vector>
#include <algorithm>
#include <cmath>
#include <sys/stat.h>
#include <cstdlib>
#include "../../nonltr/ChromListMaker.h"
#include "../../utility/AffineId.h"
#include "Runner.h"
#include "Trainer.h"
#include "ClusterFactory.h"
#include "bvec.h"
#include "Progress.h"
#ifdef _OPENMP
#include <omp.h>
#endif
Runner::Runner(int argc, char **argv)
{
	get_opts(argc, argv);
	if (k == -1) {
		auto pr = find_k();
		k = pr.first;
	}
	if (similarity < 0.6) {
		align = true;
	}
	if (sample_size == 0) {
		sample_size = 1500;
	}
	srand(10);
}

int Runner::run()
{
	largest_count = 0;
	Progress progress(files.size(), "Reading in sequences");
	for (auto i = 0; i < files.size(); i++) {
		auto f = files.at(i);
		ChromListMaker maker(f);
		auto chromList = maker.makeChromOneDigitList();

		progress++;
//		cout << "Reading in sequences from " << f << "..." << endl;
		uint64_t local_largest_count = 0;
#pragma omp parallel for reduction(max:local_largest_count)
	        for (int i = 0; i < chromList->size(); i++) {
			std::vector<uint64_t> values;
			KmerHashTable<unsigned long, uint64_t> table(k, 1);
			ChromosomeOneDigit *chrom = dynamic_cast<ChromosomeOneDigit *>(chromList->at(i));
			fill_table<uint64_t>(table, chrom, values);
			uint64_t l_count = *std::max_element(values.begin(), values.end());
			if (l_count > local_largest_count) {
				local_largest_count = l_count;
			}
		}
		if (local_largest_count > largest_count) {
			largest_count = local_largest_count;
		}
	}
	progress.end();


	if (largest_count <= std::numeric_limits<uint8_t>::max()) {
		cout << "Using 8 bit histograms" << endl;
		return do_run<uint8_t>();
	} else if (largest_count <= std::numeric_limits<uint16_t>::max()) {
		cout << "Using 16 bit histograms" << endl;
		return do_run<uint16_t>();
	} else if (largest_count <= std::numeric_limits<uint32_t>::max()){
	       	cout << "Using 32 bit histograms" << endl;
		return do_run<uint32_t>();
	} else if (largest_count <= std::numeric_limits<uint64_t>::max()) {
	       	cout << "Using 64 bit histograms" << endl;
		return do_run<uint64_t>();
	} else {
		throw "Too big sequence";
	}
}

void usage(std::string progname)
{
	std::cout << "Usage: " << progname << " *.fasta [--id 0.90] [--kmer 3] [--delta 5] [--output output.clstr] [--iterations 20] [--align] [--sample 1500] [--pivot 40] [--threads TMAX]" << std::endl << std::endl;
	#ifndef VERSION
        #define VERSION "(undefined)"
        #endif
        std::cout << "Version " << VERSION << " compiled on " << __DATE__ << " " << __TIME__;
        #ifdef _OPENMP
        std::cout << " with OpenMP " << _OPENMP;
        #else
        std::cout << " without OpenMP";
        #endif
	std::cout << std::endl;


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
		} else if (arg == "-a" || arg == "--align") {
			align = true;
		} else if ((arg == "-s" || arg == "--sample") && i + 1 < argc) {
			sample_size = strtol(argv[i+1], NULL, 10);
			if (errno) {
				perror(argv[i+1]);
				exit(EXIT_FAILURE);
			} else if (sample_size <= 0) {
				fprintf(stderr, "Sample size must be greater than 0.\n");
				exit(EXIT_FAILURE);
			}
			i++;
		} else if ((arg == "-p" || arg == "--pivot") && i + 1 < argc) {
			pivots = strtol(argv[i+1], NULL, 10);
			if (errno) {
				perror(argv[i+1]);
				exit(EXIT_FAILURE);
			} else if (sample_size <= 0) {
				fprintf(stderr, "Points per pivot must be greater than 0.\n");
				exit(EXIT_FAILURE);
			}
			i++;
		} else if ((arg == "-t" || arg == "--threads") && i + 1 < argc) {
			try {
				std::string opt = argv[i+1];
				int threads = std::stoi(opt);
				if (threads <= 0) {
					throw std::invalid_argument("");
				}
				#ifdef _OPENMP
				omp_set_num_threads(threads);
				#endif
			} catch (std::exception e) {
				cerr << "Number of threads must be greater than 0." << endl;
				exit(1);
			}

			i++;

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
		} else if ((arg == "-i" || arg == "--iter" || arg == "--iterations") && i + 1 < argc) {
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
				usage(*argv);
				exit(EXIT_FAILURE);
			}
		}
	}
	if (files.empty()) {
		usage(*argv);
		exit(EXIT_FAILURE);
	}
}

pair<int,uint64_t> Runner::find_k()
{
	unsigned long long count = 0, length = 0, largest_count = 0;
        uint64_t longest_seq = 0;
	uintmax_t num_sequences = 0;
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
			num_sequences++;

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


double global_mat[4][4] = {{1, -1, -1, -1},
			   {-1, 1, -1, -1},
			   {-1, -1, 1, -1},
			   {-1, -1, -1, 1}};
double global_sigma = -2;
double global_epsilon = -1;

void test()
{
	std::vector<
		std::pair<std::string,
			  std::string> > aligns = {
		{"GATCTCAG","GACAG"},
		{"GACAG","GATCAG"},
		{"GGAACCTT", "GGCCAATT"},
		{"GATCCATTACCG", "GATATTACCTT"},
		{"AGATGGTGCACGAACCGCGATTTGATGAATAACCTATTCGAACAGATTCCACCCCGTACTTAGATTCCACGGTAACAGTG",
		 "AGATGGTgaCggacccaTTTaagAATtAACCTAcTCGacAGAtTCCAcCtCCGtctaGATTCCACGGTacAaagTGAAGG"}
	};
	for (auto pr : aligns) {
		AffineId aid(pr.first.c_str(), 0, pr.first.length()-1,
			     pr.second.c_str(), 0, pr.second.length()-1);
		cout << aid.getAlign() << endl;
	}
	//exit(0);
}
template<class T>
int Runner::do_run()
{
	using pvec = vector<Point<T> *>;
	using pmap = map<Point<T>*, pvec*>;

	ClusterFactory<T> factory(k);
	auto points = factory.build_points(files, [&](nonltr::ChromosomeOneDigit *p){ return factory.get_divergence_point(p); });
	Trainer<T> tr(points, sample_size, largest_count, similarity, pivots, global_mat, global_sigma, global_epsilon, align ? 0 : k);
	tr.train();
	vector<uint64_t> lengths;
	for (Point<T>* p : points) {
		if (!align) {
			p->set_data_str("");
		}
		lengths.push_back(p->get_length());
	}
	// Initializing BVec
	bvec<T> bv(lengths, 1000);
	lengths.clear();
	// Inserting points into BVec
	uint64_t idx = 0;
	for (Point<T>* p : points) {
		p->set_id(idx++);
		bv.insert(p);
	}
	bv.insert_finalize();
//	cout << "bv size: " << bv.report() << endl;
	// Point<T>* mid = points[points.size()/2];
	// auto rng = bv.get_range(mid->get_length() * 0.99,
	// 			mid->get_length() / 0.99);
	// auto begin = bv.iter(rng.first);
	// auto end = bv.iter(rng.second);
	// size_t before = bv.report();
	// for (int i = 0; i < 1; i++) {
	// 		bool is_min = false;
	// 		Point<T>* p = tr.get_close(mid, begin, end, is_min);
	// 		size_t after = bv.report();
	// 		if (is_min) {
	// 			string expr = (after + 1 == before) ? "true" : "false";
	// 			if (expr == "false") {
	// 				throw expr;
	// 			}
	// 			cout << expr << endl;
	// 			cout << "is min" << endl;
	// 		} else {
	// 			cout << "is not min" << endl;
	// 		}
	// }
	factory.MS(bv, bandwidth, similarity, tr, output, iterations, delta);
	return 0;
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
