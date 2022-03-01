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
#include <unistd.h>
#include <cstring>
#include <libgen.h>
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
		sample_size = 3000;
	}
	srand(10);
}

int Runner::run()
{
	largest_count = 0;
	Progress progress(files.size(), "Reading in sequences");
	for (auto i = 0; i < files.size(); i++) {
		auto f = files.at(i);
		if (access(f.c_str(), F_OK) == -1) {
			cerr << "File \"" << f << "\" does not exist" << endl;
			exit(1);
		}
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
	std::cout << "Usage: " << progname << " *.fasta [--id 0.90] [--kmer 3] [--delta 5] [--output output.clstr] [--iterations 20] [--align] [--sample 3000] [--pivot 40] [--threads TMAX]" << std::endl << std::endl;
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
	std::string raw = R"(
If using long sequences or low (<60%) identity, use MeShClust2 (https://github.com/TulsaBioinformaticsToolsmith/MeShClust2).



The most important parameter, --id, controls the identity of the sequences.
    If the identity is below 60%, alignment is automatically used instead of k-mer measures.
    However, alignment can be forced with the --align parameter.

--kmer decides the size of the kmers. It is by default automatically decided by average sequence length,
       but if provided, MeShClust can speed up a little by not having to find the largest sequence length.
       Increasing kmer size can increase accuracy, but increases memory consumption fourfold.

--delta decides how many clusters are looked around in the final clustering stage.
	Increasing it creates more accuracy, but takes more time.

--output specifies the output file, in CD-HIT's CLSTR format

--iterations specifies how many iterations in the final stage of merging are done until convergence.

--align forces alignment to be used, which can be much slower than k-mer features, but is
	more accurate than using k-mer features to guess alignment.

--threads sets the number of threads to be used. By default OpenMP uses the number of available cores
	  on your machine, but this parameter overwrites that.

--sample selects the total number of sample pairs of sequences used for both training and testing.
	 1500 is the default value.

--pivot selects the maximum number of pairs selected from one pivot sequence. Increasing this means
	less pivots are available, but more pairs are selected for one sequence, which can lead to
	higher training accuracy. The default value is 40.

If the argument is not listed here, it is interpreted as an input file.


If you find this tool helpful, please cite:

James, Benjamin T. et al. (2018), MeShClust: an intelligent tool for clustering DNA sequences. Nucleic Acids Research, gky315.

)";
	std::cout << raw << endl;
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
			} else if (pivots <= 0) {
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
			} else if (delta < 0) {
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
	std::sort(files.begin(), files.end(), [](const std::string &a, const std::string &b) {
			char* as = strdup(a.c_str());
			char* bs = strdup(b.c_str());
			char* a_bn = basename(as);
			char* b_bn = basename(bs);
			bool ret = std::string(a_bn) < std::string(b_bn);
			free(as);
			free(bs);
			return ret;
		});
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
	ClusterFactory<T> factory(k);

	auto points = factory.build_points(files, [&](nonltr::ChromosomeOneDigit *p){ return factory.get_divergence_point(p); });
	Trainer<T> tr(points, sample_size, largest_count, similarity, pivots, global_mat, global_sigma, global_epsilon, align ? 0 : k);
	tr.train();
	vector<uint64_t> lengths;
	for (Point<T>* p : points) {
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
