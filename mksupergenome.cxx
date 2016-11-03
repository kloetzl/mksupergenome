/**
 * @file
 *
 * This is the main file. It contains functions to parse the commandline
 * arguments, read files etc.
 *
 * @brief The main file
 * @author Fabian Klötzl
 *
 * @section License
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <err.h>
#include <errno.h>
#include <getopt.h>
#include <string.h>
#include <unistd.h>

#include "config.h"
#include "global.h"
#include "io.h"
#include "process.h"
#include "sequence.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int FLAGS = flags::none;
int THREADS = 1;
double RANDOM_ANCHOR_PROP = 0.05;

void usage(void);
void version(void);

int main(int argc, char *argv[])
{
	int version_flag = 0;
	bool use_longest = true;
	std::string reference_file_name;

	static struct option long_options[] = {
		{"version", no_argument, &version_flag, 1},
		{"help", no_argument, NULL, 'h'},
		{"verbose", no_argument, NULL, 'v'},
		{"threads", required_argument, NULL, 't'},
		{"reference", required_argument, NULL, 'r'},
		{0, 0, 0, 0}};

#ifdef _OPENMP
	// Use all available processors by default.
	THREADS = omp_get_num_procs();
#endif

	// parse arguments
	while (1) {
		int c = getopt_long(argc, argv, "hvt:r:", long_options, NULL);

		if (c == -1) {
			break;
		}

		switch (c) {
			case 0: break;
			case 'h': usage(); break;
			case 'v':
				FLAGS |= FLAGS & flags::verbose ? flags::extra_verbose
												: flags::verbose;
				break;
			/*case 'r': {
				if (strcmp(optarg, "longest") == 0) break;
				use_longest = false;
				reference_file_name = optarg;
				break;
			}*/
			case 't': {
#ifdef _OPENMP
				errno = 0;
				char *end;
				long unsigned int threads = strtoul(optarg, &end, 10);

				if (errno || end == optarg || *end != '\0') {
					warnx("Expected a number for -t argument, but '%s' was "
						  "given. Ignoring -t argument.",
						  optarg);
					break;
				}

				if (threads > (long unsigned int)omp_get_num_procs()) {
					warnx(
						"The number of threads to be used, is greater then the "
						"number of available processors; Ignoring -t %lu "
						"argument.",
						threads);
					break;
				}

				THREADS = threads;
#else
				warnx(
					"This version of ugalign was built without OpenMP and thus "
					"does not support multi threading. Ignoring -t argument.");
#endif
				break;
			}
			case '?': /* intentional fall-through */
			default: usage(); break;
		}
	}

	if (version_flag) {
		version();
	}

	argc -= optind;
	argv += optind;

	// at least one file name must be given
	if (argc == 0) {
		errx(1, "At least one filename needs to be supplied.");
	}

	// at max `argc` many files have to be read.
	std::vector<genome> genomes;
	genomes.reserve(argc);

	// genome reference_genome;
	std::vector<std::string> file_names;
	file_names.reserve(argc);

	while (*argv) {
		file_names.push_back(*argv++);
	}

	if (file_names.size() < 2) {
		file_names.push_back("-"); // if no files are supplied, read from stdin
	}

	/*if (!use_longest) {
		auto it = std::find(file_names.begin(), file_names.end(),
							reference_file_name);
		if (it != file_names.end()) {
			file_names.erase(it);
		}
		reference_genome = read_genome(reference_file_name);
	}*/

	// read all genomes
	for (auto file_name : file_names) {
		genomes.push_back(read_genome(file_name));
	}

	// if no reference is given, pick the longest genome
	/*if (use_longest) {
		auto it = max_element(begin(genomes), end(genomes),
							  [](const genome &a, const genome &b) {
								  return a.joined_length < b.joined_length;
							  });
		reference_genome = *it;
		// delete longest_index
		genomes.erase(it);
	}*/

	// flatten the `genomes` array into `sequences`.
	std::vector<sequence> sequences{};
	auto inserter = std::back_inserter(sequences);
	for (auto &genome : genomes) {
		std::move(genome.contigs.begin(), genome.contigs.end(), inserter);
	}

	std::vector<sequence> supergenome{};

	auto it = max_element(begin(sequences), end(sequences),
						  [](const sequence &a, const sequence &b) {
							  return a.size() < b.size();
						  });
	auto ref = *it;
	supergenome.push_back(ref);
	sequences.erase(it);

	auto set = sequences;
	while (!set.empty()) {
		auto nm = filter(ref, set);
		std::cerr << "set: "<< set.size() << " nm: " << nm.size() << std::endl;
		if (nm.empty()) {
			break;
		}
		auto it = max_element(begin(nm), end(nm),
							  [](const sequence &a, const sequence &b) {
								  return a.size() < b.size();
							  });
		auto new_ref = *it;
		supergenome.push_back(new_ref);
		nm.erase(it);
		set = nm;
		ref = new_ref;
	}

	for (auto &seq : supergenome) {
		std::cout << seq.get_name() << "\n";
	}

	return 0;
}

/**@brief
 * Prints the usage to stdout and then exits successfully.
 */
void usage(void)
{
	const char str[] = {
		"Usage: mksupergenome [-lv] [-t INT] [-r longest|FILENAME] "
		"FILES...\n"
		"\tFILES... can be any sequence of FASTA files. If no files are "
		"supplied, stdin is used instead.\n"
		"Options:\n"
		"  -r longest|FILENAME   Use the sequence from FILENAME as reference; "
		"default: "
		"longest\n"
		"  -v, --verbose     Prints additional information\n"
#ifdef _OPENMP
		"  -t, --threads <INT> \n"
		"                    The number of threads to be used; by default, all "
		"available processors are used\n"
#endif
		"  -h, --help        Display this help and exit\n"
		"      --version     Output version information and acknowledgments\n"};

	printf("%s", str);
	exit(EXIT_SUCCESS);
}

/**
 * This function just prints the version string and then aborts
 * the program.
 */
void version(void)
{
	const char str[] = {
		"mksupergenome " VERSION "\n"
		"Copyright (C) 2016 Fabian Klötzl\n"
		"License GPLv3+: GNU GPL version 3 or later "
		"<http://gnu.org/licenses/gpl.html>\n"
		"This is free software: you are free to change and redistribute it.\n"
		"There is NO WARRANTY, to the extent permitted by law.\n\n"};
	printf("%s", str);
	exit(EXIT_SUCCESS);
}
