/**
 * @file
 * @brief Sequence utilities
 *
 * This file contains utility functions for working with DNA sequences.
 */

#include <memory>
#include <vector>

#include <ctype.h>
#include <err.h>
#include <limits.h>
#include <string.h>
// #include "global.h"
#include "sequence.h"

void normalize(std::string &) noexcept;

sequence::sequence(std::string name, std::string nucl_) noexcept
	: name{std::move(name)},
	  nucl{std::make_shared<std::string>(std::move(nucl_))}, index{0},
	  length{nucl.get()->size()}
{
	const size_t LENGTH_LIMIT = (INT_MAX - 1) / 2;
	if (this->size() > LENGTH_LIMIT) {
		warnx("The input sequence %s is too long. The technical limit is %zu.",
			  this->name.c_str(), LENGTH_LIMIT);
	}
}

/**
 * @brief Compute the reverse complement.
 * @param base - The master string.
 * @return The reverse complement.
 */
std::string reverse(const std::string &base)
{
	std::string ret{};
	ret.reserve(base.size());

	size_t len = base.size();
	auto str = new char[len + 1];

	for (size_t k = 0; k < len; k++) {
		char c = base[len - k - 1], d;

		// if (c & 8) {
		// 	d = c; // N or -
		// } else {
		d = c ^= (c & 2) ? 4 : 21; // ACGT
		// }

		str[k] = d;
	}

	ret.replace(0, len, str);
	delete[] str;

	return ret;
}

double gc_content(const std::string &seq) noexcept
{
	size_t gc = 0;
	size_t length = seq.size();

	for (size_t i = 0; i < length; i++) {
		char masked = seq[i] & 'G' & 'C';
		if (masked == ('G' & 'C')) {
			gc++;
		}
	}

	return static_cast<double>(gc) / length;
}
