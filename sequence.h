/**
 * @file
 * @brief Functions and structures for DNA sequences
 *
 */
#pragma once

#include <memory>
#include <string>
#include <vector>

class sequence
{
	std::string name{};
	std::shared_ptr<std::string> nucl{};
	size_t index{0}, length{0};

  public:
	sequence() = default;
	sequence(std::string, std::string) noexcept;

	auto size() const noexcept
	{
		return length;
	}

	auto get_name() const
	{
		if (index == 0 || length == 0) {
			return name;
		}

		return name + " (" + std::to_string(index) + ".." +
			   std::to_string(index + length) + ")";
	}

	auto get_nucl() const
	{
		return nucl.get()->substr(index, length);
	}

	auto begin() const
	{
		return nucl.get()->begin() + index;
	}

	auto end() const
	{
		return nucl.get()->begin() + index + length;
	}

	auto sub(size_t new_index, size_t new_length)
	{
		auto that = *this; // non-const
		that.index += new_index;
		that.length = new_length;
		return that;
	}

	auto c_str() const
	{
		return nucl.get()->c_str() + index;
	}
};

std::string reverse(const std::string &);
double gc_content(const std::string &) noexcept;

class genome
{
  public:
	std::string name{};
	std::vector<sequence> contigs{};
	size_t joined_length{0};

	genome() = default;

	genome(std::string _name, std::vector<sequence> _contigs) noexcept
		: name{std::move(_name)}, contigs{std::move(_contigs)}
	{
		size_t total = contigs.size() - 1;

		for (const auto &contig : contigs) {
			total += contig.size();
		}

		joined_length = total;
	}
};
