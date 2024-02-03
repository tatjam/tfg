#pragma once

#include <fstream>

static void write_string_to_file(const std::string& path, const std::string& str)
{
	std::ofstream fout(path);
	fout << str;
	fout.close();
}

static std::string read_file_to_string(const std::string& path)
{
	std::ifstream fin(path);
	std::string content((std::istreambuf_iterator<char>(fin)),
						(std::istreambuf_iterator<char>()));
	fin.close();
	return content;
}