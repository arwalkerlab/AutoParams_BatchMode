#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <locale>
#include <vector>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>
#include <iomanip>
#include <set>
#include <ctime>
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;

// Command Line Interaction
std::string GetSystemResponse(const char* cmd);
bool CheckFileExists(std::string filename);

// Check available tools on the system
bool CheckProgramExists(std::string program);

// Internal utilities
std::string trim_whitespace(const std::string& str, const std::string& whitespace = " \t");
void silent_shell(const char* cmd);
std::vector <std::string> string_split(std::string source, char delimiter);
std::string substring_replace(std::string original, std::string target, std::string replacement);
void DeleteIfEmpty(std::string);
std::string FormattedTimeStamp();
int count_element_in_array(std::vector<std::string> array,std::string element);
#endif