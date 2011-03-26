#ifndef _STRING_TOOLS_H
#define _STRING_TOOLS_H

#include <string>
#include <vector>


//return a string of strings joined by sep
std::string join_token(std::vector<std::string> const& strings, const char *sep);
std::string join_token(char **strings, unsigned int ns, const char *sep);
std::vector<std::string> split_token(std::string const& str, std::string const& delim);
std::vector<std::string> split_token(char const* string, int num_chars, char const* delim);
void convert_escapes(char * dest, const char* src);

#endif // _STRING_TOOLS_H
