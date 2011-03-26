#include <sstream>
#include <cstring>
#include <cstdio>
#include <cstdlib>

#include "cisortho/string_tools.h"

using std::string;
using std::vector;

string join_token(vector<string> const& s, const char *sep){

	std::ostringstream res;
	int ns = (int)s.size();
	
	if (ns == 0) {} // do nothing
	else if (ns == 1) { res<<s[0]; }
	else {
		res<<s[0];
		for (int i=1; i < ns; i++) res<<sep<<s[i];
	}
	return res.str();
}

//assumes zero-terminated strings
string join_token(char **s, unsigned int ns, const char *sep){

	std::ostringstream res;
	
	if (ns == 0) {} // do nothing
	else if (ns == 1) { res<<s[0]; }
	else {
		res<<s[0];
		for (unsigned int i=1; i < ns; i++) res<<sep<<s[i];
	}
	return res.str();
}

vector<string> split_token(string const& str, string const& delim){
  size_t pos = 0, pos2 = 0;
  size_t max = str.size();
  vector<string> fields(0);
  while (pos <= max){
    pos2 = std::min(max, str.find(delim, pos));
    fields.push_back(str.substr(pos, pos2 - pos));
    pos = pos2 + delim.size();
  }
  return fields;
}

vector<string> split_token(char const* cstr, int num_chars, char const* cdelim){
  return split_token(string(cstr, num_chars), string(cdelim));
}


//convert all c-language escapes in a string to valid values
void convert_escapes(char * dest, const char* src)
{
    char * dest_ptr = dest;
    char const* src_ptr = src;

    const char * single_escapes = "\a\b\f\n\r\t\v\'\"\\\?";
    const char * single_escape_chars = "abfnrtv\'\"\\?";
    while (*src_ptr != '\0')
    {
        if (*src_ptr == '\\')
        {
            //eat and find the next character
            char escape = *(++src_ptr);
            char const* ptr = strchr(single_escape_chars, escape);
            if (ptr == NULL)
            {
                fprintf(stderr, "Error: don't understand escape sequence '\\%c'\n",
                        escape);
                exit(1);
            }
            ptrdiff_t ind = ptr - single_escape_chars;
            *(dest_ptr++) = single_escapes[ind];
            ++src_ptr;
        }
        else
        {
            //copy it verbatim
            *(dest_ptr++) = *(src_ptr++);
        }
    }
    //add null terminator to dest
    *dest_ptr = '\0';
}
