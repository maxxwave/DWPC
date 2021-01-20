// Input map class

#ifndef __INPUT_MAP_H__
#define __INPUT_MAP_H__


#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <sstream>

class input_map_t {
    public:
        typedef std::string key_t;
        typedef std::string val_t;

        std::map<key_t, val_t> inp_map;

        std::string delim;
        std::string comment;

        input_map_t ()
        {
            delim = "=";
            comment = '#';
        }

        ~input_map_t () {};

        int read_file(const char* fname)
        {
            std::ifstream input(fname);
            if ( input.is_open() ) {
                std::string line;
                while ( getline( input, line)) {
                    if ( line[0] == '#' ) continue;

                    insert(line);

                }
                input.close();
            } else {
                std::cerr << "Failed to open file: " << fname << std::endl;
                return 1;
            }
            return 0;
        }

        void print()
        {
            for ( auto it = inp_map.begin(); it != inp_map.end(); it++)
                std::cout << it->first << " : " << it->second <<std::endl;
        }


        int insert(std::string &line)
        {

            int c_pos = line.find(comment);
            int n = line.substr(0,c_pos).find(delim);

            if ( n != std::string::npos) {
                std::string blank = " ";
                line.replace( n, blank.length(), blank);
                std::stringstream ss(line);
                key_t str1;
                val_t str2;
                ss >> str1 >> str2;
                inp_map[str1] = str2;

                return 1;
            }

            return 0;

        }


        template <typename T>
        T get(const char* char_key)
        {
            T val ;
            std::string key(char_key);
            auto it = inp_map.find(key);
            if ( it != inp_map.end() ) {
                std::stringstream ss(it->second);
                ss >> val;
            } else {
                std::cerr << "Warning: Input key not found. Key = " << char_key << std::endl;
            }

            return val;
        }

        template <typename T>
        T get(const char* char_key, T val)
        {
            std::string key(char_key);
            auto it = inp_map.find(key);
            if ( it != inp_map.end() ) {
                std::stringstream ss(it->second);
                ss >> val;
            } else {
                std::cerr << "Warning: Input key not found. Key = " << char_key << " , using default = " << val << std::endl;
            }

            return val;
        }

        bool exists(const char* char_key)
        {
            std::string key(char_key);
            auto it = inp_map.find(key);
            if ( it != inp_map.end() ) {
                return true;
            } else {
                return false;
            }
        }


};

#endif
