// DWPC - Domain Wall Propagation Code
// Copyright (C) 2019-2026 Razvan-Vasile Ababei and Matt Ellis
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This program is free software: you can redistribute it and/or modify it under the terms of
// the GNU General Public License as published by the Free Software Foundation, either version 3
// of the License, or (at your option) any later version. It is distributed WITHOUT ANY WARRANTY;
// see the file LICENSE for details. If you use DWPC in published work, please cite it as
// described in CITATION.cff.

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
            T val(0);
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
