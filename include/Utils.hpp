/**
 * @Author: Pengyao Ping
 * @Date:   2023-02-07 14:10:30
 * @Last Modified by:   Pengyao Ping
 * @Last Modified time: 2023-03-05 12:52:34
 */
#pragma once

#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <memory>
#include <map>
#include <iostream>
#include <stdint.h>
#include <getopt.h>
#include <vector>
#include <functional>
#include <sharg/all.hpp> // includes all necessary headers

class Utils{
    public:
        static Utils& getInstance();
        
        void logger(int log_level, const std::string& message);
                
    private:
        std::ofstream logFile;

        Utils();
        Utils(const Utils&) = delete;
        Utils& operator=(const Utils&) = delete;        
};

#endif