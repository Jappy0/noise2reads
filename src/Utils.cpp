#include "Utils.hpp"
#include "LoggingLevels.hpp"
#include <iostream>
#include <fstream>
#include <ctime>
#include <iostream>
#include <unordered_map>
#include <set>
#include <map>
#include <thread>
#include <unistd.h>
#include <functional>
#include <chrono>
#include <iomanip>
#include <sstream>
// #include <format>
#include <boost/format.hpp>

using namespace std;

Utils& Utils::getInstance() {
    static Utils instance;
    return instance;
}

void Utils::logger(int log_level, const std::string& message){
    if (log_level < LOG_LEVEL) return;

    // Get the current time and format it as a string
    time_t rawTime;
    time(&rawTime);
    struct tm * timeInfo = localtime(&rawTime);
    char timeString[20];
    strftime(timeString, 20, "%Y-%m-%d %H:%M:%S", timeInfo);

    std::string log_color;
    std::string log_prefix;
    switch (log_level) {
    case LOG_LEVEL_DEBUG:
        log_color = COLOR_DEBUG;
        log_prefix = "DEBUG: ";
        break;
    case LOG_LEVEL_INFO:
        log_color = COLOR_INFO;
        log_prefix = "INFO: ";
        break;
    case LOG_LEVEL_WARNING:
        log_color = COLOR_WARNING;
        log_prefix = "WARNING: ";
        break;
    case LOG_LEVEL_ERROR:
        log_color = COLOR_ERROR;
        log_prefix = "ERROR: ";
        break;
    case LOG_LEVEL_CRITICAL:
        log_color = COLOR_CRITICAL;
        log_prefix = "CRITICAL ERROR: ";
        break;
    default:
        log_color = COLOR_RESET;
        break;
    }
    logFile << timeString << ": " << log_prefix << message << endl;
    std::cout << timeString << ": " << log_color << log_prefix << COLOR_RESET << message << std::endl;
    // if(logFile.is_open()){
    //     logFile << timeString << ": " << log_prefix << message << endl;
    //     logFile.close();
    // }else {
    //     std::cerr << "\033[1;31mERROR: reads2graph log file opened failed.\033[0m" << std::endl;
    // }
}

Utils::Utils() {
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm now_tm = *std::localtime(&now_c);

    std::ostringstream oss;
    oss << std::put_time(&now_tm, "_%Y%m%d_%H%M%S.log");

    // Get the current time and format it as a string
    time_t rawTime;
    time(&rawTime);
    struct tm * timeInfo = localtime(&rawTime);
    char timeString[20];
    strftime(timeString, 20, "%Y-%m-%d %H:%M:%S", timeInfo);

    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL) {
        std::string filePath = std::string(cwd) + "/reads2graph" + oss.str();
        logFile.open(filePath, std::ios::app);
    } else {
        std::cerr << timeString << ": "<< "Error: unable to get current working directory." << std::endl;
    }
}

