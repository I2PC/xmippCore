/***************************************************************************
 *
 * Authors:     David Strelak (davidstrelak@gmail.com)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef CORE_UTILS_TIME_UTILS_H_
#define CORE_UTILS_TIME_UTILS_H_

#include <chrono>
#include <utility>
#include <iostream>

/**
 * This is a utility class able to measure duration of a block of code.
 * It can be used for benchmarking parts of the code by wrapping it using lambda
 * expressions.
 *
 * Unless explicitly stated, its method can be used similarly to the following text:
 * timeUtils::measureTimeMs([&]{
 *  code_to_measure
 * });
 *
 * Typical usage would be to wrap several parts of the code with the 'report'
 * methods, disable reporting, and re-enable it once the code reaches the proper
 * location:
 *
 * void foo() {
 *  timeUtils::reportTimeMs("foo", [&]{
 *      // original code
 *  });
 * }
 *
 * int main() {
 *  timeUtils::disableReport();
 *  // some code, including e.g. foo()
 *  timeUtils::enableReport();
 *  // more code, including e.g. foo()
 * }
 *
 * Reporting methods are not thread-safe.
 */
class timeUtils
{
public:

    /**
     * Function to measure execution time of some code
     * @param ToDuration 'resolution' of the time, e.g. std::chrono::seconds
     * @param F function to run
     * @param Args arguments of the function
     */
    template<typename ToDuration, typename F, typename ...Args>
    static typename ToDuration::rep measureTime(F &&func, Args &&...args) {
        using namespace std::chrono;
        auto start = high_resolution_clock::now();
        std::forward<decltype(func)>(func)(std::forward<Args>(args)...);
        return duration_cast<ToDuration>(high_resolution_clock::now() - start).count();
    }

    /**
     * Function to measure and report execution time of some code
     * @param ToDuration 'resolution' of the time, e.g. std::chrono::seconds
     * @param unit to show in the report (e.g. 's' for std::chrono::seconds)
     * @param funName identification of the block
     * @param F function to run
     * @param Args arguments of the function
     */
    template<typename ToDuration, typename F, typename... Args>
    static void reportTime(const std::string &unit,
            const std::string &funName, F &&func, Args &&...args) {
        ++indent_counter;
        auto duration = measureTime<ToDuration>(func, args...);
        if (isReportEnabled) {
            reportIndent(indent_counter);
            std::cout << funName << separator << duration << separator << unit << "\n";
        }
        --indent_counter;
    }

    /**
     * Helper function, measuring time of some code in seconds
     * @param F function to run
     * @param Args arguments of the function
     */
    template<typename F, typename... Args>
    static std::chrono::seconds::rep measureTimeS(F &&func, Args &&...args) {
        return measureTime<std::chrono::seconds>(func, args...);
    }

    /**
     * Helper function, measuring time of some code in milliseconds
     * @param F function to run
     * @param Args arguments of the function
     */
    template<typename F, typename... Args>
    static std::chrono::milliseconds::rep measureTimeMs(F &&func, Args &&...args) {
        return measureTime<std::chrono::milliseconds>(func, args...);
    }

    /**
     * Helper function, measuring time of some code in seconds and reporting
     * it to std::cout.
     * @param funName some name the block (for identification purposes)
     * @param F funcion to run
     * @param Args arguments of the function
     */
    template<typename F, typename... Args>
    static void reportTimeS(const std::string &funName, F &&func, Args &&...args) {
        reportTime<std::chrono::seconds>("s", funName, func, args...);
    }

    /**
     * Helper function, measuring time of some code in milliseconds and reporting
     * it to std::cout.
     * @param funName some name the block (for identification purposes)
     * @param F funcion to run
     * @param Args arguments of the function
     */
    template<typename F, typename... Args>
    static void reportTimeMs(const std::string &funName, F &&func, Args &&...args) {
        reportTime<std::chrono::milliseconds>("ms", funName, func, args...);
    }

    /**
     * Helper function, measuring time of some code in microseconds and reporting
     * it to std::cout.
     * @param funName some name the block (for identification purposes)
     * @param F funcion to run
     * @param Args arguments of the function
     */
    template<typename F, typename... Args>
    static void reportTimeUs(const std::string &funName, F &&func, Args &&...args) {
        reportTime<std::chrono::microseconds>("us", funName, func, args...);
    }

    static void enableReport() {
        isReportEnabled = true;
    }

    static void disableReport() {
        isReportEnabled = false;
    }

    static void setSeparator(const std::string &s) {
        separator = s;
    }

private:
    /**
     * This method is responsible for indenting the current block of code
     * in the report fucntions
     */
    static void reportIndent(int level) {
        if (0 == level) {
            return;
        }
        std::cout << std::string((level - 1) * 2, ' ');
        std::cout << " » ";
    }

    static int indent_counter;
    static bool isReportEnabled;
    static std::string separator;

}; // timeUtils

#endif /* CORE_UTILS_TIME_UTILS_H_ */
