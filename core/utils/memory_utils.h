/***************************************************************************
 *
 * Authors:    David Strelak (davidstrelak@gmail.com)
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

#ifndef CORE_UTILS_MEMORY_UTILS_H_
#define CORE_UTILS_MEMORY_UTILS_H_

#include <cstddef>
#include <stdlib.h>

namespace memoryUtils
{

    inline constexpr size_t operator"" _kB(unsigned long long int bytes) {
      return 1024 * bytes;
    }

    inline constexpr size_t operator"" _MB(unsigned long long int bytes) {
      return 1024 * 1024 * bytes;
    }

    inline constexpr size_t operator"" _GB(unsigned long long int bytes) {
      return 1024 * 1024 * 1024 * bytes;
    }

    inline constexpr double operator"" _kB(long double bytes) {
      return 1024. * bytes;
    }

    inline constexpr double operator"" _MB(long double bytes) {
      return 1024. * 1024. * bytes;
    }

    inline constexpr double operator"" _GB(long double bytes) {
      return 1024. * 1024. * 1024. * bytes;
    }


    inline void* page_aligned_alloc(size_t bytes) {
        return aligned_alloc(4096, bytes);
    }

    template<typename T>
    inline T* page_aligned_alloc(size_t elems, bool initToZero) {
        size_t bytes = elems * sizeof(T);
        auto p = (T*)page_aligned_alloc(elems * sizeof(T));
        if (initToZero) {
            memset(p, 0, bytes);
        }
        return p;
    }

    template<typename T>
    inline constexpr T kB(T bytes) {
      return  bytes / (T)1024;
    }

    template<typename T>
    inline constexpr T MB(T bytes) {
      return  bytes / ((T)1024 * T(1024));
    }

    template<typename T>
    inline constexpr T GB(T bytes) {
      return  bytes / ((T)1024 * (T)1024 * (T)1024);
    }
} // memoryUtils

#endif /* CORE_UTILS_MEMORY_UTILS_H_ */
