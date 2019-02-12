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

#include "userSettings.h"
#include <iostream>

std::list<UserSettings> UserSettings::storages;

bool UserSettings::store() {
    std::ofstream output(path);
    if (output) {
        for (const auto& rec : data) {
            output << rec.first << delim << rec.second << "\n";
        }
        output.close();
        return true;
    }

    std::cerr<<"Error opening output file "<< path <<  std::endl;
    return false;
}

bool UserSettings::reload() {
    std::ifstream file(path);
    std::string line;

    if (!file) {
        std::cerr<<"File "<< path << " does NOT exist yet or cannot be read."
            << std::endl;
        return false;
    }
    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string key, value;
        std::getline(ss, key, delim);
        std::getline(ss, value);
        data[key] = value;
    }
    file.close();
    return true;
}
