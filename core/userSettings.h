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

/**
 * This class is able to store/read user settings.
 * Use it to permanently store some data/setting.
 */

#ifndef USERSETTINGS_H_
#define USERSETTINGS_H_

#include <string>
#include <algorithm>
#include <map>
#include <list>
#include <fstream>
#include <iosfwd>
#include <sstream>
#include <typeinfo>
#include <mutex>

class UserSettings {
public:
    /**
     * Returns the only instance of this class
     */
    static UserSettings& get(const std::string &path = std::string()) {
        static std::mutex mtx;
        auto res = find_if(storages.begin(), storages.end(),
            [&path](const UserSettings& obj)
            {return obj.identifier.compare(path) == 0;});
        if (storages.end() == res) {
            // it seems that we don't have such a storage yet
            mtx.lock();
            // make sure that nobody did it meanwhile
            res = find_if(storages.begin(), storages.end(),
                [&path](const UserSettings& obj)
                {return obj.identifier.compare(path) == 0;});
            if (storages.end() == res) {
                storages.emplace_front(UserSettings(path));
                res = storages.begin();
            }
            mtx.unlock();
        }
        // now we have to have it
        return *res;
    }

    /**
     * Save data to HDD
     */
    bool store();

    /**
     * Reload data from HDD. All changes will be discarded
     */
    bool reload();

    /**
     * Insert new / update existing key-value pair.
     * Value will be converted to std::string using std::stringstream and
     * operator <<.
     * @param obj storing the key
     * @param key to use
     * @param value to store
     * @returns false if identical key-value pair already existed, true otherwise
     */
    template<typename T, typename U>
    bool insert(const T &obj, const std::string &key, const U &value) {
        auto fullKey = getFullName(obj, key);
        auto result = data.find(fullKey);
        std::stringstream ss;
        ss << value;
        std::string valStr = ss.str();
        bool newRec = result == data.end();
        bool update = !newRec && (valStr.compare(result->second) != 0);

        if (newRec || update) {
            data[fullKey] = valStr;
            wasChanged = true;
            return true;
        }
        return false;
    }

    /**
     * Finds value using key.
     * @param obj making query
     * @param key to search up
     * @param value to be set (if the record exists)
     * @returns true if record with the key exists, false otherwise
     */
    template<typename T, typename U>
    bool find(const T &obj, const std::string &key, U &value) {
        auto result = data.find(getFullName(obj, key));
        if (result != data.end()) {
            std::stringstream ss (result->second);
            ss >> value;
            return true;
        } else {
            return false;
        }
    }

    UserSettings(UserSettings&&) = default;                // Move construct

    /** Destructor, stores the data on change */
    ~UserSettings() { if (wasChanged) store(); };

protected:
    /** Constructor, loads the data from HDD */
    UserSettings(const std::string &path = std::string()) : wasChanged(false),
        identifier(path), path(path) {
        if (path.empty()) {
            this->path = std::string(getenv("HOME")) + "/.xmipp.settings";
        }
        reload();
    };

    /**
     * Generate full key name using the caller
     * @param obj of the caller
     * @param key to use
     */
    template<typename T>
    const std::string getFullName(const T &obj, const std::string &key) {
        std::stringstream ss;
        ss << typeid(obj).name() << key;
        return ss.str();
    }

private:
    // delete copy and move constructors and assign operators
    UserSettings(UserSettings const&) = delete;           // Copy construct
    UserSettings& operator=(UserSettings const&) = delete;   // Copy assign
    UserSettings& operator=(UserSettings &&) = delete;       // Move assign

    /** Actual records */
    std::map<std::string, std::string> data;

    /** If true, content of the records changed since last load */
    bool wasChanged;

    /** Delimiter used while storing data on HDD */
    char delim = '\t';

    /** Path to the file where data are stored */
    std::string path;

    /** Unique identifier of the instance */
    std::string identifier;

    static std::list<UserSettings> storages;
};

#endif /* USERSETTINGS_H_ */
