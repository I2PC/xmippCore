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

#ifndef XMIPPCORE_CORE_XMIPP_MARSAGLIA_H_
#define XMIPPCORE_CORE_XMIPP_MARSAGLIA_H_

/** Marsaglia class.
 *
 * Marsaglia random functions
 *
 * These functions are not ported to Python.
 *
 * This class has been designed for those programs that  have to have a large
 * set of random numbers but do not have time to generate them properly on the
 * fly. The idea is to have  a large file with lots of random numbers in it and
 * to store some of them in memory and retrieve (from memory), as many times as
 * needed, in a random fashion.
 *
 * As source of the numbers I recommend the "Marsaglia's Random Number CDROM"
 * available, for free, at your favorite web site. (Search for it in any search
 * engine  and you will get tons of hits.) The following is from the CD
 * description:
 *
 * "This CDROM will contain 4.8 billion random bits. They were produced by a
 * combination of several of the best deterministic random number generators
 * (RNG's), together with three sources of white noise, as well as black noise
 * (from a rap music digital recording). My intent is to provide an unassailable
 * source for those who absolutely positively have to have a large, reliable set
 * of random numbers for serious simulation (Monte Carlo) studies."
 *
 * When the random numbers are doubles, by default they are in the interval [0,1]
 *
 * This class is not ported to Python.
 */
template <typename T>
class Marsaglia
{
private:
    char* random_vector; // read the data right here
    T* T_random_vector;
    long pointer_in_memory;
    FileName fn_rand;
    long vector_size;
    long Number_of_Numbers;

public:
    /** Constructor
     * @ingroup Marsaglia
     *
     * @code
     * Marsaglia rodalnino("masaglia", 1000000, 34);
     * @endcode
     *
     * M_max (optional) is the magnitude of the maximum value of the random
     * number (exclusive), therefore must be positive
     */
    Marsaglia(const FileName& fn_in, int No_Numbers)
    {
        Init(fn_in, No_Numbers);
    }

    /// Empty constructor
    Marsaglia()
    {}

    /** You may use init for reading another set of random numbers
     */
    void Init(const FileName& fn_in, int No_Numbers)
    {
        int Type_size; // sizeof(type)

        pointer_in_memory = 0;
        Number_of_Numbers = No_Numbers; // initialize class variable
        Type_size = sizeof(T);

        std::ifstream in(fn_in.c_str());
        in.seekg(0, std::ios::end); // End of file
        std::streampos sp = in.tellg(); // Size of file
        if (sp < Number_of_Numbers * Type_size)
            REPORT_ERROR(ERR_IO_SIZE, (std::string) "Marsaglia::Init: File " + fn_in +
                         "is too small");
        else
        {
            // get a random number to set the file pointer at a random position
            randomize_random_generator();

            random_vector = new char[(Number_of_Numbers * Type_size)];
            T_random_vector = (T*) random_vector;
            in.seekg((std::streampos) FLOOR(rnd_unif(0.f, (double)(sp -
                                            (std::streamoff)(Number_of_Numbers * Type_size)))), std::ios::beg);
            in.read(random_vector, (Number_of_Numbers * Type_size));

            in.close();
        }
        if (typeid(double) == typeid(T))
            Verify_double();
    }

    /** Get a random number from the memory
     *
     * If you are at the end of the stream the pointer will be radomly moved
     * before stracting the number.
     */
    T Get_One_Number()
    {
        if (pointer_in_memory >= Number_of_Numbers)
            pointer_in_memory = (int) FLOOR(rnd_unif(0.f, (double)
                                            (Number_of_Numbers - 1)));
        return (T_random_vector[pointer_in_memory++]);
    }

    /** Calculate random vector log (use only with doubles)
     */
    void Marsaglia_log()
    {
        if (typeid(double) != typeid(T) && typeid(double) != typeid(T))
            REPORT_ERROR(ERR_TYPE_INCORRECT,
                         "Marsaglia: I do not know how to calculate integer logs");

        for (int hh = 0; hh < Number_of_Numbers; hh++)
            if (T_random_vector[hh] == 0.)
                T_random_vector[hh] = -1e+20f;
            else
                T_random_vector[hh] = log(T_random_vector[hh]);
    }

    /** Multiply random vector by constant
     */
    void mul(T mul_cte)
    {
        for (int hh = 0; hh < Number_of_Numbers; hh++)
            T_random_vector[hh] *= mul_cte;
    }

    /** Calculate mod of random vector, only make sense with integers
     */
    void operator&= (T mod_cte)
    {
        for (int hh = 0; hh < Number_of_Numbers; hh++)
            T_random_vector[hh] &= mod_cte;
    }

    /** Add a constant
     */
    void add(T add_cte)
    {
        for (int hh = 0; hh < Number_of_Numbers; hh++)
            T_random_vector[hh] += add_cte;
    }

    /** Set Maximum value (only valid for integers)
     */
    void M_max(const FileName& fn_in, T m_max)
    {
        int Type_size; // sizeof(type)
        Type_size = sizeof(T);

        std::ifstream in(fn_in.c_str());
        in.seekg(0, std::ios::end);              // End of file
        std::streampos sp = in.tellg();     // Size of file
        T power_of_2 = (T) NEXT_POWER_OF_2(m_max);
        if (power_of_2 == m_max)
            power_of_2 = (T) NEXT_POWER_OF_2(m_max + 1);
        T mask = power_of_2 - 1;
        T aux_number;

        // get a random number to set the file pointer at a random position
        in.seekg((std::streampos) FLOOR(rnd_unif(0.f, (double)(sp -
                                        (std::streamoff)(Number_of_Numbers*Type_size)))), std::ios::beg);
        for (int ii = 0; ii < Number_of_Numbers;)
        {
            aux_number = T_random_vector[ii];
            aux_number &= mask;
            if (aux_number > m_max ||
                (T_random_vector[ii] <= 0) && (aux_number == 0))
            {
                if (in.eof())
                    in.seekg((std::streampos) FLOOR(rnd_unif(0.f, (double)
                                                    (sp - (std::streamoff)(Number_of_Numbers*Type_size)))),
                             std::ios::beg);
                in.read((char*) &(T_random_vector[ii]), Type_size);
            }
            else
            {
                T_random_vector[ii] = aux_number * (T) SGN(T_random_vector[ii]);
                ii++;
            }
        }
        in.close();
    }

private:
    /** Verify double
     *
     * Be aware that Marsaglia reads blindly the data, therefore if the type
     * double is selected several of the "random" numbers may not be valid (the
     * number are created from a source of random bits and although 4 random
     * bits are one random integer, four random bits may not be a valid double).
     * If Marsaglia is "double" The constructor will run the following function
     * that will fix the problem
     */
    void Verify_double()
    {
        unsigned int* int_random_vector;
        long long MaxInteger;
        static_assert(sizeof(double) != sizeof(int), "I do not know how to make the double correction");
        MaxInteger = (long long) pow(2.0, sizeof(unsigned int) * 8.0);
        int_random_vector = (unsigned int*) random_vector;
        for (int hh = 0; hh < Number_of_Numbers; hh++)
            T_random_vector[hh] = (T)((double) int_random_vector[hh] /
                                      (double) MaxInteger);
    }
};

#endif /* XMIPPCORE_CORE_XMIPP_MARSAGLIA_H_ */
