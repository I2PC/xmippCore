/***************************************************************************
 *
 * Authors: Martin Salinas Anton (martin.salinas@cnb.csic.es)
 *
 * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA
 *
 * All comments concerning this program package may be sent to the
 * e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/
/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>

#include "xmipp_image.h"
#include "xmipp_error.h"
#include <tiffio.h>
// Original libraries included from Relion repository. Keeping them here commented for now
// while the library is being fully integrated & tested so traceability is easier
//#include <src/image.h>

class EERRenderer {
	private:

	FileName fn_movie;

	bool ready;
	bool is_legacy;
	bool is_7bit;
	bool read_data;

	std::vector<long long> frame_starts, frame_sizes;
	unsigned char* buf;

	static const char EER_FOOTER_OK[];
	static const char EER_FOOTER_ERR[];
	static const int EER_IMAGE_WIDTH, EER_IMAGE_HEIGHT, EER_IMAGE_PIXELS;
	static const unsigned int EER_LEN_FOOTER;
	static const uint16_t TIFF_COMPRESSION_EER8bit, TIFF_COMPRESSION_EER7bit;

	int eer_upsampling;
	int nframes;
	int preread_start, preread_end;
	long long file_size;
	void readLegacy(FILE *fh);
	void lazyReadFrames();

	template <typename T>
	void render16K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons);

	template <typename T>
	void render8K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons);

	template <typename T>
	void render4K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons);

	static TIFFErrorHandler prevTIFFWarningHandler;

	public:

	EERRenderer();
	~EERRenderer();

	//TODO: Implement proper copy constructors. Currently, they are disabled to prevent memory corruption.
	EERRenderer(const EERRenderer&)
	{
		REPORT_ERROR(ERR_NOT_IMPLEMENTED, "Copy constructor for EERRenderer not implemented yet.");
	}

	EERRenderer& operator=(const EERRenderer&)
	{
		REPORT_ERROR(ERR_NOT_IMPLEMENTED, "Copy assignment operator for EERRenderer not implemented yet.");
	}

	// Wrapper to the default TIFF warning handler to suppress EER private tag warnings
	static void TIFFWarningHandler(const char* module, const char* fmt, va_list ap);
	static void silenceTIFFWarnings();

	// 1-indexed
	void setFramesOfInterest(int start, int end)
	{
		if (is_legacy)
			return;

		if (read_data)
			REPORT_ERROR(ERR_LOGIC_ERROR, "Logic error in EERRenderer::setFramesOfInterest(). This must be set before rendering.");
		preread_start = start - 1;
		preread_end = end - 1;
	}

	void read(FileName _fn_movie, int eer_upsampling=2);

	int getNFrames();
	int getWidth();
	int getHeight();

	// Frame indices are 1-indexed.
	// image is cleared.
	// This function is thread-safe (except for timing).
	// It is caller's responsibility to make sure type T does not overflow.
	template <typename T>
	long long renderFrames(int frame_start, int frame_end, MultidimArray<T> &image);

	// The gain reference for EER is not multiplicative! So the inverse is taken here.
	// 0 means defect.
	template <typename T>
	static void loadEERGain(FileName fn_gain, MultidimArray<T> &gain, int eer_upsampling)
	{
		const bool is_multiplicative = (fn_gain.getExtension() == "gain");
		if (is_multiplicative)
		{
			silenceTIFFWarnings();
			fn_gain += ":tif";
		}

		Image<T> original;
		original.read(fn_gain, true, 0, false, true); // explicitly use the first page
		const int nx_in = XSIZE(original());
		const int ny_in = YSIZE(original());
		const long long size_out = EER_IMAGE_WIDTH * eer_upsampling;

		// Revert Y flip in TIFF reader
		if (is_multiplicative)
		{
			const int ylim = ny_in / 2;
			for (int y1 = 0; y1 < ylim; y1++)
			{
				const int y2 = ny_in - 1 - y1;
				for (int x = 0; x < nx_in; x++)
				{
					const T tmp = DIRECT_A2D_ELEM(original(), y1, x);
					DIRECT_A2D_ELEM(original(), y1, x) = DIRECT_A2D_ELEM(original(), y2, x);
					DIRECT_A2D_ELEM(original(), y2, x) = tmp;
				}
			}
		} 

		if (eer_upsampling == 2 && nx_in == EER_IMAGE_WIDTH && ny_in == EER_IMAGE_HEIGHT) // gain = 4K and grid = 8K
		{
			gain.initZeros(size_out, size_out);
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(gain)
				DIRECT_A2D_ELEM(gain, i, j) = DIRECT_A2D_ELEM(original(), i / 2, j / 2);
		}
		else if ((eer_upsampling == 1 && nx_in == EER_IMAGE_WIDTH && ny_in == EER_IMAGE_HEIGHT) || // gain = 4K and grid = 4K
		         (eer_upsampling == 2 && nx_in == EER_IMAGE_WIDTH * 2 && ny_in == EER_IMAGE_HEIGHT * 2)) // gain = 8K and grid = 8K
		{
			gain = original();
		}
		else if (eer_upsampling == 1 && nx_in == EER_IMAGE_WIDTH * 2 && ny_in == EER_IMAGE_HEIGHT * 2) // gain = 8K and grid = 4K
		{
			gain.initZeros(size_out, size_out);
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(original())
				DIRECT_A2D_ELEM(gain, i / 2, j / 2) += DIRECT_A2D_ELEM(original(), i, j);
		}
		else
		{
			std::cerr << "Size of input gain: X = " << nx_in << " Y = " << ny_in << " Expected: X = " << size_out << " Y = " << size_out << std::endl;
			REPORT_ERROR(ERR_PARAM_INCORRECT, "Invalid gain size in EERRenderer::upsampleEERGain()");
		}
		
		if (!is_multiplicative)
		{
			double sum = 0;
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(gain)
				sum += DIRECT_MULTIDIM_ELEM(gain, n);
			sum /= size_out * size_out;

			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(gain)
			{
				if (DIRECT_MULTIDIM_ELEM(gain, n) != 0)
				{
					DIRECT_MULTIDIM_ELEM(gain, n) = sum / DIRECT_MULTIDIM_ELEM(gain, n);
				}
			}
		}
	}

	static bool isEER(FileName fn_movie)
	{
		FileName ext = fn_movie.getExtension();
		return (ext == "eer"  || ext == "ecc");
	}
};
