/***************************************************************************
 *
 * Authors: Oier Lauzirika Zarrabeitia (olauzirika@cnb.csic.es)
 * 			Martin Salinas Anton (martin.salinas@cnb.csic.es)
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

#include <algorithm>
#include <string>
#include <tiffio.h>
#include <vector>

#include "multidim_array.h"
#include "xmipp_error.h"
#include "xmipp_filename.h"
#include "xmipp_image_base.h"
#include "metadata_vec.h"

#ifdef TIMING
	#define RCTIC(label) (EERtimer.tic(label))
	#define RCTOC(label) (EERtimer.toc(label))

	Timer EERtimer;
	int TIMING_READ_EER = EERtimer.setNew("read EER");
	int TIMING_BUILD_INDEX = EERtimer.setNew("build index");
	int TIMING_UNPACK_RLE = EERtimer.setNew("unpack RLE");
	int TIMING_RENDER_ELECTRONS = EERtimer.setNew("render electrons");
#else
	#define RCTIC(label)
	#define RCTOC(label)
#endif

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

	void readLegacy(FILE *fh)
	{
		/* Load everything first */
		RCTIC(TIMING_READ_EER);
		buf = (unsigned char*)malloc(file_size);
		if (buf == NULL)
			REPORT_ERROR(ERR_MEM_NOTENOUGH, "Failed to allocate the buffer.");
		if (fread(buf, sizeof(char), file_size, fh) != file_size)
			REPORT_ERROR(ERR_IO_SIZE, "EERRenderer::readLegacy: Failed to read the expected size from " + fn_movie);
		RCTOC(TIMING_READ_EER);

		/* Build frame index */
		RCTIC(TIMING_BUILD_INDEX);
		long long pos = file_size;
		while (pos > 0)
		{
			pos -= EER_LEN_FOOTER;
			if (strncmp(EER_FOOTER_OK, (char*)buf + pos, EER_LEN_FOOTER) == 0)
			{
				pos -= 8;
				long long frame_size = *(long long*)(buf + pos) * sizeof(long long);
				pos -= frame_size;
				frame_starts.push_back(pos);
				frame_sizes.push_back(frame_size);
	#ifdef DEBUG_EER
				printf("Frame: LAST-%5d Start at: %08d Frame size: %lld\n", frame_starts.size(), pos, frame_size);
	#endif
			}
			else // if (strncmp(EER_FOOTER_ERR, (char*)buf + pos, EER_LEN_FOOTER) == 0)
			{
				REPORT_ERROR(ERR_DOCFILE, "Broken frame in file " + fn_movie);
			}
		}
		std::reverse(frame_starts.begin(), frame_starts.end());
		std::reverse(frame_sizes.begin(), frame_sizes.end());
		RCTOC(TIMING_BUILD_INDEX);

		nframes = frame_starts.size();
		read_data = true;
	}

	void lazyReadFrames()
	{
		#pragma omp critical(EERRenderer_lazyReadFrames)
		{
			if (!read_data) // cannot return from within omp critical
			{	
				TIFF *ftiff = TIFFOpen(fn_movie.c_str(), "r");

				frame_starts.resize(nframes, 0);
				frame_sizes.resize(nframes, 0);
				buf = (unsigned char*)malloc(file_size); // This is big enough
				if (buf == NULL)
					REPORT_ERROR(ERR_MEM_NOTENOUGH, "Failed to allocate the buffer for " + fn_movie);
				long long pos = 0;

				// Read everything
				for (int frame = 0; frame < nframes; frame++)
				{
					if ((preread_start > 0 && frame < preread_start) ||
							(preread_end > 0 && frame > preread_end))
						continue;

					TIFFSetDirectory(ftiff, frame);
					const int nstrips = TIFFNumberOfStrips(ftiff);
					frame_starts[frame] = pos;

					for (int strip = 0; strip < nstrips; strip++)
					{
						const int strip_size = TIFFRawStripSize(ftiff, strip);
						if (pos + strip_size >= file_size)
							REPORT_ERROR(ERR_LOGIC_ERROR, "EER: buffer overflow when reading raw strips.");

						TIFFReadRawStrip(ftiff, strip, buf + pos, strip_size);
						pos += strip_size;
						frame_sizes[frame] += strip_size;
					}
		#ifdef DEBUG_EER
					printf("EER in TIFF: Read frame %d from %s, nstrips = %d, current pos in buffer = %lld / %lld\n", frame, fn_movie.c_str(), nstrips, pos, file_size);
		#endif
				}

				TIFFClose(ftiff);

				read_data = true;
			}
		}
	}

	template <typename T>
	void render16K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons)
	{
		for (int i = 0; i < n_electrons; i++)
		{
			int x = ((positions[i] & 4095) << 2) | (symbols[i] & 3); // 4095 = 111111111111b, 3 = 00000011b
			int y = ((positions[i] >> 12) << 2) | ((symbols[i] & 12) >> 2); //  4096 = 2^12, 12 = 00001100b
				DIRECT_A2D_ELEM(image, y, x)++;
		}
	}

	template <typename T>
	void render8K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons)
	{
		for (int i = 0; i < n_electrons; i++)
		{
			int x = ((positions[i] & 4095) << 1) | ((symbols[i] & 2) >> 1); // 4095 = 111111111111b, 2 = 00000010b
			int y = ((positions[i] >> 12) << 1) | ((symbols[i] & 8) >> 3); //  4096 = 2^12, 8 = 00001000b
				DIRECT_A2D_ELEM(image, y, x)++;
		}
	}

	template <typename T>
	void render4K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons)
	{
		for (int i = 0; i < n_electrons; i++)
		{
			int x = positions[i] & 4095; // 4095 = 111111111111b
			int y = positions[i] >> 12; //  4096 = 2^12
				DIRECT_A2D_ELEM(image, y, x)++;
		}
	}

	static TIFFErrorHandler prevTIFFWarningHandler;

	public:

	EERRenderer()
	{
		ready = false;
		read_data = false;
		buf = NULL;
		preread_start = -1;
		preread_end = -1;
		eer_upsampling = 2;
	}

	~EERRenderer()
	{
		if (buf != NULL)
			free(buf);
	}

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
	static void TIFFWarningHandler(const char* module, const char* fmt, va_list ap)
	{
		// Silence warnings for private tags
		if (strcmp("Unknown field with tag %d (0x%x) encountered", fmt) == 0)
			return;

		if (prevTIFFWarningHandler != NULL)
			prevTIFFWarningHandler(module, fmt, ap);
	}

	static void silenceTIFFWarnings()
	{
		if (prevTIFFWarningHandler == NULL)
		{
			// Thread safety issue:
			// Calling this simultaneously is safe but
			TIFFErrorHandler prev = TIFFSetWarningHandler(EERRenderer::TIFFWarningHandler);

			// we have to make sure prevTIFFWarningHandler does NOT become our own TIFFWarningHandler
			// to avoid an infinite loop.
			if (prev != EERRenderer::TIFFWarningHandler)
				prevTIFFWarningHandler = prev;
		}
	}

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

	void read(const FileName &_fn_movie, int eer_upsampling)
	{
		if (ready)
			REPORT_ERROR(ERR_LOGIC_ERROR, "Logic error: you cannot recycle EERRenderer for multiple files (now)");

		if (eer_upsampling == 1 || eer_upsampling == 2 || eer_upsampling == 3)
			this->eer_upsampling = eer_upsampling;
		else
		{
			std::cerr << "EERRenderer::read: eer_upsampling = " << eer_upsampling << std::endl;
			REPORT_ERROR(ERR_PARAM_INCORRECT, "EERRenderer::read: eer_upsampling must be 1, 2 or 3.");
		}

		fn_movie = _fn_movie;

		// First of all, check the file size
		FILE *fh = fopen(fn_movie.c_str(), "r");
		if (fh == NULL)
			REPORT_ERROR(ERR_IO_NOTOPEN, "Failed to open " + fn_movie);

		fseek(fh, 0, SEEK_END);
		file_size = ftell(fh);
		fseek(fh, 0, SEEK_SET);

		silenceTIFFWarnings();

		// Try reading as TIFF
		TIFF *ftiff = TIFFOpen(fn_movie.c_str(), "r");

		if (ftiff == NULL)
		{
			is_legacy = true;
			is_7bit = false;
			readLegacy(fh);
		}
		else
		{
			is_legacy = false;

			// Check width & size
			int width, height;
			uint16_t compression = 0;
			TIFFGetField(ftiff, TIFFTAG_IMAGEWIDTH, &width);
			TIFFGetField(ftiff, TIFFTAG_IMAGELENGTH, &height);
			TIFFGetField(ftiff, TIFFTAG_COMPRESSION, &compression);

	#ifdef DEBUG_EER
			printf("EER in TIFF: %s size = %ld, width = %d, height = %d, compression = %d\n", fn_movie.c_str(), file_size, width, height, compression);
	#endif

			// TIA can write an EER file whose first page is a sum and compressoin == 1.
			// This is not supported (yet). EPU never writes such movies.
			if (compression == EERRenderer::TIFF_COMPRESSION_EER8bit)
				is_7bit = false;
			else if (compression == EERRenderer::TIFF_COMPRESSION_EER7bit)
				is_7bit = true;
			else
				REPORT_ERROR(ERR_PARAM_INCORRECT, "Unknown compression scheme for EER " + integerToString(compression));

			if (width != EER_IMAGE_WIDTH || height != EER_IMAGE_HEIGHT)
				REPORT_ERROR(ERR_PARAM_INCORRECT, "Currently we support only 4096x4096 pixel EER movies.");

			// Find the number of frames
			nframes = TIFFNumberOfDirectories(ftiff);
			TIFFClose(ftiff);
	#ifdef DEBUG_EER
			printf("EER in TIFF: %s nframes = %d\n", fn_movie.c_str(), nframes);
	#endif
		}

		fclose(fh);
		ready = true;
	}

	int getNFrames()
	{
		if (!ready)
			REPORT_ERROR(ERR_LOGIC_ERROR, "EERRenderer::getNFrames called before ready.");

		return nframes;
	}

	int getWidth()
	{
		if (!ready)
			REPORT_ERROR(ERR_LOGIC_ERROR, "EERRenderer::getNFrames called before ready.");

		return EER_IMAGE_WIDTH << (eer_upsampling - 1);
	}

	int getHeight()
	{
		if (!ready)
			REPORT_ERROR(ERR_LOGIC_ERROR, "EERRenderer::getNFrames called before ready.");

		return EER_IMAGE_HEIGHT << (eer_upsampling - 1);
	}

	// Frame indices are 1-indexed.
	// image is cleared.
	// This function is thread-safe (except for timing).
	// It is caller's responsibility to make sure type T does not overflow.
	template <typename T>
	long long renderFrames(int frame_start, int frame_end, MultidimArray<T> &image)
	{
		if (!ready)
			REPORT_ERROR(ERR_LOGIC_ERROR, "EERRenderer::renderNFrames called before ready.");

		lazyReadFrames();

		if (frame_start < 0 || frame_start >= getNFrames() ||
				frame_end < frame_start || frame_end >= getNFrames())
		{
			std::cerr << "EERRenderer::renderFrames(frame_start = " << frame_start << ", frame_end = " << frame_end << "),  NFrames = " << getNFrames() << std::endl;
			REPORT_ERROR(ERR_PARAM_INCORRECT, "Invalid frame range was requested.");
		}

		long long total_n_electron = 0;

		std::vector<unsigned int> positions;
		std::vector<unsigned char> symbols;
		image.initZeros(getHeight(), getWidth());

		for (int iframe = frame_start; iframe <= frame_end; iframe++)
		{
			RCTIC(TIMING_UNPACK_RLE);
			if ((preread_start > 0 && iframe < preread_start) ||
						(preread_end > 0 && iframe > preread_end))
			{
				std::cerr << "EERRenderer::renderFrames(frame_start = " << frame_start + 1 << ", frame_end = " << frame_end + 1<< "),  NFrames = " << getNFrames() << " preread_start = " << preread_start + 1 << " prered_end = " << preread_end + 1<< std::endl;
				REPORT_ERROR(ERR_LOGIC_ERROR, "Tried to render frames outside pre-read region");
			}		

			long long pos = frame_starts[iframe];
			unsigned int n_pix = 0, n_electron = 0;
			const int max_electrons = frame_sizes[iframe] * 2; // at 4 bits per electron (very permissive bound!)
			if (positions.size() < max_electrons)
			{
				positions.resize(max_electrons);
				symbols.resize(max_electrons);
			}

			if (is_7bit)
			{
				unsigned int bit_pos = 0; // 4 K * 4 K * 11 bit << 2 ** 32
				unsigned char p, s;

				while (true)
				{
					// Fetch 32 bits and unpack up to 2 chunks of 7 + 4 bits.
					// This is faster than unpack 7 and 4 bits sequentially.
					// Since the size of buf is larger than the actual size by the TIFF header size,
					// it is always safe to read ahead.

					long long first_byte = pos + (bit_pos >> 3);
					const unsigned int bit_offset_in_first_byte = bit_pos & 7; // 7 = 00000111 (same as % 8)
					const unsigned int chunk = (*(unsigned int*)(buf + first_byte)) >> bit_offset_in_first_byte;

					p = (unsigned char)(chunk & 127); // 127 = 01111111
					bit_pos += 7; // TODO: we can remove this for further speed.
					n_pix += p;
					if (n_pix >= EER_IMAGE_PIXELS) break;
					if (p == 127) continue; // this should be rare.
					
					s = (unsigned char)((chunk >> 7) & 15) ^ 0x0A; // 15 = 00001111; See below for 0x0A
					bit_pos += 4;
					positions[n_electron] = n_pix;
					symbols[n_electron] = s;
					n_electron++;
					n_pix++;

					p = (unsigned char)((chunk >> 11) & 127); // 127 = 01111111
					bit_pos += 7;
					n_pix += p;
					if (n_pix >= EER_IMAGE_PIXELS) break;
					if (p == 127) continue;
					
					s = (unsigned char)((chunk >> 18) & 15) ^ 0x0A; // 15 = 00001111; See below for 0x0A
					bit_pos += 4;
					positions[n_electron] = n_pix;
					symbols[n_electron] = s;
					n_electron++;
					n_pix++;
				}
			}
			else
			{
				// unpack every two symbols = 12 bit * 2 = 24 bit = 3 byte
				// high <- |bbbbBBBB|BBBBaaaa|AAAAAAAA| -> low
				// With SIMD intrinsics at the SSSE3 level, we can unpack 10 symbols (120 bits) simultaneously.
				unsigned char p1, p2, s1, s2;

				const long long pos_limit = frame_starts[iframe] + frame_sizes[iframe];
				// Because there is a footer, it is safe to go beyond the limit by two bytes.
				while (pos < pos_limit)
				{
					// symbol is bit tricky. 0000YyXx; Y and X must be flipped.
					p1 = buf[pos];
					s1 = (buf[pos + 1] & 0x0F) ^ 0x0A; // 0x0F = 00001111, 0x0A = 00001010

					p2 = (buf[pos + 1] >> 4) | (buf[pos + 2] << 4);
					s2 = (buf[pos + 2] >> 4) ^ 0x0A;

					// Note the order. Add p before checking the size and placing a new electron.
					n_pix += p1;
					if (n_pix >= EER_IMAGE_PIXELS) break;
					if (p1 < 255)
					{
						positions[n_electron] = n_pix;
						symbols[n_electron] = s1;
						n_electron++;
						n_pix++;
					}

					n_pix += p2;
					if (n_pix >= EER_IMAGE_PIXELS) break;
					if (p2 < 255)
					{
						positions[n_electron] = n_pix;
						symbols[n_electron] = s2;
						n_electron++;
						n_pix++;
					}
	#ifdef DEBUG_EER_DETAIL
					printf("%d: %u %u, %u %u %d\n", pos, p1, s1, p2, s2, n_pix);
	#endif
					pos += 3;
				}
			}

			if (n_pix != EER_IMAGE_PIXELS)
			{
				std::cerr << "WARNING: The number of pixels is not right in " + fn_movie + " frame " + integerToString(iframe + 1) + ". Probably this frame is corrupted. This frame is skipped." << std::endl;
				continue;
			}

			RCTOC(TIMING_UNPACK_RLE);

			RCTIC(TIMING_RENDER_ELECTRONS);
			if (eer_upsampling == 3)
				render16K(image, positions, symbols, n_electron);
			else if (eer_upsampling == 2)
				render8K(image, positions, symbols, n_electron);
			else if (eer_upsampling == 1)
				render4K(image, positions, symbols, n_electron);
			else
				REPORT_ERROR(ERR_LOGIC_ERROR, "Invalid EER upsample");
			RCTOC(TIMING_RENDER_ELECTRONS);

			total_n_electron += n_electron;
	#ifdef DEBUG_EER
			printf("Decoded %lld electrons / %d pixels from frame %5d.\n", n_electron, n_pix, iframe);
	#endif
		}
	#ifdef DEBUG_EER
		printf("Decoded %lld electrons in total.\n", total_n_electron);
	#endif

	#ifdef TIMING
		EERtimer.printTimes(false);
	#endif

		return total_n_electron;
	}

	static bool isEER(const FileName &fn_movie)
	{
		FileName ext = fn_movie.getExtension();
		return (ext == "eer"  || ext == "ecc");
	}
};

const char EERRenderer::EER_FOOTER_OK[]  = "ThermoFisherECComprOK000";
const char EERRenderer::EER_FOOTER_ERR[] = "ThermoFisherECComprERR00";
const int EERRenderer::EER_IMAGE_WIDTH = 4096;
const int EERRenderer::EER_IMAGE_HEIGHT = 4096;
const int EERRenderer::EER_IMAGE_PIXELS = EERRenderer::EER_IMAGE_WIDTH * EERRenderer::EER_IMAGE_HEIGHT;
const unsigned int EERRenderer::EER_LEN_FOOTER = 24;
const uint16_t EERRenderer::TIFF_COMPRESSION_EER8bit = 65000;
const uint16_t EERRenderer::TIFF_COMPRESSION_EER7bit = 65001;

TIFFErrorHandler EERRenderer::prevTIFFWarningHandler = NULL;

int ImageBase::readEER(size_t select_img) {
	int upsampling = 1;
	StringVector info;
	DataType datatype;
	size_t found = filename.find_first_of("#");
	FileName infolist = filename.substr(found + 1);
	filename = filename.substr(0, found);
	splitString(infolist, ",", info, false);

	if (info.size() < 3)
			REPORT_ERROR(ERR_ARG_MISSING, (String) "Cannot open file " + filename +
										". Not enough header arguments.");

	const int fractioning = std::atoi(info[0].c_str());
	if (fractioning < 1)
	{
		REPORT_ERROR(ERR_PARAM_INCORRECT, "Incorrect fractioning value (must be greater than zero)");
	}
	if (select_img > fractioning) {
		REPORT_ERROR(ERR_PARAM_INCORRECT, (String) "Incorrect frame number selected (" + std::to_string(select_img) + "). Number of frames is " + std::to_string(fractioning) + ".");
	}

	int _xDim,_yDim,_zDim;
	size_t _nDim;
  	if (info[1] == "4K")
	{
		_xDim = _yDim = 4096;
	} 
	else if (info[1] == "8K")
	{
		_xDim = _yDim = 8192;
	}
	else {
		REPORT_ERROR(ERR_PARAM_INCORRECT, "Incorrect output size. Valid sizes are: 4K, 8K.");
	}

	_zDim = 1;
	_nDim = select_img > 0 ? 1 : fractioning;
	setDimensions(_xDim, _yDim, _zDim, _nDim);
	mdaBase->coreAllocateReuse();

  	if(info[2] == "uint8")
		datatype = DT_UChar;
	else if (info[2] == "uint16")
		datatype = DT_UShort;
	else
		REPORT_ERROR(ERR_PARAM_INCORRECT, "Incorrect output data type. Valid types are: uint8, uint16.");
	
	MDMainHeader.setValue(MDL_SAMPLINGRATE_X,(double) -1);
	MDMainHeader.setValue(MDL_SAMPLINGRATE_Y,(double) -1);
	MDMainHeader.setValue(MDL_DATATYPE,(int)datatype);

	MD.clear();
	MD.push_back(std::unique_ptr<MDRowVec>(new MDRowVec(MDL::emptyHeaderVec())));

	if( dataMode < DATA )
		return 0;

	EERRenderer renderer;
	renderer.read(dataFName, upsampling);

	MultidimArray<int> buffer;
	const auto nEerFrames = renderer.getNFrames();
	const auto step = nEerFrames / fractioning;
	if (select_img > 0)
	{
		// Render single frame
		if (select_img > fractioning)
		{
			REPORT_ERROR(ERR_LOGIC_ERROR, "Requested frame greater than the fractioning");
		}
		else
		{
			
			buffer.resizeNoCopy(_yDim, _xDim);
			const auto first = (select_img-1)*step;
			const auto last = first + step - 1;
			renderer.renderFrames(first, last, buffer);
		}
	}
	else
	{
		// Render the whole movie
		buffer.resizeNoCopy(fractioning, 1, _yDim, _xDim);
		MultidimArray<int> frameAlias;
		for(size_t i = 0; i < fractioning; ++i)
		{
			frameAlias.aliasImageInStack(buffer, i);
			const auto first = i*step;
			const auto last = first + step - 1;
			renderer.renderFrames(first, last, frameAlias);
		}
	}
	
	setPage2T(
		0UL, reinterpret_cast<char*>(MULTIDIM_ARRAY(buffer)),
		DT_Int,
		MULTIDIM_SIZE(buffer)
	);

	return 0;
}
