#pragma once

/*
zipstream Library License:
--------------------------

The zlib/libpng License Copyright (c) 2003 Jonathan de Halleux.

This software is provided 'as-is', without any express or implied warranty. In
no event will the authors be held liable for any damages arising from the use
of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim
   that you wrote the original software. If you use this software in a
   product, an acknowledgment in the product documentation would be
   appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

3. This notice may not be removed or altered from any source distribution

Author: Jonathan de Halleux, dehalleux@pelikhan.com, 2003

Altered by: Andreas Zieringer 2003 for OpenSG project
            made it platform independent, gzip conform, fixed gzip footer

Altered by: Geoffrey Hutchison 2005 for Open Babel project
            minor namespace modifications, VC++ compatibility

Altered by: Mathieu Malaterre 2008, for GDCM project
            when reading deflate bit stream in DICOM special handling of \0 is
            needed also when writing deflate back to disk, the add_footer must
            be called

Altered by: Timo Bingmann 2016-2019, make it actually work.
            Reversed DICOM weirdness I think.
*/

/*
Code is from:

http://www.codeproject.com/KB/stl/zipstream.aspx

 TODO:

zran.c
    index a zlib or gzip stream and randomly access it
    - illustrates the use of Z_BLOCK, inflatePrime(), and
      inflateSetDictionary() to provide random access

  So I might after all be able to implement seeking :)
*/

#include "zip_stream_fwd.hpp"

#include <zlib.h>

#include <istream>
#include <ostream>
#include <streambuf>
#include <vector>
#include <cstdint>

//! default gzip buffer size, change this to suite your needs
static const size_t zstream_default_buffer_size = 4096;

//! Compression strategy, see zlib doc.
enum class ZipStrategy { Filtered = 1, HuffmanOnly = 2, Default = 0 };

//! Header/Footer Formats
enum class ZipFormat { None, CrcFooter, GZip, Default = CrcFooter };

/******************************************************************************/
// basic_zip_streambuf

/**
 * \brief A stream decorator that takes raw input and zips it to a ostream.
 *
 * The class wraps up the inflate method of the zlib library 1.1.4
 * http://www.gzip.org/zlib/
 */
template <typename CharT, typename Traits = std::char_traits<CharT> >
class basic_zip_streambuf : public std::basic_streambuf<CharT, Traits> {
public:
    typedef std::basic_ostream<CharT, Traits>& ostream_reference;
    typedef unsigned char byte_type;
    typedef std::vector<byte_type> byte_vector_type;
    typedef std::vector<CharT> char_vector_type;
    typedef typename std::basic_streambuf<CharT, Traits>::int_type int_type;

    /**
     * Construct a zip stream More info on the following parameters can be found
     * in the zlib documentation.
     */
    basic_zip_streambuf(ostream_reference ostream, int level, ZipStrategy strategy, int window_size,
                        int memory_level, size_t buffer_size);

    ~basic_zip_streambuf();

    int sync() final;
    int_type overflow(int_type c) final;

    /** flushes the zip buffer and output buffer.
     *
     * This method should be called at the end of the compression. Calling flush
     * multiple times, will lower the compression ratio.
     */
    std::streamsize flush();

    //! returns a reference to the output stream
    ostream_reference get_ostream() const;

    //! returns the latest zlib error status
    int get_zerr() const;

    //! returns the crc of the input data compressed so far.
    uint32_t get_crc() const;

    //! returns the size (bytes) of the input data compressed so far.
    uint32_t get_in_size() const;

    //! returns the size (bytes) of the compressed data so far.
    unsigned long get_out_size() const;

private:
    bool zip_to_stream(char* buffer, std::streamsize buffer_size);

    ostream_reference ostream_;
    z_stream zip_stream_;
    int err_;
    byte_vector_type output_buffer_;
    char_vector_type buffer_;
    uint32_t crc_;
};

/******************************************************************************/
// basic_unzip_streambuf

/**
 * \brief A stream decorator that takes compressed input and unzips it to a
 * istream.
 *
 * The class wraps up the deflate method of the zlib library 1.1.4
 * http://www.gzip.org/zlib/
 */
template <typename CharT, typename Traits = std::char_traits<CharT> >
class basic_unzip_streambuf : public std::basic_streambuf<CharT, Traits> {
public:
    typedef std::basic_istream<CharT, Traits>& istream_reference;
    typedef unsigned char byte_type;
    typedef CharT char_type;
    typedef std::vector<byte_type> byte_vector_type;
    typedef std::vector<char_type> char_vector_type;
    typedef typename basic_unzip_streambuf<CharT, Traits>::int_type int_type;

    /** Construct a unzip stream. More info on the following parameters can be
     * found in the zlib documentation.
     */
    basic_unzip_streambuf(istream_reference istream, int window_size, size_t read_buffer_size,
                          size_t input_buffer_size);

    ~basic_unzip_streambuf();

    int_type underflow();

    //! returns the compressed input istream
    istream_reference get_istream();

    //! returns the zlib stream structure
    z_stream& get_zip_stream();

    //! returns the latest zlib error state
    int get_zerr() const;

    //! returns the crc of the uncompressed data so far
    uint32_t get_crc() const;

    //! returns the number of uncompressed bytes
    unsigned long get_out_size() const;

    //! returns the number of read compressed bytes
    uint32_t get_in_size() const;

private:
    void put_back_from_zip_stream();

    std::streamsize unzip_from_stream(char_type* buffer, std::streamsize buffer_size);

    size_t fill_input_buffer();

    istream_reference istream_;
    z_stream zip_stream_;
    int err_;
    byte_vector_type input_buffer_;
    char_vector_type buffer_;
    uint32_t crc_;
};

/******************************************************************************/
// basic_zip_ostream

template <typename CharT, typename Traits>
class basic_zip_ostream final
    : public basic_zip_streambuf<CharT, Traits>
    , public std::basic_ostream<CharT, Traits> {
public:
    typedef char char_type;
    typedef std::basic_ostream<CharT, Traits>& ostream_reference;
    typedef std::basic_ostream<CharT, Traits> ostream_type;

    explicit basic_zip_ostream(ostream_reference ostream, ZipFormat format = ZipFormat::Default,
                               int level = Z_BEST_COMPRESSION,
                               ZipStrategy strategy = ZipStrategy::Default,
                               /* windowBits is passed < 0 to suppress zlib header */
                               int window_size = -15, int memory_level = 8,
                               size_t buffer_size = zstream_default_buffer_size);

    ~basic_zip_ostream();

    //! returns type of format
    ZipFormat format() const;

    //! flush inner buffer and zipper buffer
    basic_zip_ostream<CharT, Traits>& zflush();

    //! flush buffer and write footer
    void finished();

private:
    basic_zip_ostream<CharT, Traits>& add_header();
    basic_zip_ostream<CharT, Traits>& add_footer();

    ZipFormat format_;
    bool added_footer_;
};

/******************************************************************************/
// basic_zip_istream

template <typename CharT, typename Traits>
class basic_zip_istream final
    : public basic_unzip_streambuf<CharT, Traits>
    , public std::basic_istream<CharT, Traits> {
public:
    typedef std::basic_istream<CharT, Traits>& istream_reference;
    typedef std::basic_istream<CharT, Traits> istream_type;

    explicit basic_zip_istream(istream_reference istream,
                               /* windowBits is passed < 0 to suppress zlib header */
                               int window_size = -15,
                               size_t read_buffer_size = zstream_default_buffer_size,
                               size_t input_buffer_size = zstream_default_buffer_size);

    //! returns true if it is a gzip file
    bool is_gzip() const;

    /*! return crc check result
     *
     * This must be called after the reading of compressed data is finished!
     * This method compares it to the crc of the uncompressed data.
     *
     *    \return true if crc check is succesful
     */
    bool check_crc();

    //! return data size check
    bool check_data_size() const;

    //! return the crc value in the file
    long get_gzip_crc() const;

    //! return the data size in the file
    long get_gzip_data_size() const;

    void read_footer();

protected:
    int check_header();

    bool is_gzip_;
    uint32_t gzip_crc_;
    uint32_t gzip_data_size_;
};

/******************************************************************************/
