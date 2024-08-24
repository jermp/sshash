#pragma once

#include <string>

template <class CharT, class Traits = std::char_traits<CharT> >
class basic_zip_ostream;

template <class CharT, class Traits = std::char_traits<CharT> >
class basic_zip_istream;

//! A typedef for basic_zip_ostream<char>
using zip_ostream = basic_zip_ostream<char>;
//! A typedef for basic_zip_istream<char>
using zip_istream = basic_zip_istream<char>;
