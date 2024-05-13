#pragma once

#include <string>
#include "DFN.hpp"

using namespace std;

namespace DFN_Library
{

bool ImportDFN(const string& file_path, DFN& dfn);

bool ImportFractures(const string& file_name, DFN& dfn);

}
