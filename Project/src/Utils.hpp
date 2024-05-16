#pragma once

#include <string>
#include "DFN.hpp"

using namespace std;

namespace DFN_Library
{

bool ImportDFN(const string& file_path, DFN& dfn, Piano& Plane);

bool ImportFractures(const string& file_name, DFN& dfn, Piano& Plane);

void ParametriPiano(vector<array<double,3>> parametri, array<double, 4>& Coefficienti_piano);
}
