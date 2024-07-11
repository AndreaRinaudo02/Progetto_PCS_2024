#pragma once

#include <string>
#include "DFN.hpp"

using namespace std;

namespace DFN_Library
{

bool ImportDFN(const string& file_path, DFN& dfn, Piano& Plane);

bool ImportFractures(const string& file_name, DFN& dfn, Piano& Plane);

void ParametriPiano(vector<array<double,3>> parametri, array<double, 4>& Coefficienti_piano);

void Calcola_tracce(DFN& dfn, Piano& Plane);

void Fratture_vicine(DFN& dfn, vector<array<unsigned int, 2>>& coppie_vicine);

void Crea_bolle(DFN& dfn, vector<array<double , 4>>& bolle);

void StampaTracce(const string& file_name, DFN& dfn);

void RettaIntersezione(Piano &Plane,
                       vector<array<unsigned int, 2>>& coppie_vicine,
                       map<array<unsigned int, 2>, array<array<double, 3>, 2>>& Retta);

void IntersezioneLati(map<array<unsigned int, 2>,array<array<double, 3>, 2>>& Retta, DFN& dfn);

void TagliaTracce(DFN& dfn, PolygonalMesh& mesh);

void StampaSottopoligoni(const string& file_name, DFN& dfn, PolygonalMesh& mesh);

}
