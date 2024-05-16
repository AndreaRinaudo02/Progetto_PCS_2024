#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include "Utils.hpp"

using namespace std;

namespace DFN_Library
{

struct DFN
{
    unsigned int NumberFractures = 0;                                      //numero di fratture
    vector<unsigned int> FractureId = {};                                  //Id fratture
    vector<vector<array<double, 3>>> FractureCoordinates = {};             //coordinate dei vertici (Id frattura, numero vertice, coordinata)
    map<unsigned int, vector<array<double, 3>>> FracturesVertices = {};    //dizionario Id-coordinate vertici

    unsigned int NumberTraces = 0;                                         //numero di tracce
    vector<unsigned int> TracesId = {};                                    //Id tracce
    vector<array<array<double, 3>,2>> TracesCoordinates = {};              //coordinate dei vertici (Id traccia, estremo, coordinata)
    map<unsigned int, array<array<double, 3>,2>> TracesVertices = {};      //mappa Id-coordinate vertici
    map<unsigned int, array<unsigned int, 2>> TracesFractures = {};        //mappa Id-Id fratture generanti
};

struct Piano
{
    vector<unsigned int> PlaneId={};                                       //memorizza gli Id dei piani (coincidono con quelli delle fratture)
    map<unsigned int, array<double, 4>> Plane = {};                        //mappa Id-coefficienti del piano
};

}
