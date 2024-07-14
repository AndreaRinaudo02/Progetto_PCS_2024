#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <list>
#include <Eigen/Eigen>
#include <map>
#include "Utils.hpp"

using namespace std;
using namespace Eigen;

namespace DFN_Library
{

struct DFN
{
    unsigned int NumberFractures = 0;                                      //numero di fratture
    vector<unsigned int> FractureId = {};                                  //Id fratture
    vector<vector<array<double, 3>>> FractureCoordinates = {};             //coordinate dei vertici (Id frattura, numero vertice, coordinata)
    vector<vector<array<double, 3>>> FracturesVertices = {};               //mappa Id-coordinate vertici

    unsigned int NumberTraces = 0;                                         //numero di tracce
    vector<unsigned int> TracesId = {};                                    //Id tracce
    vector<array<array<double, 3>,2>> TracesVertices = {};                 //mappa Id-coordinate vertici
    vector<array<unsigned int, 2>> TracesFractures = {};                   //mappa Id-Id fratture generanti
    vector<vector<unsigned int>> FractureTraces = {};                      //mappa Id frattura-Id tracce
    map<array<unsigned int, 2>, bool> Tips = {};                               //mappa Id traccia/Id fratture-passante
    map<array<unsigned int, 2>, array<array<double,3>,2>> Retta = {};          //mappa Id traccia/Id fratture-Retta (P,V)
    map<array<unsigned int, 2>, vector<unsigned int>> LatiIntersecati = {};  //mappa Id traccia/Id fratture-Lati intersecati

    vector<vector<vector<array<double,3>>>> Sottopoligoni = {};                //mappa Id frattura/coordinate sottopoligoni
};

struct Piano
{
    vector<unsigned int> PlaneId={};                                       //memorizza gli Id dei piani (coincidono con quelli delle fratture)
    vector<array<double, 4>> Plane = {};                                   //mappa Id-coefficienti del piano
};

struct PolygonalMesh
{
    unsigned int NumberCell0D = {};                               // Numero di celle 0D
    vector<unsigned int> Cell0DId = {};                           // ID delle celle 0D
    vector<array<double, 3>> Cell0DCoordinates = {};              // Coordinate delle celle 0D

    unsigned int NumberCell1D = {};                               // Numero di celle 1D
    vector<unsigned int> Cell1DId = {};                           // ID delle celle 1D
    vector<array<unsigned int, 2>> Cell1DVertices = {};           // Due indici per rappresentare l'origine e la fine dello spigolo

    unsigned int NumberCell2D = {};                               // Numero di celle 2D
    vector<unsigned int> Cell2DId = {};                           // ID delle celle 2D
    vector<vector<unsigned int>> Cell2DVertices = {};             // Indici dei vertici delle celle 2D
    vector<vector<unsigned int>> Cell2DEdges = {};                // Indici degli spigoli delle celle 2D
};

}
