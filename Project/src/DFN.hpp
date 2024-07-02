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
    map<unsigned int, array<array<double, 3>,2>> TracesVertices = {};      //mappa Id-coordinate vertici
    map<unsigned int, array<unsigned int, 2>> TracesFractures = {};        //mappa Id-Id fratture generanti
    map<unsigned int, vector<unsigned int>> FractureTraces = {};           //mappa Id frattura-Id tracce
    map<array<unsigned int, 2>, bool> Tips = {};                               //mappa Id traccia/Id fratture-passante
    map<array<unsigned int, 2>, array<array<double,3>,2>> Retta = {};          //mappa Id traccia/Id fratture-Retta (P,V)
    map<array<unsigned int, 2>, array<unsigned int ,2>> LatiIntersecati = {};  //mappa Id traccia/Id fratture-Lati intersecati
};

struct Piano
{
    vector<unsigned int> PlaneId={};                                       //memorizza gli Id dei piani (coincidono con quelli delle fratture)
    map<unsigned int, array<double, 4>> Plane = {};                        //mappa Id-coefficienti del piano
};

struct PolygonalMesh
{
    unsigned int NumberCell0D = 0;                                // Numero di celle 0D
    vector<unsigned int> Cell0DId = {};                           // ID delle celle 0D
    vector<Vector2d> Cell0DCoordinates = {};                      // Coordinate delle celle 0D
    map<unsigned int, list<unsigned int>> Cell0DMarkers = {};     // Marcatori delle celle 0D

    unsigned int NumberCell1D = 0;                                // Numero di celle 1D
    vector<unsigned int> Cell1DId = {};                           // ID delle celle 1D
    vector<vector<unsigned int>> Cell1DVertices = {};             // Due indici per rappresentare l'origine e la fine dello spigolo
    map<unsigned int, list<unsigned int>> Cell1DMarkers = {};     // Marcatori delle celle 1D

    unsigned int NumberCell2D = 0;                                // Numero di celle 2D
    vector<unsigned int> Cell2DId = {};                           // ID delle celle 2D
    vector<vector<unsigned int>> Cell2DVertices = {};             // Indici dei vertici delle celle 2D
    vector<vector<unsigned int>> Cell2DEdges = {};                // Indici degli spigoli delle celle 2D
    map<unsigned int, list<unsigned int>> Cell2DMarkers = {};     // Marcatori delle celle 2D
};

}
