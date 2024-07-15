#include "Utils.hpp"
#include "DFN.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <gtest/gtest.h>

using namespace std;
using namespace DFN_Library;

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    int n = RUN_ALL_TESTS();

    if(n > 0)
    {
        cerr << "Test non superati" << endl;
        return 1;
    }

    DFN dfn;                             //crea l'oggetto "dfn" che conterrà tutte le informazioni sul DFN
    Piano piano;                         //crea l'oggetto "piano" che conterrà le informazioni utili a calcolare le fratture
    vector<PolygonalMesh> Mesh = {};     //crea un vettore di oggetti "mesh" che conterrà le informazioni sulla mesh delle fratture

    string file_path = "DFN";            //indica la cartella contenente i file con i dati sui DFN


    if(!ImportDFN(file_path, dfn, piano))
    {

        return 2;
    }

    Mesh.resize(dfn.NumberFractures);

    TagliaTracce(dfn, Mesh);

    return 0;
}
