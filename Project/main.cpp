#include "Utils.hpp"
#include "DFN.hpp"
#include <string>
#include <iostream>
#include <fstream>

using namespace std;
using namespace DFN_Library;

int main()
{
    DFN dfn;     //crea l'oggetto "dfn" che conterrà tutte le informazioni sul DFN


    string file_path = "DFN";      //indica la cartella contenente i file con i dati sui DFN


    if(!ImportDFN(file_path, dfn))
    {

        return 2;
    }

    return 0;
}
