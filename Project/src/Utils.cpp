#include <string>
#include <fstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Utils.hpp"
#include "DFN.hpp"

using namespace std;
namespace DFN_Library

{

bool ImportDFN(const string& file_path, DFN& dfn, Piano& Plane)
{
    if(!ImportFractures(file_path + "/FR3_data.txt", dfn, Plane))
    {
        return false;
    }
    return true;
}

bool ImportFractures(const string &file_name, DFN& dfn, Piano &Plane)

{
    ifstream file;
    file.open(file_name);     //apre il file

    if(!file.is_open())
    {
        cerr << "Error: Unable to open file " << file_name << endl;
        return false;
    }

    unsigned int n;
    unsigned int m;
    char ch;
    double coordinata;
    string riga;

    getline(file, riga);       //scarta intestazione

    file >> m;                 //numero di fratture

    dfn.NumberFractures = m;   //aggiorna numero di fratture

    dfn.FractureCoordinates.resize(m);     //riscala i vettori
    dfn.FractureId.resize(m);
    Plane.PlaneId.resize(m);

    getline(file, riga);       //scarta riga

    for (unsigned int i=0; i < m; i++)   //ciclo for (per iterare ogni frattura)
    {

        getline(file, riga);    //scarta intestazione
        getline(file, riga);    //legge Id e numero di vertici
        stringstream ss(riga);

        ss >> n >> ch;                     //legge Id frattura, scarta ";"
        dfn.FractureId[i] = n;            //aggiunge Id all'elenco
        dfn.FracturesVertices[n] = {};     //aggiunge Id come chiave della mappa
        Plane.PlaneId[i]=n;

        ss >> n;      //numero di vertici della frattura

        dfn.FractureCoordinates[i].resize(n);     //riscala il vettore

        getline(file, riga);    //scarta intestazione

        for (int k=0; k < 3; k++)   //ciclo for (per ogni componente dei vertici)
        {
            getline(file, riga);    // riga con le coordinate
            stringstream ss(riga);

            for (unsigned int j =0; j < n; j++)                // ciclo for (per ogni vertice)
            {
                ss >> coordinata >> ch;                            //estrae coordinata, scarta ";"
                dfn.FractureCoordinates[i][j][k] = coordinata;     //memorizza la coordinata
            }

        }

        dfn.FracturesVertices[dfn.FractureId[i]] = dfn.FractureCoordinates[i];    //chiave: Id, valore: vettore di array con le coordinate del vertice

    }


    file.close();     //chiude il file

    return true;
}

double* ParametriPiano(vector<array<double,3>> parametri, Piano& plane)
{
    double a, b, c, d;
    double* array = new double[4];
    a = ((parametri[1][1]-parametri[0][1])*(parametri[1][2]-parametri[0][2])) - ((parametri[2][1]-parametri[0][1])*(parametri[2][2]-parametri[0][2]));
    b = ((parametri[1][0]-parametri[0][0])*(parametri[2][2]-parametri[0][2])) - ((parametri[1][2]-parametri[0][2])*(parametri[2][0]-parametri[0][0]));
    c = ((parametri[1][0]-parametri[0][0])*(parametri[2][1]-parametri[0][1])) - ((parametri[1][1]-parametri[0][1])*(parametri[2][0]-parametri[0][0]));
    d= -(a*parametri[0][0] + b*parametri[0][1] + c*parametri[0][2]);
    array[0]=a;
    array[1]=b;
    array[2]=c;
    array[3]=d;
    return array;

}
}