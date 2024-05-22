#include <string>
#include <fstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
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

bool ImportFractures(const string &file_name, DFN& dfn, Piano& Plane)

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

        array<double,4> Coefficienti_piano = {};
        ParametriPiano(dfn.FractureCoordinates[i], Coefficienti_piano);           //calcola i coefficienti del piano contenente la frattura
        Plane.Plane[dfn.FractureId[i]] = Coefficienti_piano;                      //chiave: Id, valore: array (dimensione 4) con i parametri del piano in ordine (a, b, c, d)
    }


    file.close();     //chiude il file

    Calcola_tracce(dfn, Plane);

    return true;
}

void ParametriPiano(vector<array<double,3>> parametri, array<double, 4>& Coefficienti_piano)
{
    double a, b, c, d;

    a = ((parametri[1][1] - parametri[0][1]) * (parametri[1][2] - parametri[0][2])) - ((parametri[2][1] - parametri[0][1]) * (parametri[2][2] - parametri[0][2]));
    b = ((parametri[1][0] - parametri[0][0]) * (parametri[2][2] - parametri[0][2])) - ((parametri[1][2] - parametri[0][2]) * (parametri[2][0] - parametri[0][0]));
    c = ((parametri[1][0] - parametri[0][0]) * (parametri[2][1] - parametri[0][1])) - ((parametri[1][1] - parametri[0][1]) * (parametri[2][0] - parametri[0][0]));
    d= -(a * parametri[0][0] + b * parametri[0][1] + c * parametri[0][2]);

    Coefficienti_piano[0] = a;
    Coefficienti_piano[1] = b;
    Coefficienti_piano[2] = c;
    Coefficienti_piano[3] = d;
}

void Calcola_tracce(DFN& dfn, Piano& piano)
{
    /*
    1 calcolare fratture vicine (opzionale)   OK
    2 intersezione tra i due piani (direzione retta di intersezione) (v)
    3 intersezione con terzo piano perp. a retta di intersezione (punto iniziale della retta) (P + tv)
    4 calcolo rette lati delle fratture
    5 intersezione tra rette lati e retta di intersezione
    6 calcolo estremi frattura
    */

    //1
    vector<array<unsigned int, 2>> coppie_vicine = {};
    Fratture_vicine(dfn, coppie_vicine);

    for (int i = 0; i < coppie_vicine.size(); i++)
    {
        cout << "coppie: " << coppie_vicine[i][0]<< " " <<coppie_vicine[i][1] << endl;
    }
}

void Fratture_vicine(DFN& dfn, vector<array<unsigned int, 2>>& coppie_vicine)
{
    for(unsigned int i = 0; i < dfn.NumberFractures - 1; i++)
    {
        vector<array<double, 3>> coordinate_1 = dfn.FracturesVertices[i];
        array<double, 3> centro_1 = {};

        double x = 0;
        double y = 0;
        double z = 0;
        for(int a = 0; a < coordinate_1.size(); a++)
        {

            x += coordinate_1[a][0];
            y += coordinate_1[a][1];
            z += coordinate_1[a][2];
        }

        centro_1[0] = x/coordinate_1.size();
        centro_1[1] = y/coordinate_1.size();
        centro_1[2] = z/coordinate_1.size();

        vector<double> distanze = {};

        for(int c = 0; c < coordinate_1.size(); c++)
        {
            x = centro_1[0] - coordinate_1[c][0];
            y = centro_1[1] - coordinate_1[c][1];
            z = centro_1[2] - coordinate_1[c][2];
            x = pow(x,2);
            y = pow(y,2);
            z = pow(z,2);

            distanze.push_back(x/4 + y/4 + z/4);  //approssima la radice quadrata
        }

        double distanza_1 = *(max_element(distanze.begin(), distanze.end()));

        for(unsigned int j = i+1; j < dfn.NumberFractures; j++)
        {

            distanze.clear();
            vector<array<double, 3>> coordinate_2 = dfn.FracturesVertices[j];
            array<double, 3> centro_2 = {};


            x = 0;
            y = 0;
            z = 0;

            for(int b = 0; b < coordinate_2.size(); b++)
            {

                x += coordinate_2[b][0];
                y += coordinate_2[b][1];
                z += coordinate_2[b][2];
            }

            centro_2[0] = x/coordinate_2.size();
            centro_2[1] = y/coordinate_2.size();
            centro_2[2] = z/coordinate_2.size();



            for(int d = 0; d < coordinate_2.size(); d++)
            {
                x = centro_2[0] - coordinate_2[d][0];
                y = centro_2[1] - coordinate_2[d][1];
                z = centro_2[2] - coordinate_2[d][2];
                x = pow(x,2);
                y = pow(y,2);
                z = pow(z,2);

                distanze.push_back(x/4 + y/4 + z/4);  //approssima la radice quadrata
            }

            double distanza_2 = *(max_element(distanze.begin(), distanze.end()));

            x = centro_1[0] - centro_2[0];
            y = centro_1[1] - centro_2[1];
            z = centro_1[2] - centro_2[2];
            x = pow(x,2);
            y = pow(y,2);
            z = pow(z,2);

            double distanza_centri = x/4 + y/4 + z/4;  //approssima la radice quadrata
            cout << i << " " << j <<"     "<< distanza_1 << " " << distanza_2 << " " << distanza_centri << endl;

            if (distanza_centri <= distanza_1 + distanza_2 + 1)  //bilancia la sottostima della norma al quadrato
            {
                array<unsigned int, 2> coppia = {i,j};
                coppie_vicine.push_back(coppia);
            }
        }
    }
}

}
