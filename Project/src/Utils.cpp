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
    if(!ImportFractures(file_path + "/FR82_data.txt", dfn, Plane))
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
    vector<array<unsigned int, 2>> coppie_vicine = {};      //verranno memorizzate coppie di Id di fratture vicine per risparmiare tempo durante il calcolo delle tracce
    Fratture_vicine(dfn, coppie_vicine);

    /*
    for (unsigned int i = 0; i < coppie_vicine.size(); i++)
    {
        cout << "coppie: " << coppie_vicine[i][0]<< " " <<coppie_vicine[i][1] << endl;
    }
    */

}

void Fratture_vicine(DFN& dfn, vector<array<unsigned int, 2>>& coppie_vicine)
{
    vector<array<double , 4>> bolle = {};     //contiene i dati riguardanti le bolle intorno alle fratture
    bolle.resize(dfn.NumberFractures);
    Crea_bolle(dfn, bolle);
    double x,y,z;

    for(unsigned int a = 0; a < dfn.NumberFractures - 1; a++)
    {
        for(unsigned int b = a + 1; b < dfn.NumberFractures; b++)
        {
            x = bolle[a][0] - bolle[b][0];
            y = bolle[a][1] - bolle[b][1];
            z = bolle[a][2] - bolle[b][2];

            double distanza = abs(x) + abs(y) + abs(z);     //calcola la distanza tra i centri di due bolle (in norma 1)

            if(distanza < bolle[a][3] + bolle[b][3])        //una coppia viene considerata vicina se i due centri distano meno rispetto alla somma dei raggi delle due bolle in esse contenuti
            {
                array<unsigned int, 2> coppia = {a,b};
                coppie_vicine.push_back(coppia);
            }
        }
    }
}

void Crea_bolle(DFN& dfn, vector<array<double, 4>>& bolle)
{
    double x;
    double y;
    double z;

    for(unsigned int i = 0; i < dfn.NumberFractures; i++)
    {
        array<double, 4> bolla = {};        //contiene i dati della bolla (coordinate del cenntro (x, y, z) e il raggio)

        vector<array<double, 3>> coordinate = dfn.FracturesVertices[i];

        x = 0;
        y = 0;
        z = 0;

        for(unsigned int k =0; k < coordinate.size(); k++)
        {
            x += coordinate[k][0];
            y += coordinate[k][1];
            z += coordinate[k][2];
        }

        bolla[0] = x/coordinate.size();    //calcola il centro della bolla facendo la media delle coordinate dei vertici
        bolla[1] = y/coordinate.size();
        bolla[2] = z/coordinate.size();

        vector<double> distanze = {};

        for(unsigned int j = 0; j < coordinate.size(); j++)
        {
            x = bolla[0] - coordinate[j][0];
            y = bolla[1] - coordinate[j][1];
            z = bolla[2] - coordinate[j][2];

            double distanza = abs(x) + abs(y) + abs(z);   //calcola la distanza dal centro di ogni verice (calcolandola in norma 1)
            distanze.push_back(distanza);
        }

        double raggio = *max_element(distanze.begin(), distanze.end());   //il raggio della bolla è pari alla distanza tra il vertice più lontano dal centro e il centro

        bolla[3] = raggio;

        bolle[i] = bolla;
    }
}


}
