#include <string>
#include <fstream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include "Utils.hpp"
#include "DFN.hpp"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
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
    vector<array<unsigned int, 2>> coppie_vicine = {}; //verranno memorizzate coppie di Id di fratture vicine per risparmiare tempo durante il calcolo delle tracce
    Fratture_vicine(dfn, coppie_vicine);

    //2 - 3

    map<array<unsigned int, 2>, array<array<double, 3>, 2>> Retta={};
    RettaIntersezione(piano, coppie_vicine, Retta);

    //4 - 5
    map<unsigned int, array<array<double, 3>, 2>> intersezioni = {};
    IntersezioneLati(Retta, dfn, intersezioni);


    /*
    cout << piano.Plane[coppie_vicine[0][1]][0]<< " "<< piano.Plane[coppie_vicine[0][1]][1]<< " "<< piano.Plane[coppie_vicine[0][1]][2]<< " "<< piano.Plane[coppie_vicine[0][1]][3]<< endl;
    cout << Retta[coppie_vicine[0]][0][0] << " " << Retta[coppie_vicine[0]][0][1] << " "<<  Retta[coppie_vicine[0]][0][2]<< endl;
    cout << Retta[coppie_vicine[0]][1][0] << " " << Retta[coppie_vicine[0]][1][1] << " "<<  Retta[coppie_vicine[0]][1][2]<< endl;

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

bool StampaTracce(const string& file_name, DFN& dfn)
{
    ofstream file;
    file.open(file_name);     //apre il file

    if(!file.is_open())
    {
        cerr << "Error: Unable to open file " << file_name << endl;
        return false;
    }

    file << "# Number of Traces" << endl;
    file << dfn.NumberTraces << endl;
    file << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
    for (unsigned int i : dfn.TracesId)
    {
        file << "  " << i << "; " << dfn.TracesFractures[i][0] << "; " << dfn.TracesFractures[i][1]
        << "; " << dfn.TracesCoordinates[i][0][0] << "; " << dfn.TracesCoordinates[i][0][1] << "; "
        << dfn.TracesCoordinates[i][0][2] << "; " << dfn.TracesCoordinates[i][1][0] << "; " << dfn.TracesCoordinates[i][1][1]
        << "; " << dfn.TracesCoordinates[i][1][2] << endl;

    }

    file.close();

    return true;
}

/*
bool TracceTips(const string& file_name, DFN& dfn)
{
    ofstream file;
    file.open(file_name);     //apre il file

    if(!file.is_open())
    {
        cerr << "Error: Unable to open file " << file_name << endl;
        return false;
    }

    file << "# FractureId; NumTraces" << endl;
    file << " " << FractureTips << " " <<dfn.NumberTipsTraces << endl;
    file << # TraceId; Tips; Length << endl;
    for (unsigned int i:NumberTipsTraces)
    {
        file << TraceTips[i] << " " << Tips[i] << " " << LenghtTips[i] << endl;
    }

    file.close();

    return true;
}
*/

void RettaIntersezione(Piano &plane,
                       vector<array<unsigned int, 2>>& coppie_vicine,
                       map<array<unsigned int, 2>,array<array<double, 3>, 2>>& Retta)
{
    array<double, 3> V={};
    array<double, 3> P={};

    for (auto coppia : coppie_vicine)
    {
        V[0]= plane.Plane[coppia[0]][1]*plane.Plane[coppia[1]][2] - plane.Plane[coppia[1]][1]*plane.Plane[coppia[0]][2];
        V[1]= plane.Plane[coppia[1]][0]*plane.Plane[coppia[0]][2] - plane.Plane[coppia[0]][0]*plane.Plane[coppia[1]][2];
        V[2]= plane.Plane[coppia[0]][0]*plane.Plane[coppia[1]][1] - plane.Plane[coppia[0]][1]*plane.Plane[coppia[1]][0];

        if (V[0]!=0 || V[1]!=0 || V[2]!=0)
        {
            Matrix3d C;
            C << plane.Plane[coppia[0]][0], plane.Plane[coppia[0]][1], plane.Plane[coppia[0]][2],
                plane.Plane[coppia[1]][0], plane.Plane[coppia[1]][1], plane.Plane[coppia[1]][02],
                V[0], V[1], V[2];

            Vector3d B;
            B << -plane.Plane[coppia[0]][3], -plane.Plane[coppia[1]][3], 0;

            Vector3d X = C.colPivHouseholderQr().solve(B);
            P[0]=X[0];
            P[1]=X[1];
            P[2]=X[2];


            array<array<double,3>,2> PV={P,V};
            Retta[coppia]=PV;
        }
    }
}

/*
void IntersezioneLati(map<array<unsigned int, 2>,array<array<double, 3>, 2>>& Retta, DFN& dfn,
                      map<unsigned int, array<array<double, 3>, 2>>& intersezioni)
{
    for(const auto& pair : Retta)
    {
        for(int i = 0; i<2; i++)
        {
            unsigned int Id = pair.first[i];
            //cout << Id <<endl;
            for (unsigned int j=0; j < dfn.FractureCoordinates[Id].size(); j++)
            {
                array<array<double, 3>, 2> lato = {};

                if (j == dfn.FractureCoordinates[Id].size()-1)
                {
                    lato = {dfn.FractureCoordinates[Id][j], dfn.FractureCoordinates[Id][0]};
                }
                else
                {
                    lato = {dfn.FractureCoordinates[Id][j], dfn.FractureCoordinates[Id][j+1]};
                }
                array<double, 3> retta1 = {lato[0][0]-pair.second[0][0], lato[0][1]-pair.second[0][1], lato[0][2]-pair.second[0][2]};
                array<double, 3> retta2 = {lato[1][0]-pair.second[0][0], lato[1][1]-pair.second[0][1], lato[1][2]-pair.second[0][2]};

                array<double, 3> prod_vett1 = {retta1[1]*pair.second[1][2]-retta1[2]*pair.second[1][1], retta1[2]*pair.second[1][0]-retta1[0]*pair.second[1][2], retta1[0]*pair.second[1][1]-retta1[1]*pair.second[1][0]};
                array<double, 3> prod_vett2 = {retta2[1]*pair.second[1][2]-retta2[2]*pair.second[1][1], retta2[2]*pair.second[1][0]-retta2[0]*pair.second[1][2], retta2[0]*pair.second[1][1]-retta2[1]*pair.second[1][0]};

                double prod_scal = prod_vett1[0] * prod_vett2[0] + prod_vett1[1] * prod_vett2[1] + prod_vett1[2]* prod_vett2[2];
                //cout << prod_scal << endl;
                if (prod_scal < 0)
                {
                    double t,s;
                    array<double, 3> d = {lato[0][0]-lato[1][0], lato[0][1]-lato[1][1], lato[0][2]-lato[1][2]};

                    MatrixXd C(3,2);
                    C << pair.second[1][0], -d[0], pair.second[1][1], -d[1], pair.second[1][2], -d[2];

                    Vector3d B;
                    B << lato[0][0] - pair.second[0][0], lato[0][1] - pair.second[0][1], lato[0][2] - pair.second[0][2] ;

                    VectorXd X = C.colPivHouseholderQr().solve(B);
                    t =X[0];
                    s =X[1];
                    //cout << t << " " << s << endl;

                }
            }
        }
    }
}
*/
}
