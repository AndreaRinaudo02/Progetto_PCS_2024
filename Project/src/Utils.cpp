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
    if(!ImportFractures(file_path + "/FR362_data.txt", dfn, Plane))
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
        dfn.FractureId[i] = n;             //aggiunge Id all'elenco
        dfn.FracturesVertices[n] = {};     //aggiunge Id come chiave della mappa
        Plane.PlaneId[i]=n;

        ss >> n;      //numero di vertici della frattura

        dfn.FractureCoordinates[i].resize(n);     //riscala il vettore

        getline(file, riga);        //scarta intestazione

        for (int k=0; k < 3; k++)   //ciclo for (per ogni componente dei vertici)
        {
            getline(file, riga);    // riga con le coordinate
            stringstream ss(riga);

            for (unsigned int j =0; j < n; j++)                    // ciclo for (per ogni vertice)
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

    a = ((parametri[1][1] - parametri[0][1]) * (parametri[2][2] - parametri[0][2])) - ((parametri[2][1] - parametri[0][1]) * (parametri[1][2] - parametri[0][2]));
    b = -((parametri[1][0] - parametri[0][0]) * (parametri[2][2] - parametri[0][2])) - ((parametri[1][2] - parametri[0][2]) * (parametri[2][0] - parametri[0][0]));
    c = ((parametri[1][0] - parametri[0][0]) * (parametri[2][1] - parametri[0][1])) - ((parametri[1][1] - parametri[0][1]) * (parametri[2][0] - parametri[0][0]));

    double tol = 1e-12;

    if(abs(a)<tol)
    {
        a = 0;
    }

    if(abs(b)<tol)
    {
        b = 0;
    }
    if(abs(c)<tol)
    {
        c = 0;
    }

    d = - (a * parametri[0][0] + b * parametri[0][1] + c * parametri[0][2]);

    if(abs(d)<tol)
    {
        d = 0;
    }

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

    map<array<unsigned int, 2>, array<array<double, 3>, 2>> Retta={};   //memorizza i dati della retta di interezione in uma mappa coppia-Punto/direttrice
    RettaIntersezione(piano, coppie_vicine, Retta);

    //4 - 5 - 6
    IntersezioneLati(Retta, dfn);        //calcola le tracce

    const string file_name = "results.txt";       //stampa su file
    StampaTracce(file_name, dfn);
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

void StampaTracce(const string& file_name, DFN& dfn)
{
    ofstream file;
    file.open(file_name);     //apre il file

    file << "# Number of Traces" << endl;      //inizio stampa prima parte
    file <<"   "<< dfn.NumberTraces << endl;

    if(dfn.NumberTraces >0)
    {
        file << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;

        for (unsigned int i : dfn.TracesId)      //stampa prima parte
        {
            file << "  " << i << "; " << dfn.TracesFractures[i][0] << "; " << dfn.TracesFractures[i][1]
            << "; " << dfn.TracesVertices[i][0][0] << "; " << dfn.TracesVertices[i][0][1] << "; "
            << dfn.TracesVertices[i][0][2] << "; " << dfn.TracesVertices[i][1][0] << "; " << dfn.TracesVertices[i][1][1]
            << "; " << dfn.TracesVertices[i][1][2] << endl;
        }

        file << endl;

        file << "--------------------------------------------------------" << endl;

        file << endl;

        for (auto pair : dfn.FractureTraces)        //inizio stampa seconda parte
        {
            file << "# FractureId; NumTraces" << endl;
            file << "   " << pair.first << "; " << pair.second.size() << endl;
            file << "# TraceId; Tips; Length" << endl;

            vector<unsigned int> tracce = pair.second;               //occorre ordinare le tracce e misurarne la lunghezza
            map<unsigned int, double> lunghezza = {};                //mappa traccia-lunghezza

            for (unsigned int k = 0; k < tracce.size(); k = k+1)
            {
                array<array<double, 3>,2> vertici = dfn.TracesVertices[tracce[k]];

                double x = vertici[0][0] - vertici[1][0];
                double y = vertici[0][1] - vertici[1][1];
                double z = vertici[0][2] - vertici[1][2];

                double distanza = sqrt(x*x + y*y +z*z);      //calcola la lunghezza in norma 2

                lunghezza[tracce[k]] = distanza;


                unsigned int prossimo = tracce[k];           //algoritmo di sorting che ordina le tracce in base a Tips
                int j = k;

                array<unsigned int, 2> coppia1 = {prossimo,pair.first};
                array<unsigned int, 2> coppia2 = {tracce[j-1],pair.first};

                while ((j > 0) && dfn.Tips[coppia1] < dfn.Tips[coppia2])
                {
                    tracce[j] = tracce[j-1];
                    j = j-1;
                    coppia2 = {tracce[j-1], pair.first};
                }
                tracce[j] = prossimo;
            }

            for (unsigned int k = 0; k < tracce.size(); k = k+1)      //algoritmo di sorting che ordina le tracce, a parità di Tips, in base alla lunghezza
            {
                unsigned int prossimo = tracce[k];
                int j = k;

                array<unsigned int, 2> coppia1 = {prossimo,pair.first};
                array<unsigned int, 2> coppia2 = {tracce[j-1],pair.first};

                while ((j > 0) && lunghezza[tracce[j-1]] < lunghezza[prossimo] && dfn.Tips[coppia1] == dfn.Tips[coppia2])
                {
                    tracce[j] = tracce[j-1];
                    j = j-1;
                    coppia2 = {tracce[j-1], pair.first};
                }
                tracce[j] = prossimo;
            }

            for(unsigned int Id_traccia : tracce)        //stampa seconda parte
            {
                array<unsigned int, 2> coppia = {Id_traccia,pair.first};
                file << "   " << Id_traccia << "; " << dfn.Tips[coppia] << "; " << lunghezza[Id_traccia] << endl;
            }
            file << endl;
        }

        file.close();
    }
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
    array<double, 3> V={};       //memorizza le componenti della direttrice
    array<double, 3> P={};       //memorizza le componenti di un punto sulla retta

    for (auto coppia : coppie_vicine)
    {
        //la direttrice si trova intersecando i due piani contenenti le fratture

        V[0]= plane.Plane[coppia[0]][1]*plane.Plane[coppia[1]][2] - plane.Plane[coppia[1]][1]*plane.Plane[coppia[0]][2];
        V[1]= plane.Plane[coppia[1]][0]*plane.Plane[coppia[0]][2] - plane.Plane[coppia[0]][0]*plane.Plane[coppia[1]][2];
        V[2]= plane.Plane[coppia[0]][0]*plane.Plane[coppia[1]][1] - plane.Plane[coppia[0]][1]*plane.Plane[coppia[1]][0];


        double tol = 1e-12;

        if(abs(V[0]) < tol)
        {
            V[0] = 0;
        }

        if(abs(V[1]) < tol)
        {
            V[1] = 0;
        }

        if(abs(V[2]) < tol)
        {
            V[2] = 0;
        }

        if (V[0]!=0 || V[1]!=0 || V[2]!=0)  //esclude fratture parallele
        {
            //calcola le coordinate P intersecando i due piani precedenti con un terzo piano perpendicolare a entrambi

            Matrix3d C;
            C << plane.Plane[coppia[0]][0], plane.Plane[coppia[0]][1], plane.Plane[coppia[0]][2],
                plane.Plane[coppia[1]][0], plane.Plane[coppia[1]][1], plane.Plane[coppia[1]][2],
                V[0], V[1], V[2];

            Vector3d B;
            B << -plane.Plane[coppia[0]][3], -plane.Plane[coppia[1]][3], -0;

            Vector3d X = C.colPivHouseholderQr().solve(B);
            P[0]=X[0];
            P[1]=X[1];
            P[2]=X[2];

            if(abs(P[0]) < tol)
            {
                P[0] = 0;
            }

            if(abs(P[1]) < tol)
            {
                P[1] = 0;
            }

            if(abs(P[2]) < tol)
            {
                P[2] = 0;
            }

            array<array<double,3>,2> PV={P,V};
            Retta[coppia]=PV;
        }
    }
}

void IntersezioneLati(map<array<unsigned int, 2>,array<array<double, 3>, 2>>& Retta, DFN& dfn)
{
    unsigned int Id_traccia = 0;
    double prod_scal;

    for(const auto& pair : Retta)     //itera le coppie (e le rispettive rette di intersezione)
    {
        //double s,t

        array<vector<double>, 2> Parametri_t = {};     //P + tV identifica i punti di intrsezione tra la retta e i bordi delle fratture

        for(int i = 0; i<2; i++)      //itera sulle due fratture (della coppia)
        {
            unsigned int Id = pair.first[i];

            for (unsigned int j=0; j < dfn.FractureCoordinates[Id].size(); j++)   //itera sui vertici della frattura
            {
                array<array<double, 3>, 2> lato = {};      //memorizza gli estremi dei lati

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

                prod_scal = prod_vett1[0] * prod_vett2[0] + prod_vett1[1] * prod_vett2[1] + prod_vett1[2]* prod_vett2[2];

                /*
                if(pair.first[0] == 1 && pair.first[1] == 2)
                {
                    cout << pair.first[0] << " " << pair.first[1]<<endl;
                    cout << endl;
                    cout << "Frattura: " << Id << endl;
                    cout << endl;
                    cout << "P: "<<pair.second[0][0] << " " << pair.second[0][1] << " " << pair.second[0][2] <<endl;
                    cout << "V: "<<pair.second[1][0] << " " << pair.second[1][1] << " " << pair.second[1][2] <<endl;
                    cout <<endl;
                    cout << "Lato: " << j << endl;
                    cout << endl;
                    cout << "V1-P: " << retta1[0] << " " << retta1[1] << " " << retta1[2] <<endl;
                    cout << "V2-P: " << retta2[0] << " " << retta2[1] << " " << retta2[2] <<endl;
                    cout << "Prod. vett. 1: " << prod_vett1[0] << " "<< prod_vett1[1] << " "<< prod_vett1[2] << endl;
                    cout << "Prod. vett. 2: " << prod_vett2[0] << " "<< prod_vett2[1] << " "<< prod_vett2[2] << endl;
                    cout << "Prodotto scalare: " << prod_scal<< endl;
                    cout << endl;
                }
                */

                if (prod_scal < 0)    //indica se il lato viene intersecato dalla retta
                {
                    array<double, 3> d = {lato[0][0]-lato[1][0], lato[0][1]-lato[1][1], lato[0][2]-lato[1][2]};

                    MatrixXd C(3,2);
                    C << pair.second[1][0], -d[0], pair.second[1][1], -d[1], pair.second[1][2], -d[2];

                    Vector3d B;
                    B << lato[0][0] - pair.second[0][0], lato[0][1] - pair.second[0][1], lato[0][2] - pair.second[0][2] ;

                    VectorXd X = C.colPivHouseholderQr().solve(B);

                    Parametri_t[i].push_back(X[0]);
                }
            }
        }

        if(Parametri_t[0].size() !=0 && Parametri_t[1].size() != 0)    //interessano solo i casi di intersezione
        {
            if(max(Parametri_t[0][0], Parametri_t[0][1]) >= min(Parametri_t[1][0], Parametri_t[1][0]))      //assicura l'esistenza di una traccia
            {
                double T1_traccia, T2_traccia;
                T1_traccia = max(Parametri_t[0][0], Parametri_t[1][0]);     //coefficienti t degli estremi della traccia
                T2_traccia = min(Parametri_t[0][1], Parametri_t[1][1]);

                dfn.FractureTraces[pair.first[0]].push_back(Id_traccia);    //associa alla frattura la traccia

                dfn.FractureTraces[pair.first[1]].push_back(Id_traccia);

                dfn.TracesFractures[Id_traccia] = pair.first;               //associa alla traccia la coppia di fratture che la genera

                double tol = 1e-15;        //tolleranza

                if(abs(T1_traccia - Parametri_t[0][0]) < tol && abs(T2_traccia - Parametri_t[0][1]) < tol)     //controlla se la traccia è passante (con una tolleranza minima)
                {
                    array<unsigned int, 2> coppia = {Id_traccia,pair.first[0]};
                    dfn.Tips[coppia] = false;

                }

                else
                {
                    array<unsigned int, 2> coppia = {Id_traccia,pair.first[0]};
                    dfn.Tips[coppia] = true;

                }

                if(abs(T1_traccia - Parametri_t[1][0]) < tol && abs(T2_traccia - Parametri_t[1][1]) < tol)
                {
                    array<unsigned int, 2> coppia = {Id_traccia,pair.first[1]};
                    dfn.Tips[coppia] = false;
                }

                else
                {
                    array<unsigned int, 2> coppia = {Id_traccia,pair.first[1]};
                    dfn.Tips[coppia] = true;
                }

                double x1,x2,y1,y2,z1,z2;

                x1 = pair.second[0][0] + T1_traccia * pair.second[1][0];
                y1 = pair.second[0][1] + T1_traccia * pair.second[1][1];
                z1 = pair.second[0][2] + T1_traccia * pair.second[1][2];


                array<double, 3> estremo1 = {x1,y1,z1};      //primo estremo della traccia

                x2 = pair.second[0][0] + T2_traccia * pair.second[1][0];
                y2 = pair.second[0][1] + T2_traccia * pair.second[1][1];
                z2 = pair.second[0][2] + T2_traccia * pair.second[1][2];

                array<double, 3> estremo2 = {x2,y2,z2};     //secondo estremo della traccia

                array<array<double, 3>, 2> estremi = {estremo1, estremo2};


                dfn.TracesVertices[Id_traccia] = estremi;       //associa alla traccia i suoi estremi


                dfn.TracesId.push_back(Id_traccia);        //memorizza l'Id della traccia


                Id_traccia = Id_traccia+1;
            }
        }

    }

    dfn.NumberTraces = dfn.TracesId.size();        //calcola il numero di tracce totali
}

}
