#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <chrono>
#include <omp.h>

using namespace std;


vector<int> Iskanje_v_mreži(vector<vector<int>> vektor2D, int ID) {
    
    vector<int> yx{ 0,0 };
    for (int ii = 0; ii < size(vektor2D) - 2; ii++) { // -2 da izključimo dodane robove
        for (int ji = 0; ji < size(vektor2D[ii]) - 2; ji++) {

            //cout << ji+1 << " " << ii+1 << endl;
            if (vektor2D[ii + 1][ji + 1] == ID) {
                yx[0] = ii + 1;
                yx[1] = ji + 1;
                break;
            }
        }
        if (vektor2D[yx[0]][yx[1]] == ID) break;
    }
    return yx;
}

int main() {

    //- Branje datoteke
    int Št_točk, Št_celic, Št_pogojev, Št_pogN;
    string bes;
    char c;
    double za1, za2;

    int ID_točk;
    vector<int> Seznam_ID_točk;
    double x, y;
    vector<double> Seznam_x, Seznam_y;
    vector<double> xy{ 0,0 };
    vector<vector<double>> Seznam_xy;
    vector<vector<int>> mreža;

    const double dx = 1.25;
    const double k = 24; // W/mK

    int ID_celic;
    vector<int> Seznam_ID_celic;
    int t1, t2, t3, t4;
    vector<int> Seznam_t1, Seznam_t2, Seznam_t3, Seznam_t4;

    int pogojN;
    vector<vector<int>> Seznam_pogoj;

    vector<char> rob_pog_kr;
    vector<double> rob_pog_temperatura;
    vector<double> rob_pog_toplotnitok;
    vector<double> rob_pog_prestop;


    //- BRANJE DATOTEKE
    ifstream file;
    file.open("primer2mreza.txt", ios::in);

    if (file.is_open()) {

        //- VOZLIŠČA
        file >> bes >> Št_točk;

        for (int i = 0; i < Št_točk; i++) {
            file >> ID_točk >> c >> x >> c >> y;

            x = round(x * 100) / 100;
            y = round(y * 100) / 100;

            Seznam_ID_točk.push_back(ID_točk);
            Seznam_x.push_back(x);
            Seznam_y.push_back(y);

            xy[0] = x;
            xy[1] = y;
            Seznam_xy.push_back(xy);
            //cout << Seznam_ID_točk[i] << "  " << Seznam_x[i] << "  " << Seznam_y[i] << "\n";
        }

        //- CELICE
        file >> bes >> Št_celic;

        for (int i = 0; i < Št_celic; i++) {
            file >> ID_celic >> c >> t1 >> c >> t2 >> c >> t3 >> c >> t4;
            Seznam_ID_celic.push_back(ID_celic);
            Seznam_t1.push_back(t1);
            Seznam_t2.push_back(t2);
            Seznam_t3.push_back(t3);
            Seznam_t4.push_back(t4);
            //cout << Seznam_ID_celic[i] << "  " << Seznam_t1[i] << "  " << Seznam_t2[i] << "  " << Seznam_t3[i] << "  " << Seznam_t4[i] << "\n";
        }

        //- ROBNI POGOJI
        file >> bes >> bes >> Št_pogojev;

        for (int i = 0; i < Št_pogojev; i++) {

            vector<int> Seznam_pogojN;

            file >> bes >> bes >> bes;
            //- POGOJ: temperatura
            if (bes == "temperatura") {

                file >> bes >> za1;
                //cout << bes << "  " << za1 << endl;
                rob_pog_temperatura.push_back(za1);
                rob_pog_toplotnitok.push_back(0);
                rob_pog_prestop.push_back(0);
                rob_pog_kr.push_back('t');

                file >> Št_pogN;

                for (int i = 0; i < Št_pogN; i++) {
                    file >> pogojN;
                    Seznam_pogojN.push_back(pogojN);
                    //cout << Seznam_pogojN[i] << endl;
                }
                Seznam_pogoj.push_back(Seznam_pogojN);
            }
            //- POGOJ: toplotni tok
            else if (bes == "toplotni") {

                file >> bes >> bes >> bes >> za1;
                //cout << bes << "  " << za1 << endl;
                rob_pog_temperatura.push_back(0);
                rob_pog_toplotnitok.push_back(za1);
                rob_pog_prestop.push_back(0);
                rob_pog_kr.push_back('k');

                file >> Št_pogN;

                for (int i = 0; i < Št_pogN; i++) {
                    file >> pogojN;
                    Seznam_pogojN.push_back(pogojN);
                    //cout << Seznam_pogojN[i] << endl;
                }
                Seznam_pogoj.push_back(Seznam_pogojN);
            }
            //- POGOJ: prestop
            else if (bes == "prestop") {

                file >> bes >> za1 >> bes >> za2;
                //cout << bes << "  " << za1 << "  " << za2 << endl;
                rob_pog_temperatura.push_back(za1);
                rob_pog_toplotnitok.push_back(0);
                rob_pog_prestop.push_back(za2);
                rob_pog_kr.push_back('p');

                file >> Št_pogN;

                for (int i = 0; i < Št_pogN; i++) {
                    file >> pogojN;
                    Seznam_pogojN.push_back(pogojN);
                    //cout << Seznam_pogojN[i] << endl;
                }
                Seznam_pogoj.push_back(Seznam_pogojN);
            }
            //- Neznan POGOJ
            else cout << "Neznana sintaksa pogoja:  " << bes << endl;
        }
        file.close();


        //- ZAPIS MREŽE
        int točke_x = round((*max_element(Seznam_x.begin(), Seznam_x.end()) - *min_element(Seznam_x.begin(), Seznam_x.end())) / dx) + 1;
        int točke_y = round((*max_element(Seznam_y.begin(), Seznam_y.end()) - *min_element(Seznam_y.begin(), Seznam_y.end())) / dx) + 1;
        mreža.assign(točke_y + 2, vector<int>(točke_x + 2, -1/*0*/));

        for (int i = 0; i < Št_točk; i++) { // y vrednosti so zrcaljene
            mreža[round((Seznam_y[i] - *min_element(Seznam_y.begin(), Seznam_y.end())) / dx) + 1][round((Seznam_x[i] - *min_element(Seznam_x.begin(), Seznam_x.end())) / dx) + 1] = Seznam_ID_točk[i] /*1*/;
        }// mreža[y][x]

        /*for (int i = 0; i < size(mreža); i++) {
            for (int j = 0; j < size(mreža[i]); j++) {
                cout << mreža[i][j] << " ";
            }
            cout << endl;
        }*/


        //- USTVARJANJE SISTEMA ENAČB
        vector<vector<double>> M(Št_točk, vector<double>(Št_točk, 0));
        vector<double> b(Št_točk, 0);

        vector<int> yx;
        int yi, xi;

        for (int i : Seznam_ID_točk) {

            yx = Iskanje_v_mreži(mreža, i);
            yi = yx[0];
            xi = yx[1];

            if (!(mreža[yi + 1][xi + 1] < 0 || mreža[yi + 1][xi - 1] < 0 || mreža[yi - 1][xi + 1] < 0 || mreža[yi - 1][xi - 1] < 0)) {

                M[i][i] = -4;
                M[i][mreža[yi + 1][xi]] = 1;
                M[i][mreža[yi - 1][xi]] = 1;
                M[i][mreža[yi][xi + 1]] = 1;
                M[i][mreža[yi][xi - 1]] = 1;
                //cout << 1 << "-" << mreža[yi+1][xi+1] << "-" << mreža[yi - 1][xi+1] << "-" << mreža[yi + 1][xi-1] << "-" << mreža[yi - 1][xi-1] << endl;
            }
            
            else {
                //cout << 2 << "-" << mreža[yi + 1][xi + 1] << "-" << mreža[yi - 1][xi + 1] << "-" << mreža[yi + 1][xi - 1] << "-" << mreža[yi - 1][xi - 1] << endl;
                for (int j = 0; j < Št_pogojev; j++) {

                    if (rob_pog_kr[j] == 't') {
                        if (find(Seznam_pogoj[j].begin(), Seznam_pogoj[j].end(), i) != Seznam_pogoj[j].end()) {
                            //cout << rob_pog_kr[j] << "  " << j + 1 << " " << i << endl;
                            
                            M[i][i] = 1;
                            b[i] = rob_pog_temperatura[j];

                            break;
                        }
                    }
                    else if (rob_pog_kr[j] == 'k') {
                        if (find(Seznam_pogoj[j].begin(), Seznam_pogoj[j].end(), i) != Seznam_pogoj[j].end()) {
                            //cout << rob_pog_kr[j] << "  " << j + 1 << " " << i << endl;

                            //- TOČKA NA ZUNANJEM KOTU
                            if (mreža[yi + 1][xi] < 0 && mreža[yi][xi + 1] < 0) {
                                M[i][i] = -2;
                                M[i][mreža[yi - 1][xi]] = 1;
                                M[i][mreža[yi][xi - 1]] = 1;
                            }
                            else if (mreža[yi + 1][xi] < 0 && mreža[yi][xi - 1] < 0) {
                                M[i][i] = -2;
                                M[i][mreža[yi - 1][xi]] = 1;
                                M[i][mreža[yi][xi + 1]] = 1;
                            }
                            else if (mreža[yi - 1][xi] < 0 && mreža[yi][xi + 1] < 0) {
                                M[i][i] = -2;
                                M[i][mreža[yi + 1][xi]] = 1;
                                M[i][mreža[yi][xi - 1]] = 1;
                            }
                            else if (mreža[yi - 1][xi] < 0 && mreža[yi][xi - 1] < 0) {
                                M[i][i] = -2;
                                M[i][mreža[yi + 1][xi]] = 1;
                                M[i][mreža[yi][xi + 1]] = 1;
                            }
                            //- TOČKA NA ROBU
                            else if (mreža[yi + 1][xi] < 0) {
                                M[i][i] = -4; // M[i][i] = M[i][mreža[yi][xi]] *-* i = mreža[yi][xi]
                                M[i][mreža[yi - 1][xi]] = 2;
                                M[i][mreža[yi][xi + 1]] = 1;
                                M[i][mreža[yi][xi - 1]] = 1;
                                //cout << " : " << mreža[yi][xi] << " " << mreža[yi + 1][xi] << " " << mreža[yi - 1][xi] << " " << mreža[yi][xi + 1] << " " << mreža[yi][xi - 1] << endl;
                            }
                            else if (mreža[yi - 1][xi] < 0) {
                                M[i][i] = -4;
                                M[i][mreža[yi + 1][xi]] = 2;
                                M[i][mreža[yi][xi + 1]] = 1;
                                M[i][mreža[yi][xi - 1]] = 1;
                                //cout << " : " << mreža[yi][xi] << " " << mreža[yi + 1][xi] << " " << mreža[yi - 1][xi] << " " << mreža[yi][xi + 1] << " " << mreža[yi][xi - 1] << endl;
                            }
                            else if (mreža[yi][xi + 1] < 0) {
                                M[i][i] = -4;
                                M[i][mreža[yi][xi - 1]] = 2;
                                M[i][mreža[yi + 1][xi]] = 1;
                                M[i][mreža[yi - 1][xi]] = 1;
                                //cout << " : " << mreža[yi][xi] << " " << mreža[yi + 1][xi] << " " << mreža[yi - 1][xi] << " " << mreža[yi][xi + 1] << " " << mreža[yi][xi - 1] << endl;
                            }
                            else if (mreža[yi][xi - 1] < 0) {
                                M[i][i] = -4;
                                M[i][mreža[yi][xi + 1]] = 2;
                                M[i][mreža[yi + 1][xi]] = 1;
                                M[i][mreža[yi - 1][xi]] = 1;
                                //cout << " : " << mreža[yi][xi] << " " << mreža[yi + 1][xi] << " " << mreža[yi - 1][xi] << " " << mreža[yi][xi + 1] << " " << mreža[yi][xi - 1] << endl;
                            }
                            //- TOČKA NA NOTRANJEM KOTU
                            else if (mreža[yi + 1][xi + 1] < 0) {
                                M[i][i] = -6;
                                M[i][mreža[yi - 1][xi]] = 2;
                                M[i][mreža[yi][xi - 1]] = 2;
                                M[i][mreža[yi + 1][xi]] = 1;
                                M[i][mreža[yi][xi + 1]] = 1;
                            }
                            else if (mreža[yi + 1][xi - 1] < 0) {
                                M[i][i] = -6;
                                M[i][mreža[yi - 1][xi]] = 2;
                                M[i][mreža[yi][xi + 1]] = 2;
                                M[i][mreža[yi + 1][xi]] = 1;
                                M[i][mreža[yi][xi - 1]] = 1;
                            }
                            else if (mreža[yi - 1][xi + 1] < 0) {
                                M[i][i] = -6;
                                M[i][mreža[yi + 1][xi]] = 2;
                                M[i][mreža[yi][xi - 1]] = 2;
                                M[i][mreža[yi - 1][xi]] = 1;
                                M[i][mreža[yi][xi + 1]] = 1;
                            }
                            else if (mreža[yi - 1][xi - 1] < 0) {
                                M[i][i] = -6;
                                M[i][mreža[yi + 1][xi]] = 2;
                                M[i][mreža[yi][xi + 1]] = 2;
                                M[i][mreža[yi - 1][xi]] = 1;
                                M[i][mreža[yi][xi - 1]] = 1;
                            }
                            else {
                                cout << "MANJKA" << endl;
                            }

                            b[i] = -2 * rob_pog_toplotnitok[j] * dx / k;

                            break;
                        }
                    }
                    else if (rob_pog_kr[j] == 'p') {
                        if (find(Seznam_pogoj[j].begin(), Seznam_pogoj[j].end(), i) != Seznam_pogoj[j].end()) {
                            //cout << rob_pog_kr[j] << "  " << j + 1 << endl;

                            //- TOČKA NA ZUNANJEM KOTU
                            if (mreža[yi + 1][xi] < 0 && mreža[yi][xi + 1] < 0) {
                                M[i][i] = -2 * (rob_pog_prestop[j] * dx / k + 1);
                                M[i][mreža[yi - 1][xi]] = 1;
                                M[i][mreža[yi][xi - 1]] = 1;
                            }
                            else if (mreža[yi + 1][xi] < 0 && mreža[yi][xi - 1] < 0) {
                                M[i][i] = -2 * (rob_pog_prestop[j] * dx / k + 1);
                                M[i][mreža[yi - 1][xi]] = 1;
                                M[i][mreža[yi][xi + 1]] = 1;
                            }
                            else if (mreža[yi - 1][xi] < 0 && mreža[yi][xi + 1] < 0) {
                                M[i][i] = -2 * (rob_pog_prestop[j] * dx / k + 1);
                                M[i][mreža[yi + 1][xi]] = 1;
                                M[i][mreža[yi][xi - 1]] = 1;
                            }
                            else if (mreža[yi - 1][xi] < 0 && mreža[yi][xi - 1] < 0) {
                                M[i][i] = -2 * (rob_pog_prestop[j] * dx / k + 1);
                                M[i][mreža[yi + 1][xi]] = 1;
                                M[i][mreža[yi][xi + 1]] = 1;
                            }
                            //-TOČKA NA ROBU 
                            else if (mreža[yi + 1][xi] < 0) {
                                M[i][i] = -2 * (rob_pog_prestop[j] * dx / k + 2);
                                M[i][mreža[yi - 1][xi]] = 2;
                                M[i][mreža[yi][xi + 1]] = 1;
                                M[i][mreža[yi][xi - 1]] = 1;
                                //cout << " : " << mreža[yi][xi] << " " << mreža[yi + 1][xi] << " " << mreža[yi - 1][xi] << " " << mreža[yi][xi + 1] << " " << mreža[yi][xi - 1] << endl;
                            }
                            else if (mreža[yi - 1][xi] < 0) {
                                M[i][i] = -2 * (rob_pog_prestop[j] * dx / k + 2);
                                M[i][mreža[yi + 1][xi]] = 2;
                                M[i][mreža[yi][xi + 1]] = 1;
                                M[i][mreža[yi][xi - 1]] = 1;
                                //cout << " : " << mreža[yi][xi] << " " << mreža[yi + 1][xi] << " " << mreža[yi - 1][xi] << " " << mreža[yi][xi + 1] << " " << mreža[yi][xi - 1] << endl;
                            }
                            else if (mreža[yi][xi + 1] < 0) {
                                M[i][i] = -2 * (rob_pog_prestop[j] * dx / k + 2);
                                M[i][mreža[yi][xi - 1]] = 2;
                                M[i][mreža[yi + 1][xi]] = 1;
                                M[i][mreža[yi - 1][xi]] = 1;
                                //cout << " : " << mreža[yi][xi] << " " << mreža[yi + 1][xi] << " " << mreža[yi - 1][xi] << " " << mreža[yi][xi + 1] << " " << mreža[yi][xi - 1] << endl;
                            }
                            else if (mreža[yi][xi - 1] < 0) {
                                M[i][i] = -2 * (rob_pog_prestop[j] * dx / k + 2);
                                M[i][mreža[yi][xi + 1]] = 2;
                                M[i][mreža[yi + 1][xi]] = 1;
                                M[i][mreža[yi - 1][xi]] = 1;
                                //cout << " : " << mreža[yi][xi] << " " << mreža[yi + 1][xi] << " " << mreža[yi - 1][xi] << " " << mreža[yi][xi + 1] << " " << mreža[yi][xi - 1] << endl;
                            }
                            //- TOČKA NA NOTRANJEM KOTU
                            else if (mreža[yi + 1][xi + 1] < 0) {
                                M[i][i] = -2 * (rob_pog_prestop[j] * dx / k + 3);
                                M[i][mreža[yi - 1][xi]] = 2;
                                M[i][mreža[yi][xi - 1]] = 2;
                                M[i][mreža[yi + 1][xi]] = 1;
                                M[i][mreža[yi][xi + 1]] = 1;
                            }
                            else if (mreža[yi + 1][xi - 1] < 0) {
                                M[i][i] = -2 * (rob_pog_prestop[j] * dx / k + 3);
                                M[i][mreža[yi - 1][xi]] = 2;
                                M[i][mreža[yi][xi + 1]] = 2;
                                M[i][mreža[yi + 1][xi]] = 1;
                                M[i][mreža[yi][xi - 1]] = 1;
                            }
                            else if (mreža[yi - 1][xi + 1] < 0) {
                                M[i][i] = -2 * (rob_pog_prestop[j] * dx / k + 3);
                                M[i][mreža[yi + 1][xi]] = 2;
                                M[i][mreža[yi][xi - 1]] = 2;
                                M[i][mreža[yi - 1][xi]] = 1;
                                M[i][mreža[yi][xi + 1]] = 1;
                            }
                            else if (mreža[yi - 1][xi - 1] < 0) {
                                M[i][i] = -2 * (rob_pog_prestop[j] * dx / k + 3);
                                M[i][mreža[yi + 1][xi]] = 2;
                                M[i][mreža[yi][xi + 1]] = 2;
                                M[i][mreža[yi - 1][xi]] = 1;
                                M[i][mreža[yi][xi - 1]] = 1;
                            }
                            else {
                                cout << "MANJKA" << endl;
                            }

                            b[i] = -2 * rob_pog_prestop[j] * dx * rob_pog_temperatura[j] / k;

                            break;
                        }
                    }
                    //- Zunanje else
                    else {
                        cout << "Neznan robni pogoj" << endl;
                    }
                }
            }
        }

        
        //- IZRAČUN SISTEMA ENAČB
        vector<double> T(Št_točk, 100);
        double d;
        int Št_iteracij = 600;

        /*auto čas0 = chrono::high_resolution_clock::now();

        for (int x = 0; x < Št_iteracij; x++) {

            for (int j = 0; j < Št_točk; j++) {
                d = b[j];

                for (int i = 0; i < Št_točk; i++) {

                    if (j != i) {
                        d = d - M[j][i] * T[i];
                    }
                }
                T[j] = d / M[j][j];
            }
            //cout << T[1670] << endl; // Zrcaljeno čez x os
        }
        
        auto čas1 = chrono::high_resolution_clock::now();
        chrono::duration<double> čas_razlika01 = čas1 - čas0;
        cout << "Cas singular//: " << čas_razlika01.count() << " s" << endl;*/


        T.assign(Št_točk, 100);
        auto čas2 = chrono::high_resolution_clock::now();

        for (int i = 0; i < Št_iteracij; i++) {

            for (int i = 0; i < Št_točk; i++) {

                double sum = 0;

                for (int j = 0; j < Št_točk; j++) {

                    if (i != j) sum += M[i][j] * T[j];
                }
                T[i] = (b[i] - sum) / M[i][i];
            }
            cout << T[1670] << endl;
        }

        auto čas3 = chrono::high_resolution_clock::now();
        chrono::duration<double> čas_razlika23 = čas3 - čas2;
        cout << "Cas parallel/*: " << čas_razlika23.count() << " s" << endl;


        /*for (int i = 0; i < 100; i++) {
            for (int j = 0; j < 100; j++) {
                cout << M[i][j] << " ";
            }
            cout << endl;
        }/*
        for (int i = 0; i < Št_točk; i++) {
            cout << i << " - " << b[i] << endl;
        }*/


        auto čas4 = chrono::high_resolution_clock::now();/////////
        //- IZPIS MATRIK
        ofstream newFile;
        newFile.open("M_b.txt", ios::out);

        if (newFile.is_open()) {
            newFile << size(M) << endl;

            for (int i = 0; i < size(M); i++) {
                for (int j = 0; j < size(M[i]); j++) {
                    newFile << M[i][j] << " ";
                }
                newFile << endl;
            }
            newFile << endl;
            newFile << size(b) << endl;

            for (int i = 0; i < size(b); i++) {
                newFile << b[i] << endl;
            }
        }
        else cout << "Beležke ni bilo možno napisati\n";

        //- IZPIS REŠITEV
        ofstream rešitev;
        rešitev.open("T.txt", ios::out);

        if (rešitev.is_open()) {
            rešitev << size(T) << endl;

            for (int i = 0; i < size(T); i++) {
                rešitev << T[i] << endl;
            }
            rešitev.close();
        }
        else cout << "Beležke ni bilo možno napisati\n";

        auto čas5 = chrono::high_resolution_clock::now();////////
        chrono::duration<double> čas_razlika45 = čas5 - čas4;/////////////
        cout << "Cas pisanja datoteke: " << čas_razlika45.count() << " s" << endl;///////////

        //- .VTK ZAPIS
        std::ofstream VTK("Prikaz.vtk");

        VTK << "# vtk DataFile Version 3.0" << endl;
        VTK << "Prikaz" << endl;
        VTK << "ASCII" << endl;
        VTK << "DATASET UNSTRUCTURED_GRID" << endl;

        VTK << "POINTS " << Št_točk << " float" << endl;
        for (int i = 0; i < Št_točk; i++) {
            VTK << Seznam_x[i] << " " << Seznam_y[i] << " 0" << endl;
        }
        VTK << endl;

        VTK << "CELLS " << Št_celic << " " << Št_celic * 5 << endl;
        for (int i = 0; i < Št_celic; i++) {
            VTK << 4 << " " << Seznam_t1[i] << " " << Seznam_t2[i] << " " << Seznam_t3[i] << " " << Seznam_t4[i] << endl;
        }
        VTK << endl;

        VTK << "CELL_TYPES " << Št_celic << endl;
        for (int i = 0; i < Št_celic; i++) {
            VTK << 9 << endl;
        }
        VTK << endl;

        VTK << "POINT_DATA " << Št_točk << endl;
        VTK << "SCALARS Temperatura  float 1" << endl;
        VTK << "LOOKUP_TABLE default" << endl;

        for (int i = 0; i < Št_točk; i++) {
            VTK << T[i] << endl;
        }

        VTK.close();
    }
    
    //- Datoteka ni odprta
    else cout << "Datoteka ni najdena";

    return 0;
}

