#include <iostream>
#include <fstream>
#include <vector>

// Funkcija za izračun odvoda prvega reda z diferenčno shemo za robno vrednost na levi strani
double forwardDifference(std::vector<std::pair<double, double>>& data, int idx, double dx) {
    return (-3 * data[idx].second + 4 * data[idx + 1].second - data[idx + 2].second) / (2 * dx);
}

// Funkcija za izračun odvoda prvega reda z diferenčno shemo za robno vrednost na desni strani
double backwardDifference(std::vector<std::pair<double, double>>& data, int idx, double dx) {
    return (3 * data[idx].second - 4 * data[idx - 1].second + data[idx - 2].second) / (2 * dx);
}

// Funkcija za izračun odvoda prvega reda z centralno diferenčno shemo
double centralDifference(std::vector<std::pair<double, double>>& data, int idx, double dx) {
    return (data[idx + 1].second - data[idx - 1].second) / (2 * dx);
}

int main() {
    std::ifstream inputFile("poly.txt");
    if (!inputFile.is_open()) {
        std::cerr << "Unable to open file poly.txt\n";
        return 1;
    }

    int numOfValues;
    inputFile >> numOfValues;

    std::vector<std::pair<double, double>> data;
    double x, fx;
    for (int i = 0; i < numOfValues; ++i) {
        inputFile >> x >> fx;
        data.emplace_back(x, fx);
    }

    inputFile.close();

    double dx = 0.0294117647058822; // ∆x
    std::ofstream outputFile("diff_poly.txt");
    if (!outputFile.is_open()) {
        std::cerr << "Unable to create file diff_poly.txt\n";
        return 1;
    }

    // Izračun odvodov in zapis rezultatov v datoteko
    for (int i = 0; i < numOfValues; ++i) {
        double derivative;
        if (i == 0) {
            derivative = forwardDifference(data, i, dx);
        } else if (i == numOfValues - 1) {
            derivative = backwardDifference(data, i, dx);
        } else {
            derivative = centralDifference(data, i, dx);
        }
        outputFile << "Odvod f'(" << data[i].first << ") = " << derivative << std::endl;
    }

    outputFile.close();

    std::cout << "Differential results saved in diff_poly.txt\n";

    return 0;
}
