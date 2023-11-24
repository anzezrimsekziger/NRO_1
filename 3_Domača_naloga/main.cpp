/******************************************************************************

Welcome to GDB Online.
  GDB online is an online compiler and debugger tool for C, C++, Python, PHP, Ruby, 
  C#, OCaml, VB, Perl, Swift, Prolog, Javascript, Pascal, COBOL, HTML, CSS, JS
  Code, Compile, Run and Debug online from anywhere in world.

*******************************************************************************/
#include <iostream>
#include <cmath>

// Funkcija za izračun arctan(x) s Taylorjevo vrsto
double arctan(double x) {
    double sum = 0.0;
    double term = x;
    double power = x * x;
    int n = 1;

    while (fabs(term) > 1e-8) {
        sum += term;
        term = -term * power / (2 * n + 1);
        power *= x * x;
        n++;
    }

    return sum;
}

// Funkcija za izračun e^(3x)
double exp_3x(double x) {
    return exp(3 * x);
}

// Funkcija za integracijo s trapezno metodo
double trapezoidal_integration(double (*f)(double), double a, double b, int n) {
    double h = (b - a) / n;
    double sum = 0.5 * (f(a) + f(b));

    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        sum += f(x);
    }

    return h * sum;
}

int main() {
    double pi = acos(-1.0); // vrednost PI
    double a = 0.0; // spodnja meja integrala
    double b = pi / 4.0; // zgornja meja integrala
    int n = 1000; // število podintervalov

    double result = trapezoidal_integration([](double x) {
        return exp_3x(x) * arctan(x / 2.0);
    }, a, b, n);

    std::cout << "Ocena integrala: " << result << std::endl;

    return 0;
}

