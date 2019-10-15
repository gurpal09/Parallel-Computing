//============================================================================
// Name        : coutte.cpp
// Author      : Gurpal Singh
// Date	       : 02/04/2017
// Description : Exact and Numerical solution of 1D Coutte flow
//
// To Compile  : g++ -Wall coutte.cpp -std=c++11 -o -coutte.exe -O3
// To Execute  : ./coutte.exe <Reynolds number> <Pressure Gradient>
//============================================================================

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <string>
using namespace std;

#define NY 21 // The number of points in y-direction with the boundaries included

typedef double REAL;
typedef int INT;

int main(INT argc, char *argv[])
{

    if (argc != 3) {
        perror("Command-line usage: executableName <Reynolds number> <Pressure gradient>");
        exit(1);
    }

    const REAL Re = atof(argv[1]);                   // Reynolds Number
    const REAL dPdx = atof(argv[2]);                 // Pressure Gradient
    const REAL h = 1.0f;                             // Plate Height
    const REAL uPlate = 1.0f;                        // Plate Velocity of Top Plate
    const REAL rho = 1.0f;                           // Density of fluid
    const REAL nu = (uPlate * h) / Re;               // Kinematic Viscosity
    const REAL mu = nu * rho;                        // Dynamic Viscosity
    const REAL dy = h / (NY - 1);                    // Space between points in y direction
    const REAL dt = (0.5f * dy * dy) / nu;           // time step
    const REAL timeEnd = (h * h) / nu;               // Time to reach Steady-State
    const REAL nTimeSteps = (INT)ceil(timeEnd / dt); // Number of Time steps required

    REAL *unew = new double[NY](); // Array to hold velocity for n+1 time steps
    REAL *u = new double[NY]();    // Array to hold velocity for n time step
    REAL *y = new double[NY]();    // Array to keep track of y coordinate
    REAL *tmp;                     // Pointer to a Real, needed to swap pointers

    // Applying Boundary Conditions
    u[0] = 0.f;
    unew[0] = 0.f;
    u[NY - 1] = uPlate;
    unew[NY - 1] = uPlate;

    // Time March using the finite difference formula
    for (int i = 1; i <= nTimeSteps; i++) {
        for (INT j = 1; j < NY - 1; j++) {
            unew[j] = dt * (-dPdx / rho + mu * (u[j + 1] - 2.0f * u[j] + u[j - 1]) / (dy * dy))
                      + u[j];
        }
        // Swapping the pointers
        tmp = u;
        u = unew;
        unew = tmp;
    }

    // Exact Solution and saving it in unew
    for (INT i = 0; i < NY; i++) {
        y[i] = i * dy;
        unew[i] = uPlate * y[i] / h + dPdx / (2.f * mu) * (y[i] * y[i] - h * y[i]);
    }

    // Save the output to a file an manipulate into 2 columns
    ofstream myfile;
    myfile << "\n";
    myfile.open("CoutteData.txt");
    myfile << "Data for Reynolds Number:" << Re << left << setw(10)
           << " and Pressure Gradient:" << dPdx << endl;
    myfile << std::fixed;
    myfile << std::setprecision(5);
    myfile << left << setw(20) << "Exact solution" << setw(20) << "Numerical Solution" << endl;
    myfile << "\n";
    for (int i = 0; i < 21; i++) {
        myfile << left << setw(20) << unew[i] << setw(20) << u[i] << endl;
    }
    myfile.close();
    cout << "The data was written to the file successfully" << endl;

    // Freeing the allocated memory
    delete[] unew;
    delete[] u;
    delete[] y;

    return 0;
}
