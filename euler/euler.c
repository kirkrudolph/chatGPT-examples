// euler.h

#ifndef EULER_H
#define EULER_H

#include <stdlib.h>

// Solves the given initial value problem using the forward Euler method
// with a fixed step size. The function f should compute the derivative
// of the state vector y at the given time t and state y, and store the
// result in the given derivative vector dydt. The initial state y0 and
// the time interval [t0, t1] should be given as input. The solution
// will be stored in the given array y. The number of steps n should
// be chosen such that t1 - t0 = n * h, where h is the fixed step size.
void forward_euler(void (*f)(double t, double* y, double* dydt), double t0, double t1, double* y0, double* y, size_t n);
void backward_euler(void (*f)(double t, double* y, double* dydt), double t0, double t1, double* y0, double* y, size_t n);

#endif // EULER_H

// euler.c

//#include "euler.h"

#include <math.h>

void forward_euler(void (*f)(double t, double* y, double* dydt), double t0, double t1, double* y0, double* y, size_t n)
{
    double h = (t1 - t0) / n;
    double t = t0;
    for (size_t i = 0; i < n; i++)
    {
        // Compute derivative at current time and state
        double dydt[2];
        f(t, y0, dydt);

        // Update state using forward Euler method
        for (size_t j = 0; j < 2; j++)
        {
            y0[j] += h * dydt[j];
            y[i * 2 + j] = y0[j];
        }

        t += h;
    }
}

void backward_euler(void (*f)(double t, double* y, double* dydt), double t0, double t1, double* y0, double* y, size_t n)
{
    double h = (t1 - t0) / n;
    double t = t0;
    for (size_t i = 0; i < n; i++)
    {
        // Compute derivative at current time and state
        double dydt[2];
        f(t, y0, dydt);

        // Update state using backward Euler method
        for (size_t j = 0; j < 2; j++)
        {
            y0[j] += h * dydt[j];
        }

        // Compute derivative at updated time and state
        f(t + h, y0, dydt);

        // Update state using backward Euler method
        for (size_t j = 0; j < 2; j++)
        {
            y0[j] += h * dydt[j];
            y[i * 2 + j] = y0[j];
        }

        t += h;
    }
}

void f(double t, double* y, double* dydt)
{
    dydt[0] = y[0];
}

#include <iostream>

int main(){

    { // Forward Euler
        double t0 = 0;
        double t1 = 1;
        double y0[2] = { 1, 0 };
        double y[2 * 10]; // Solution will be stored in this array

        forward_euler(f, t0, t1, y0, y, 10);

        std::cout << "Forward Euler: " << std::endl;
        for (int i = 0; i < 20; i +=2){
            std::cout << y[i]<< "\t" << y[i+1] << std::endl;
        }
    }

    { // Backward Euler 
        double t0 = 0;
        double t1 = 1;
        double y0[2] = { 1, 0 };
        double y[2 * 10]; // Solution will be stored in this array
        backward_euler(f, t0, t1, y0, y, 10);

        std::cout << "Backward Euler: " << std::endl;
        for (int i = 0; i < 20; i +=2){
            std::cout << y[i]<< "\t" << y[i+1] << std::endl;
        }
    }

    return 0;
}