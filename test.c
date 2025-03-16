#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>

#define PI 3.14159265358979323846

// Cyclically Shifts frequency domain array before discrete ifft
void ifft_shift(double complex *X, int sz) 
{
    int mid = sz / 2; 
    double complex temp[sz]; 

    for (int i = 0; i < sz; i++) 
    {
        temp[i] = X[i];
    }
    int shifted_index = 0;
    for (int i = 0; i < sz; i++) 
    {
        shifted_index = (i + mid) % sz;
        X[i] = temp[shifted_index];
    }
}

// Cyclically Shifts frequency domain array before discrete fft
void fft_shift(double complex *X, int sz) 
{
    int mid = sz / 2; 
    double complex temp[sz]; 

    for (int i = 0; i < sz; i++) 
    {
        temp[i] = X[i];
    }
    int shifted_index = 0;
    for (int i = 0; i < sz; i++) 
    {
        shifted_index = (i + mid) % sz;
        X[shifted_index] = temp[i];
    }
}

// Mixed-radix FFT
void fft(double complex *X, double complex *Y, int sz) {
    if (sz <= 1) 
    {
        Y[0] = X[0];
        return;
    }

    // Factorize sz into smaller sizes
    int factors[sz];
    int n_factors = 0;
    int n = sz;
    for (int i = 2; i <= n; i++) 
    {
        while (n % i == 0) 
        {
            factors[n_factors++] = i;
            n /= i;
        }
    }

    // Recursive FFT
    double complex *temp = (double complex *)malloc(sz * sizeof(double complex));
    for (int i = 0; i < sz; i++) 
    {
        temp[i] = X[i];
    }

    int stride = 1;
    for (int i = 0; i < n_factors; i++) 
    {
        int factor = factors[i];
        int next_stride = stride * factor;

        for (int j = 0; j < sz; j += next_stride) 
        {
            for (int k = 0; k < stride; k++) 
            {
                double complex sum = 0;
                for (int l = 0; l < factor; l++) 
                {
                    double complex W = cexp(-I * 2.0 * PI * l * k / factor);
                    sum += temp[j + k + l * stride] * W;
                }
                Y[j + k] = sum;
            }
        }

        for (int j = 0; j < sz; j++) 
        {
            temp[j] = Y[j];
        }

        stride = next_stride;
    }

    free(temp);
    fft_shift(Y, sz);
}

// Mixed-radix IFFT
void ifft(double complex *X, double complex *Y, int sz) {
    ifft_shift(X, sz);

    if (sz <= 1) {
        Y[0] = X[0];
        return;
    }

    // Factorize sz into smaller sizes
    int factors[sz];
    int n_factors = 0;
    int n = sz;
    for (int i = 2; i <= n; i++) {
        while (n % i == 0) {
            factors[n_factors++] = i;
            n /= i;
        }
    }

    // Recursive IFFT
    double complex *temp = (double complex *)malloc(sz * sizeof(double complex));
    for (int i = 0; i < sz; i++) 
    {
        temp[i] = X[i];
    }

    int stride = 1;
    for (int i = 0; i < n_factors; i++) 
    {
        int factor = factors[i];
        int next_stride = stride * factor;

        for (int j = 0; j < sz; j += next_stride) 
        {
            for (int k = 0; k < stride; k++) 
            {
                double complex sum = 0;
                for (int l = 0; l < factor; l++) 
                {
                    double complex W = cexp(I * 2.0 * PI * l * k / factor);
                    sum += temp[j + k + l * stride] * W;
                }
                Y[j + k] = sum / factor;
            }
        }

        for (int j = 0; j < sz; j++) 
        {
            temp[j] = Y[j];
        }

        stride = next_stride;
    }

    free(temp);
}

int main() {
    int sz = 8;
    double complex X[] = {1, 1, 1, 1, 0, 0, 0, 0};
    double complex Y[sz];

    fft(X, Y, sz);

    printf("FFT:\n");
    for (int i = 0; i < sz; i++) {
        printf("%f + %fi\n", creal(Y[i]), cimag(Y[i]));
    }

    ifft(Y, X, sz);

    printf("\nIFFT:\n");
    for (int i = 0; i < sz; i++) {
        printf("%f + %fi\n", creal(X[i]), cimag(X[i]));
    }

    return 0;
}