#include <stdio.h>
#include <complex.h>

void print_array(complex double *arr, int sz, const char *label) {
    printf("%s: ", label);
    for (int i = 0; i < sz; i++) {
        printf("(%lf + %lfi) ", creal(arr[i]), cimag(arr[i]));
    }
    printf("\n");
}

// Cyclically Shifts frequency domain array before discrete ifft
void ifft_shift(complex double *X, int sz) {
    int mid = sz / 2; 

    complex double temp[sz]; 

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

void fft_shift(complex double *X, int sz) 
{
    int mid = sz / 2; 

    complex double temp[sz]; 

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

int main() {
    complex double X_even[] = {1, 2, 3, 4, 5, 6, 7, 8};  // Even length array
    complex double X_odd[] = {1, 2, 3, 4, 5, 6, 7};      // Odd length array
    int sz_even = sizeof(X_even) / sizeof(X_even[0]);
    int sz_odd = sizeof(X_odd) / sizeof(X_odd[0]);

    printf("=== Testing Even-Length Array ===\n");
    print_array(X_even, sz_even, "Original X_even");

    fft_shift(X_even, sz_even);
    print_array(X_even, sz_even, "After fft_shift");

    ifft_shift(X_even, sz_even);
    print_array(X_even, sz_even, "After ifft_shift (Should match original)");

    printf("\n=== Testing Odd-Length Array ===\n");
    print_array(X_odd, sz_odd, "Original X_odd");

    fft_shift(X_odd, sz_odd);
    print_array(X_odd, sz_odd, "After fft_shift");

    ifft_shift(X_odd, sz_odd);
    print_array(X_odd, sz_odd, "After ifft_shift (Should match original)");

    return 0;
}