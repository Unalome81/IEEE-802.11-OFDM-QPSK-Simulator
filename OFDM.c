#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <complex.h>  // Required for complex numbers

#define N_FFT 64

#define PI 3.14159265358979323846

int len_Tx_Signal = 0;

unsigned char message[] = "The Supreme Lord Shree Krishna said: I taught this eternal science of Yog to the Sun God, Vivasvan, who passed it on to Manu; and Manu, in turn, instructed it to Ikshvaku.";

//"The Supreme Lord Shree Krishna said: I taught this eternal science of Yog to the Sun God, Vivasvan, who passed it on to Manu; and Manu, in turn, instructed it to Ikshvaku.";

double RRC_Filter_Tx[21] = {-0.000454720514876223, 0.00353689555574986, -0.00714560809091226, 0.00757906190517828, 0.00214368242727367, -0.0106106866672496, 0.0300115539818315, -0.0530534333362480, -0.0750288849545787, 0.409168714634052, 0.803738600397980, 0.409168714634052, -0.0750288849545787, -0.0530534333362480, 0.0300115539818315, -0.0106106866672496, 0.00214368242727367, 0.00757906190517828, -0.00714560809091226, 0.00353689555574986, -0.000454720514876223};


// Cyclically Shifts frequency domain array before discrete ifft
void ifft_shift(complex double *X) 
{
    int mid = N_FFT / 2;

    for(int i = 0; i < mid; ++i)
    {
        int temp = X[i];
        X[i] = X[i + mid];
        X[i + mid] = temp;
    }
}


// Discete Inverse Fourier Transform : O(N^2)
void ifft(complex double *X, complex double *Y) // X = Frequency Domain, Y = Time Domain
{
    ifft_shift(X);

    complex double W;
    for(int i = 0; i < N_FFT; ++i)
    {
        Y[i] = 0;
        for(int j = 0; j < N_FFT; ++j)
        {
            W = cexp((I * 2.0 * PI * i * j) / N_FFT);
            Y[i] += X[j] * W;
        }
        Y[i] /= N_FFT;
    }
}


void Slice_Repeater(complex double *X, complex double *Y, int starti, int lo, int hi, int r) // ind is starting index of Y i.e Output, lo and hi are the starting and ending indices of X(input), r is the number of times we are repeating
{
    int index = starti; 

    for(int i = 0; i < r; ++i)
    {
        for(int j = lo; j < hi; ++j)
        {
            Y[index] = X[j];
            index++;
        }        
    }
}

void Display(complex double *X, int sz)
{
    for (int i = 0; i < sz; i++) 
    {
        printf("%d  " "%lf + %lf*I \n",i+1, creal(X[i]), cimag(X[i]));
    }
    printf("\n");
}

void Preamble_Generator(double scale, double complex *P_k, double complex *virtual_subcarrier, double complex *Preamble, int type)
{
    //complex double Preamble[160];
    // Scaling S_k
    for (int i = 0; i < 53; i++) 
    {
        P_k[i] *= scale;
    }    

    double complex preamble_freq[N_FFT] = {0};
    double complex preamble_time[N_FFT] = {0};

    // Fill the first 6 elements with the first 6 virtual subcarriers (which are zeros)
    for (int i = 0; i < 6; i++) {
        preamble_freq[i] = virtual_subcarrier[i];
    }
    
    // Fill the next 53 elements with P_k values
    for (int i = 0; i < 53; i++) {
        preamble_freq[6 + i] = P_k[i];
    }
    
    // Fill the last 5 elements with the remaining virtual subcarriers (indices 6 to 10)
    for (int i = 0; i < 5; i++) {
        preamble_freq[6 + 53 + i] = virtual_subcarrier[6 + i];
    }

    // Print the short preamble in frequency domain
    // printf("Preamble (Frequency Domain):\n");
    //Display(preamble_freq, N_FFT);

    ifft(preamble_freq, preamble_time);

    // printf("Preamble Shift (Time Domain):\n");
    // Display(preamble_time, N_FFT);
    
    // if(type == 0)
    // {
    //     printf("Short Preamble:\n");
    // }
    // else 
    // {
    //     printf("Long Preamble:\n");
    // }

    if(type == 0) // Short Preamble
        Slice_Repeater(preamble_time, Preamble, 0, 0, 16, 10);
    else
    {
        Slice_Repeater(preamble_time, Preamble, 0, 32, N_FFT, 1);
        Slice_Repeater(preamble_time, Preamble, 32, 0, N_FFT, 2);
    }

    //Display(Preamble, 160); 
}

complex double* Allocate_Array_1D(int size) {
    complex double *array = (complex double *)calloc(size, sizeof(complex double));
    if (array == NULL) {
        printf("Memory allocation failed for 1D array!\n");
        exit(1);
    }
    return array;
}

// Function to allocate a 2D array of complex double
complex double** Allocate_Array_2D(int rows, int cols) {
    complex double **array = (complex double **)malloc(rows * sizeof(complex double *));
    if (array == NULL) {
        printf("Memory allocation failed for 2D array (row pointers)!\n");
        exit(1);
    }

    for (int i = 0; i < rows; i++) {
        array[i] = (complex double *)calloc(cols, sizeof(complex double));
        if (array[i] == NULL) {
            printf("Memory allocation failed for 2D array (row %d)!\n", i);
            exit(1);
        }
    }
    return array;
}

void Decimal_To_Binary(int msg, int* msg_bits) 
{
    for(int i = 0; i < 8; ++i)
    {
        msg_bits[i] = 0;
    }    

    for(int i = 7; msg > 0; --i)
    {
        msg_bits[i] = msg % 2;
        msg /= 2;
    }
}

void QPSK_Modulator(complex double **Data_Payload, complex double **Data_Payload_Mod, int rows)
{
    for(int i = 0; i < rows; ++i)
    {
        for(int j = 0; j < 96; j += 2)
        {
            int a = Data_Payload[i][j], b = Data_Payload[i][j + 1];

            if(a == 0 && b == 0)
                Data_Payload_Mod[i][j / 2] = (1 + 1*I) / sqrt(2);
            else if(a == 0 && b == 1)
                Data_Payload_Mod[i][j / 2] = (-1 + 1*I) / sqrt(2);
            else if(a == 1 && b == 0)
                Data_Payload_Mod[i][j / 2] = (-1 - 1*I) / sqrt(2);
            else if(a == 1 && b == 1)
                Data_Payload_Mod[i][j / 2] = (1 - 1*I) / sqrt(2);
        }
    }
}

complex double *Transmitter()
{
    double complex Short_Preamble[160], Long_Preamble[160];

    // Scale Factor
    double scale = sqrt(13.0 / 6.0);

    double complex virtual_subcarrier[11] = {0};

    double complex S_k[53] = 
    {
    0, 0, 1 + 1*I, 0, 0, 0, -1 - 1*I, 0, 0, 0, 
    1 + 1*I, 0, 0, 0, -1 - 1*I, 0, 0, 0, -1 - 1*I, 0, 0, 0, 
    1 + 1*I, 0, 0, 0, 0, 0, 0, 0, -1 - 1*I, 0, 0, 0, 
    -1 - 1*I, 0, 0, 0, 1 + 1*I, 0, 0, 0, 1 + 1*I, 0, 0, 0, 
    1 + 1*I, 0, 0, 0, 1 + 1*I, 0, 0
    };

    Preamble_Generator(scale, S_k, virtual_subcarrier, Short_Preamble, 0);

    double complex L_k[53] ={1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1};


    scale = 1;
    Preamble_Generator(scale, L_k, virtual_subcarrier, Long_Preamble, 1);    

    // Preparing Payload

    int M = 4; // QPSK modulation
    int bits_per_symb = log2(M); // Bits per symbol
    int n_bits = 96; // Number of bits in one frame

    int message_size = sizeof(message)/sizeof(message[0]) - 1; // removing null character from message size
    int rows = ceil(8 * message_size / 96.0);

    int total_bits = rows * 96;

    complex double *Data = Allocate_Array_1D(total_bits);

    if (Data == NULL) 
    {
        printf("Memory allocation failed!\n");
        return NULL;
    }

    // Encoder 
    int msg_bits[8] = {0};

    int index = 0;

    for(int i = 0; i < message_size; ++i)
    {
        Decimal_To_Binary(message[i], msg_bits);

        for(int j = 0; j < 8; ++j)
        {
            Data[index] = msg_bits[j];
            index++;
        }
    }

    // Creating Data Payloads

    complex double **Data_Payload = Allocate_Array_2D(rows, 96);

    for(int i = 0; i < rows; ++i)
    {
        int row_start = i * 96;
        for(int j = 0; j < 96; ++j)
        {
            Data_Payload[i][j] = Data[row_start + j];
        }
    }

    // Combining consecutive bits for QPSK modulation: Creating 1 X 48 complex numbers from 1 * 2 real numbers

    complex double **Data_Payload_Mod = Allocate_Array_2D(rows, 48);

    QPSK_Modulator(Data_Payload, Data_Payload_Mod, rows);

    // Data Frames

    int pilot[] = {1,1,1,-1};
    complex double **Data_Frames = Allocate_Array_2D(rows, 64);

    for(int i = 0; i < rows; ++i)
    {
        Slice_Repeater(virtual_subcarrier, Data_Frames[i], 0, 0, 6, 1); // Filled = 6

        Slice_Repeater(Data_Payload_Mod[i], Data_Frames[i], 6, 0, 5, 1);    // Filled = 11
        Data_Payload_Mod[i][11] = pilot[0];                                 // Filled = 12

        Slice_Repeater(Data_Payload_Mod[i], Data_Frames[i], 12, 5, 18, 1);  // Filled = 25
        Data_Payload_Mod[i][25] = pilot[1];                                 // Filled = 26

        Slice_Repeater(Data_Payload_Mod[i], Data_Frames[i], 26, 18, 24, 1); // Filled = 32
        
        Data_Payload_Mod[i][32] = 0;                                        // Filled = 33

        Slice_Repeater(Data_Payload_Mod[i], Data_Frames[i], 33, 24, 30, 1); // Filled  = 39
        Data_Payload_Mod[i][39] = pilot[2];                                 // FIlled = 40

        Slice_Repeater(Data_Payload_Mod[i], Data_Frames[i], 40, 30, 43, 1); // Filled  = 53
        Data_Payload_Mod[i][53] = pilot[3];                                 // Filled  = 54

        Slice_Repeater(Data_Payload_Mod[i], Data_Frames[i], 54, 43, 48, 1); // Filled  = 59
        Slice_Repeater(virtual_subcarrier, Data_Frames[i], 59, 6, 11, 1); // Filled = 64     
    }

    complex double **Data_Frames_IFFT = Allocate_Array_2D(rows, 64);

    for(int i = 0; i < rows; ++i)
    {
        ifft(Data_Frames[i], Data_Frames_IFFT[i]);
    }

    complex double **Data_Frames_Time_TX = Allocate_Array_2D(rows, 80);

    for(int i = 0; i < rows; ++i)
    {
        Slice_Repeater(Data_Frames_IFFT[i], Data_Frames_Time_TX[i], 0, 48, 64, 1);
        Slice_Repeater(Data_Frames_IFFT[i], Data_Frames_Time_TX[i], 16, 0, 64, 1);
    }

    int Data_Frame_Size = rows * 80 + 160 + 160;
    complex double *Data_Frame_TX = Allocate_Array_1D(Data_Frame_Size);

    Slice_Repeater(Short_Preamble, Data_Frame_TX, 0, 0, 160, 1);   // Filled = 160
    Slice_Repeater(Long_Preamble, Data_Frame_TX, 160, 0, 160, 1);  // Filled = 320

    int start = 320;

    for(int i = 0; i < rows; ++i)
    {
        Slice_Repeater(Data_Frames_Time_TX[i], Data_Frame_TX, start, 0, 80, 1);
        start += 80;
    }

    complex double *Data_Frame_TX_Oversamp = Allocate_Array_1D(2*Data_Frame_Size);

    for(int i = 0; i < Data_Frame_Size; ++i)
    {
        Data_Frame_TX_Oversamp[2*i] = Data_Frame_TX[i];
        Data_Frame_TX_Oversamp[2*i + 1] = 0;
    }  

    // Root Raised Cosine filter (Rx)

    double rolloff_factor_rx = 0.5;
    int length_rrc_rx = 10; 
    int oversampling_rate_rx = 2; 
    
    // Convoulation of the Tx Frame with RRC Filter

    int len_inp_sig = 2*Data_Frame_Size;
    int len_coeff = 21; 
    int len_out_sig = len_inp_sig + len_coeff - 1; 

    complex double *inp_signal_padded = Allocate_Array_1D(len_inp_sig + 2 * (len_coeff - 1));
    Slice_Repeater(Data_Frame_TX_Oversamp, inp_signal_padded, 20, 0, len_inp_sig, 1);

    complex double *Tx_signal = Allocate_Array_1D(len_out_sig);

    for(int i = 0; i < len_out_sig; ++i)
    {
        for(int j = 0; j < len_coeff; ++j)
        {
            Tx_signal[i] += RRC_Filter_Tx[j] * inp_signal_padded[i + j - 1];
        }
    }

    int len_out_signal_repeated = 10 * len_out_sig;
    complex double *Tx_signal_repeated = Allocate_Array_1D(len_out_signal_repeated);

    Slice_Repeater(Tx_signal, Tx_signal_repeated, 0, 0, len_out_sig, 10);

    len_Tx_Signal = len_out_signal_repeated;

    return Tx_signal;
}

double gaussian_noise(double mean, double variance) {
    double u1, u2, z0;
    
    u1 = ((double) rand() + 1.0) / ((double) RAND_MAX + 1.0);
    u2 = ((double) rand() + 1.0) / ((double) RAND_MAX + 1.0);

    z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);

    return mean + sqrt(variance) * z0;
}

complex double* Transmission_Over_Air(complex double* TX_signal, double snr)
{
    complex double* TX_OTA_signal = Allocate_Array_1D(len_Tx_Signal); 

    double Tx_signal_power = 0.0;
    for (int i = 0; i < len_Tx_Signal; i++) 
    {
        Tx_signal_power += sqrt(cabs(TX_signal[i]));
    }
    Tx_signal_power /= len_Tx_Signal;  // Mean power

    double snr_linear = pow(10, Tx_signal_power / 10);

    double noise_power = Tx_signal_power / snr_linear;

    for(int i = 0 ; i < len_Tx_Signal; ++i)
    {
        TX_OTA_signal[i] += sqrt(noise_power) * (gaussian_noise(0, 1) + I*gaussian_noise(0, 1));
    }

    return TX_OTA_signal;
}

void Receiver()
{
    
}


int main() 
{
    complex double* TX_signal = Transmitter();

    printf("%d", len_Tx_Signal);

    complex double* Tx_OTA_signal = Transmission_Over_Air(TX_signal, 100);

    return 0;    
}