#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <complex.h>  // Required for complex numbers
#include <time.h>

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Global Arrays and variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

#define N_FFT 64

#define PI 3.14159265358979323846

#define fc_hz 5e9  // Carrier frequency (5 GHz)
#define fs_hz 20e6  // Sampling frequency (20 MHz)
#define ts_sec (1/fs_hz)  // Time period (50 ns)
#define num_snr 35

unsigned char message[] = "The Supreme Lord Shree Krishna said: I taught this eternal science of Yog to the Sun God, Vivasvan, who passed it on to Manu; and Manu, in turn, instructed it to Ikshvaku.";
//"The Supreme Lord Shree Krishna said: I taught this eternal science of Yog to the Sun God, Vivasvan, who passed it on to Manu; and Manu, in turn, instructed it to Ikshvaku.";


int len_Tx_Signal_repeated = 0, data_frames_number = 0, Data_Frame_Size = 0; // data_frames_number = number of payloads

int len_RRC_Coeff = 21, len_RRC_rx = 10;

double complex *Data = NULL; // This is where all the binary data is stored

double complex **Data_Payload_Mod = NULL;

double complex RRC_Filter_Tx[21] = {-0.000454720514876223, 0.00353689555574986, -0.00714560809091226, 0.00757906190517828, 0.00214368242727367, -0.0106106866672496, 0.0300115539818315, -0.0530534333362480, -0.0750288849545787, 0.409168714634052, 0.803738600397980, 0.409168714634052, -0.0750288849545787, -0.0530534333362480, 0.0300115539818315, -0.0106106866672496, 0.00214368242727367, 0.00757906190517828, -0.00714560809091226, 0.00353689555574986, -0.000454720514876223};

double complex Long_preamble_slot_Frequency[N_FFT]; // Size = N_FFT

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper Functions for Debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

void Display(double complex *X, int sz)
{
    for (int i = 0; i < sz; i++) 
    {
        printf("%d  %lf + %lf*I \n",i+1, creal(X[i]), cimag(X[i]));
    }
    printf("\n");
}

void Display_Double(double* X, int sz)
{
    for (int i = 0; i < sz; i++) 
    {
        printf("%d  %lf \n" , i+1, X[i]);
    }
    printf("\n");
}

void write_complex_array_to_file(double complex* A, int len_A, char* fname) 
{
    FILE* file = fopen(fname, "w");  
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    for (int i = 0; i < len_A; i++) {
        if (!isinf(creal(A[i])) && !isinf(cimag(A[i])))
        {
            //fprintf(file, "%.15e + %.15ei", creal(A[i]), cimag(A[i]));

            fprintf(file, "%.15e", creal(A[i]));
            if (i < len_A - 1) 
            {
                fprintf(file, "\t");  // Tab separator
            }
        }
        else 
        {
            printf("\n%d", i + 1);
        }
    }
    fprintf(file, "\n");

    fclose(file);
}

void write_double_array_to_file(double * A, int len_A, char* fname) 
{
    FILE* file = fopen(fname, "w");  
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    for (int i = 0; i < len_A; i++) 
    {

        fprintf(file, "%.2e", A[i]);
        if (i < len_A - 1) 
        {
            fprintf(file, "\t");  // Tab separator
        }
    }
    fprintf(file, "\n");

    fclose(file);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helper functions for stack/heap memory allocation, fft, ifft, convolution and  slicing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

double complex* Allocate_Array_1D(int size) 
{
    double complex *array = (double complex *)calloc(size, sizeof(double complex));
    if (array == NULL) {
        printf("Memory allocation failed for 1D array!\n");
        exit(1);
    }
    return array;
}

// Function to allocate a 2D array of double complex
double complex** Allocate_Array_2D(int rows, int cols) {
    double complex **array = (double complex **)malloc(rows * sizeof(double complex *));
    if (array == NULL) {
        printf("Memory allocation failed for 2D array (row pointers)!\n");
        exit(1);
    }

    for (int i = 0; i < rows; i++) {
        array[i] = (double complex *)calloc(cols, sizeof(double complex));
        if (array[i] == NULL) {
            printf("Memory allocation failed for 2D array (row %d)!\n", i);
            exit(1);
        }
    }
    return array;
}

void Deallocate_Array_2D(double complex **array, int rows) 
{
    if (array == NULL) 
        return;

    for (int i = 0; i < rows; i++) 
    {
        if (array[i] != NULL) 
        {
            free(array[i]); 
            array[i] = NULL;
        }
    }

    free(array);
    array = NULL; 
}

void Slice_Repeater(double complex *X, double complex *Y, int starti, int lo, int hi, int r) // ind is starting index of Y i.e Output, lo and hi are the starting and ending indices of X(input), r is the number of times we are repeating
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

// Cyclically Shifts frequency domain array before discrete ifft
void ifft_shift(double complex *X, int sz) {
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

// Discete Inverse Fourier Transform : O(N^2)
// void ifft(double complex *X, double complex *Y, int sz) // X = Frequency Domain, Y = Time Domain
// {
//     ifft_shift(X, sz); 

//     double complex W;

//     for (int i = 0; i < sz; ++i)
//     {
//         Y[i] = 0;
//         for (int j = 0; j < sz; ++j)
//         {
//             W = cexp((I * 2.0 * PI * i * j) / sz); 
//             Y[i] += X[j] * W;
//         }
//         Y[i] /= sz; 
//     }
// }

// void fft(double complex *X, double complex *Y, int sz) // X = Time Domain, Y = Frequency Domain
// {
//     double complex W;

//     for (int i = 0; i < sz; ++i)
//     {
//         Y[i] = 0;
//         for (int j = 0; j < sz; ++j)
//         {
//             W = cexp((-I * 2.0 * PI * j * i) / sz);
//             Y[i] += X[j] * W;
//         }
//     }

//     fft_shift(Y, sz); 
// }

void fft(double complex *X, double complex *Y, int sz) // X = Time Domain, Y = Frequency Domain
{
    if (sz <= 1) {
        Y[0] = X[0];
        return;
    }

    double complex *X_even = malloc(sz / 2 * sizeof(double complex));
    double complex *X_odd = malloc(sz / 2 * sizeof(double complex));
    double complex *Y_even = malloc(sz / 2 * sizeof(double complex));
    double complex *Y_odd = malloc(sz / 2 * sizeof(double complex));

    for (int i = 0; i < sz / 2; i++) {
        X_even[i] = X[i * 2];
        X_odd[i] = X[i * 2 + 1];
    }

    fft(X_even, Y_even, sz / 2);
    fft(X_odd, Y_odd, sz / 2);

    for (int k = 0; k < sz / 2; k++) {
        double complex W = cexp(-I * 2.0 * PI * k / sz) * Y_odd[k];
        Y[k] = Y_even[k] + W;
        Y[k + sz / 2] = Y_even[k] - W;
    }

    free(X_even);
    free(X_odd);
    free(Y_even);
    free(Y_odd);
}

void ifft(double complex *X, double complex *Y, int sz) // X = Frequency Domain, Y = Time Domain
{
    double complex *X_conj = malloc(sz * sizeof(double complex));
    double complex *Y_tmp = malloc(sz * sizeof(double complex));

    for (int i = 0; i < sz; i++) {
        X_conj[i] = conj(X[i]);
    }

    fft(X_conj, Y_tmp, sz);

    for (int i = 0; i < sz; i++) {
        Y[i] = conj(Y_tmp[i]) / sz;
    }

    free(X_conj);
    free(Y_tmp);
}


// double complex* Convolution (double complex *Inp, double complex *H, int len_Inp, int len_H)
// {
//     int len_Out = len_Inp + len_H - 1;

//     double complex *Inp_Padded = Allocate_Array_1D(len_Out);
//     Slice_Repeater(Inp, Inp_Padded, 0, 0, len_Inp, 1);
//     //for (int i = len_Inp; i < len_Out; i++) Inp_Padded[i] = 0;

//     double complex *H_Padded = Allocate_Array_1D(len_Out);
//     Slice_Repeater(H, H_Padded, 0, 0, len_H, 1);
//     //for (int i = len_H; i < len_Out; i++) H_Padded[i] = 0;

//     double complex *Inp_freq = Allocate_Array_1D(len_Out);
//     double complex *H_freq = Allocate_Array_1D(len_Out);
//     double complex *Out_freq = Allocate_Array_1D(len_Out);

//     fft(Inp_Padded, Inp_freq, len_Out);
//     fft(H_Padded, H_freq, len_Out);

//     for (int i = 0; i < len_Out; i++) 
//     {
//         Out_freq[i] = Inp_freq[i] * H_freq[i];
//     }

//     double complex *Out = Allocate_Array_1D(len_Out);

//     ifft(Out_freq, Out, len_Out);
    
//     free(Inp_Padded);
//     free(H_Padded);
//     free(Inp_freq);
//     free(H_freq);
//     free(Out_freq);

//     return Out;
// }


double complex* Convolution(double complex *Inp, double complex *H, int len_Inp, int len_H) {
    int len_Out = len_Inp + len_H - 1;  // Output length

    double complex *Out = Allocate_Array_1D(len_Out);

  
    for (int i = 0; i < len_Out; i++) 
    {
        Out[i] = 0.0 + 0.0 * I;
    }

    for (int i = 0; i < len_Inp; i++) 
    {
        for (int j = 0; j < len_H; j++) 
        {
            Out[i + j] += Inp[i] * H[j];  
        }
    }

    return Out; 
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transmitter and it's Processing Blocks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

void Preamble_Generator(double scale, double complex *P_k, double complex *virtual_subcarrier, double complex *Preamble, int type)
{
    //double complex Preamble[160];
    // Scaling S_k
    for (int i = 0; i < 53; i++) 
    {
        P_k[i] *= scale;
    }    

    double complex preamble_freq[N_FFT] = {0};
    double complex preamble_time[N_FFT] = {0};

    Slice_Repeater(virtual_subcarrier, preamble_freq, 0, 0, 6, 1);
    Slice_Repeater(P_k, preamble_freq, 6, 0, 54, 1);
    Slice_Repeater(virtual_subcarrier, preamble_freq, 59, 6, 11, 1);
    
    if(type == 1) // Long Preamble SLot Frequency is used in receiver
    {
        Slice_Repeater(preamble_freq, Long_preamble_slot_Frequency, 0, 0, N_FFT, 1);
    }

    ifft(preamble_freq, preamble_time, N_FFT);


    if(type == 0) // Short Preamble
        Slice_Repeater(preamble_time, Preamble, 0, 0, 16, 10);
    else//long preamble
    {
        Slice_Repeater(preamble_time, Preamble, 0, 32, N_FFT, 1);
        Slice_Repeater(preamble_time, Preamble, 32, 0, N_FFT, 2);
    }
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

void QPSK_Modulator(double complex **Data_Payload, double complex **Data_Payload_Mod, int data_frames_number)
{
    for(int i = 0; i < data_frames_number; ++i)
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

void Data_Generator()
{
    int message_size = sizeof(message)/sizeof(message[0]) - 1;

    data_frames_number = ceil(8 * message_size / 96.0); 

    int total_bits = data_frames_number * 96;

    int total_chars = total_bits / 8;

    Data = Allocate_Array_1D(total_bits);

    // Encoder 
    int msg_bits[8] = {0};

    int index = 0;

    for(int i = 0; i < total_chars; ++i) //convert each character to 8bit binary and store in data array
    {
        if(i < message_size)
            Decimal_To_Binary(message[i], msg_bits);
        else 
            Decimal_To_Binary(' ', msg_bits);

        for(int j = 0; j < 8; ++j)
        {
            Data[index] = msg_bits[j];
            index++;
        }
    }    
}

double complex* Transmitter() // Functions Used: Data_Generator(), Preamble_Generator(), IFFT(), QPSK_Modulator(), Slice_Repeater(), Convoulation()
{
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data Generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    Data_Generator(); // This function initalizes Data array, data frames number and data frames size 


    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Preamble Generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

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

    Preamble_Generator(scale, S_k, virtual_subcarrier, Short_Preamble, 0); //0 is short, 1 is long

    double complex L_k[53] ={1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1};


    scale = 1;
    Preamble_Generator(scale, L_k, virtual_subcarrier, Long_Preamble, 1);    

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Converting Data to Data Payload %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    double complex **Data_Payload = Allocate_Array_2D(data_frames_number, 96);

    for(int i = 0; i < data_frames_number; ++i)
    {
        int row_start = i * 96;
        for(int j = 0; j < 96; ++j)
        {
            Data_Payload[i][j] = Data[row_start + j];
        }
    }

   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QPSK Modulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    Data_Payload_Mod = Allocate_Array_2D(data_frames_number, 48);

    QPSK_Modulator(Data_Payload, Data_Payload_Mod, data_frames_number);

    Deallocate_Array_2D(Data_Payload, data_frames_number);    

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data Frames Generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    int pilot[] = {1,1,1,-1};
    double complex **Data_Frames = Allocate_Array_2D(data_frames_number, 64);

    for(int i = 0; i < data_frames_number; ++i)
    {
        Slice_Repeater(virtual_subcarrier, Data_Frames[i], 0, 0, 6, 1); // Filled = 6

        Slice_Repeater(Data_Payload_Mod[i], Data_Frames[i], 6, 0, 5, 1);    // Filled = 11
        Data_Frames[i][11] = pilot[0];                                 // Filled = 12

        Slice_Repeater(Data_Payload_Mod[i], Data_Frames[i], 12, 5, 18, 1);  // Filled = 25
        Data_Frames[i][25] = pilot[1];                                 // Filled = 26

        Slice_Repeater(Data_Payload_Mod[i], Data_Frames[i], 26, 18, 24, 1); // Filled = 32
        
        Data_Frames[i][32] = 0;                                        // Filled = 33

        Slice_Repeater(Data_Payload_Mod[i], Data_Frames[i], 33, 24, 30, 1); // Filled  = 39
        Data_Frames[i][39] = pilot[2];                                 // FIlled = 40

        Slice_Repeater(Data_Payload_Mod[i], Data_Frames[i], 40, 30, 43, 1); // Filled  = 53
        Data_Frames[i][53] = pilot[3];                                 // Filled  = 54

        Slice_Repeater(Data_Payload_Mod[i], Data_Frames[i], 54, 43, 48, 1); // Filled  = 59
        Slice_Repeater(virtual_subcarrier, Data_Frames[i], 59, 6, 11, 1); // Filled = 64    
    }

    double complex **Data_Frames_IFFT = Allocate_Array_2D(data_frames_number, 64);

    for(int i = 0; i < data_frames_number; ++i)
    {
        ifft(Data_Frames[i], Data_Frames_IFFT[i], 64);
    }

    Deallocate_Array_2D(Data_Frames, data_frames_number);

    double complex **Data_Frames_Time_TX = Allocate_Array_2D(data_frames_number, 80);

    for(int i = 0; i < data_frames_number; ++i)
    {
        Slice_Repeater(Data_Frames_IFFT[i], Data_Frames_Time_TX[i], 0, 48, 64, 1);
        Slice_Repeater(Data_Frames_IFFT[i], Data_Frames_Time_TX[i], 16, 0, 64, 1);
    }

    Deallocate_Array_2D(Data_Frames_IFFT, data_frames_number);

    Data_Frame_Size = data_frames_number * 80 + 160 + 160;
    double complex *Data_Frame_TX = Allocate_Array_1D(Data_Frame_Size);

    Slice_Repeater(Short_Preamble, Data_Frame_TX, 0, 0, 160, 1);   // Filled = 160
    Slice_Repeater(Long_Preamble, Data_Frame_TX, 160, 0, 160, 1);  // Filled = 320

    int start = 320;

    for(int i = 0; i < data_frames_number; ++i)
    {
        Slice_Repeater(Data_Frames_Time_TX[i], Data_Frame_TX, start, 0, 80, 1);
        start += 80;
    }

    Deallocate_Array_2D(Data_Frames_Time_TX, data_frames_number);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    int oversampling_rate_tx = 2; 

    double complex *Data_Frame_TX_Oversamp = Allocate_Array_1D(oversampling_rate_tx * Data_Frame_Size);

    for(int i = 0; i < Data_Frame_Size; ++i)
    {
        Data_Frame_TX_Oversamp[2*i] = Data_Frame_TX[i];
        Data_Frame_TX_Oversamp[2*i + 1] = 0;
    }  

    free(Data_Frame_TX);
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data Frame Oversampled convoulation with RRC FIlter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    int len_inp_sig = oversampling_rate_tx * Data_Frame_Size;
  
    int len_out_sig = len_inp_sig + len_RRC_Coeff - 1; 

    double complex *Tx_signal = Convolution(Data_Frame_TX_Oversamp, RRC_Filter_Tx, len_inp_sig, len_RRC_Coeff);

    len_Tx_Signal_repeated = 10 * len_out_sig;
    double complex *Tx_signal_repeated = Allocate_Array_1D(len_Tx_Signal_repeated);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Repeating Tx Signal 10 times %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Slice_Repeater(Tx_signal, Tx_signal_repeated, 0, 0, len_out_sig, 10);

    free(Tx_signal);
    free(Data_Frame_TX_Oversamp);

    return Tx_signal_repeated;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Transmission over the air and it's Processing Blocks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

double gaussian_noise(double mean, double variance) 
{
    double u1, u2, z0;
    
    u1 = ((double) rand() + 1.0) / ((double) RAND_MAX + 1.0);
    u2 = ((double) rand() + 1.0) / ((double) RAND_MAX + 1.0);

    z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159 * u2);

    return mean + sqrt(variance) * z0;
}


void Transmission_Over_Air(double complex* TX_signal, complex double* TX_OTA_signal, double snr, int len_Tx_Signal)
{
    double Tx_signal_power = 0.0;
    for (int i = 0; i < len_Tx_Signal; i++) 
    {
        Tx_signal_power += (cabs(TX_signal[i])*cabs(TX_signal[i]));
    }

    Tx_signal_power /= len_Tx_Signal; 

    double snr_linear = pow(10, snr / 10);

    double noise_power = Tx_signal_power / snr_linear;

    for(int i = 0 ; i < len_Tx_Signal; ++i)
    {
        double noise = sqrt(noise_power) * (gaussian_noise(0, 1) + I*gaussian_noise(0, 1));

        TX_OTA_signal[i] = TX_signal[i] + noise;
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Receiver and it's Processing Blocks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

double complex* Packet_Detection(double complex *Rx_Signal, int len_RX_Signal, int* len_Corr_Out)
{
    int delay_param= 16;
    int window_length = 32;
    int lenCorr_arr = len_RX_Signal - delay_param + 1 - window_length;

    double complex* corr_out = Allocate_Array_1D(lenCorr_arr);

    for (int i = 0;i < lenCorr_arr; ++i)
    {
        double complex corr_temp_arr = 0;
        double complex peak_temp_arr = 0;

        for (int k = 0; k < window_length; ++k)
        {
            corr_temp_arr += Rx_Signal[i + k] * Rx_Signal[i + k + delay_param];
            peak_temp_arr += cabs(Rx_Signal[i + k + delay_param]) *  cabs(Rx_Signal[i + k + delay_param]);
        }  
        corr_out[i] = ( cabs(corr_temp_arr) * cabs(corr_temp_arr) ) / (peak_temp_arr * peak_temp_arr);
    }

    *len_Corr_Out = lenCorr_arr;

    return corr_out;   
}

int Packet_Selection(double complex* Corr_Out, int len_Corr_Out)
{
    double packet_threshold = 0.75;
    double complex *packet_idx_arr = Allocate_Array_1D(len_Corr_Out);
    int idx_count = 0;

    for (int i=0;i<len_Corr_Out;++i)
    {
        if(cabs(Corr_Out[i])>packet_threshold)
        {
            packet_idx_arr[idx_count]=i;
            idx_count += 1;
        }

    }

    double complex *packet_idx_arr_sliced = Allocate_Array_1D(idx_count);

    Slice_Repeater(packet_idx_arr,packet_idx_arr_sliced,0,0,idx_count,1);


    double complex *temp = Allocate_Array_1D(idx_count+1);

    int a = 0, b = 0;

    for(int i = 0; i < idx_count + 1; ++i)
    {
        if(i < idx_count)
            a = packet_idx_arr_sliced[i];
        else 
            a = -1;

        if(i - 1 >= 0)
            b = packet_idx_arr_sliced[i - 1];
        else 
            b = -1;

        temp[i] = a-b;
    }

    double complex *packet_front = Allocate_Array_1D(idx_count+1);
    int packet_front_count = 0;

    for(int j=0;j<idx_count+1;++j)
    {
        if(creal(temp[j])>300)
        {
            packet_front[packet_front_count] = j;
            packet_front_count += 1;

        }
    }


    double complex *packet_front_sliced = Allocate_Array_1D(packet_front_count);

    Slice_Repeater(packet_front,packet_front_sliced,0,0,packet_front_count,1);

    double complex *packet_front_idx = Allocate_Array_1D(packet_front_count);

    for (int i = 0;i<packet_front_count;++i)
    {
        packet_front_idx[i] = packet_idx_arr_sliced[(int)creal(packet_front_sliced[i])];

    }


    int threshold_length=230, packet_idx = 0;

    for (int x=0;x<packet_front_count-1;++x)
    {
        if(cabs(Corr_Out[(int)creal(packet_front_idx[x])+230])>packet_threshold)
        {
            packet_idx = packet_front_idx[x] + len_RRC_rx+1;
            break;
        }
    }

    free(packet_idx_arr);
    free(packet_idx_arr_sliced);
    free(temp);
    free(packet_front);
    free(packet_front_sliced);
    free(packet_front_idx);

    return packet_idx;
}

void Coarse_CFO_Estimation(double complex* rx_frame, double complex* rx_frame_after_coarse, int rx_frame_size)
{
    int Short_preamble_slot_length = 16;
    
    double complex vec1[Short_preamble_slot_length];
    double complex vec2[Short_preamble_slot_length];
    
    int start1 = Short_preamble_slot_length * 5;
    for (int i = 0; i < Short_preamble_slot_length; i++) 
    {
        vec1[i] = rx_frame[start1 + i];
    }
    
    int start2 = Short_preamble_slot_length * 6;
    for (int i = 0; i < Short_preamble_slot_length; i++) 
    {
        vec2[i] = rx_frame[start2 + i];
    }

    double complex prod_consq_frame_coarse = 0;
    for (int i = 0; i < Short_preamble_slot_length; i++) 
    {
        prod_consq_frame_coarse += vec1[i] * conj(vec2[i]);
    }

    double freq_coarse_est = (-1.0 / (2 * PI * Short_preamble_slot_length * ts_sec)) *atan2(cimag(prod_consq_frame_coarse), creal(prod_consq_frame_coarse));
    
    for (int i = 0; i < rx_frame_size; i++) 
    {
        rx_frame_after_coarse[i] = rx_frame[i] * cexp(-I * 2 * PI * freq_coarse_est * ts_sec * i);
    }
}

void Fine_CFO_Estimation(double complex* rx_frame_after_coarse, double complex* rx_frame_after_fine, int rx_frame_size)
{
    int Short_preamble_slot_length = 16;

    int start1 = Short_preamble_slot_length * 12;
    int start2 = Short_preamble_slot_length * 16;
    int length= Short_preamble_slot_length * 4;

    double complex prod_consq_frame_fine = 0 + 0 * I;

    for (int i = 0; i < length; i++) 
    {
        prod_consq_frame_fine += rx_frame_after_coarse[start1 + i] * conj(rx_frame_after_coarse[start2 + i]);
    }

    double freq_fine_est = (-1.0 / (2 * PI * 64 * ts_sec)) * atan2(cimag(prod_consq_frame_fine), creal(prod_consq_frame_fine));

    for (int i = 0; i < rx_frame_size; i++) 
    {
        double complex exp_term = cexp(-1.0 * I * 2 * PI * freq_fine_est * ts_sec * i);
        rx_frame_after_fine[i] = rx_frame_after_coarse[i] * exp_term;
    }    
}

void Channel_Estimation(double complex* rx_frame_after_fine, double complex* H_est, int rx_frame_size)
{
    int Short_preamble_slot_length = 16;

    double complex* Long_preamble_1 = Allocate_Array_1D(N_FFT);
    double complex* Long_preamble_2 = Allocate_Array_1D(N_FFT);

    Slice_Repeater(rx_frame_after_fine, Long_preamble_1, 0, Short_preamble_slot_length * 12, Short_preamble_slot_length * 16, 1);
    Slice_Repeater(rx_frame_after_fine, Long_preamble_2, 0, Short_preamble_slot_length * 16, Short_preamble_slot_length * 20, 1);

    double complex *Long_preamble_1_After_FFT = Allocate_Array_1D(N_FFT);
    double complex *Long_preamble_2_After_FFT = Allocate_Array_1D(N_FFT);
    
    fft(Long_preamble_1, Long_preamble_1_After_FFT, N_FFT);
    fft(Long_preamble_2, Long_preamble_2_After_FFT, N_FFT);

    for (int i = 0; i < N_FFT; i++) 
    {
        H_est[i] = 0.5 * (Long_preamble_1_After_FFT[i] + Long_preamble_2_After_FFT[i]) * conj(Long_preamble_slot_Frequency[i]);
    }
}

void AGC_Receiver(double complex** Rx_Payload_No_Pilot, double complex** Rx_Payload_Final)
{
    for(int i = 0; i < data_frames_number; ++i)
    {
       for(int j = 0; j < 48; ++j)
       {
           double complex z = Rx_Payload_No_Pilot[i][j];

           if(creal(z) > 0)
               Rx_Payload_Final[i][j] = 1.0/sqrt(2.0);
           else 
               Rx_Payload_Final[i][j] = -1.0/sqrt(2.0);

           if(cimag(z) > 0)
               Rx_Payload_Final[i][j] += I*1.0/sqrt(2.0);
           else 
               Rx_Payload_Final[i][j] += I*-1.0/sqrt(2.0);
       }
    }
}

void QPSK_Demodulator(double complex **Rx_Payload, double complex **Data_Payload_Demod, int data_frames_number)
{
    for(int i = 0; i < data_frames_number; ++i)
    {
        for(int j = 0; j < 48; j++)
        {
            double a = creal(Rx_Payload[i][j]), b = cimag(Rx_Payload[i][j]);

            int c = 0, d = 0;

            if(a > 0 && b > 0)
            {
                c = 0;
                d = 0;
            }
            else if(a < 0 && b > 0)
            {
                c = 0;
                d = 1;
            }
            else if(a < 0 && b < 0)
            {
                c = 1;
                d = 0;
            }
            else 
            {
                c = 1;
                d = 1;
            }
            
            Data_Payload_Demod[i][2 * j] = c;
            Data_Payload_Demod[i][2 * j + 1] = d;
        }
    }
}

int Binary_To_Decimal(int* msg_bits) 
{
    int res = 0;

    for(int i = 0; i < 8; ++i)
    {
        res = (res * 2) + msg_bits[i];
    }    
    return res;
}

void Message_Generator(double complex* Data, char* msg_RX, int total_bits)
{
    int index = 0;

    int msg_bits[8];

    int msg_index = 0;
    for(int i = 0; i < total_bits; i += 8)
    {
        for(int j = 0; j < 8; ++j)
        {
            msg_bits[j] = Data[i + j];
        }

        // Decoder
        
        msg_RX[msg_index++] = Binary_To_Decimal(msg_bits);
    }
}

void Receiver(double complex* Tx_OTA_signal, int len_Tx_Signal, int data_frames_number, double* Res)
{
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Capturing Packets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    int len_Rx_Signal = len_Tx_Signal * 0.307; // Number of packets Captured

    // len_Rx_Signal = 3000;  // To be removed

    int rx_start = rand() % (len_Tx_Signal - len_Rx_Signal);

    // rx_start = 0; // To be removed

    double complex* Rx_Signal = Allocate_Array_1D(len_Rx_Signal);

    Slice_Repeater(Tx_OTA_signal, Rx_Signal, 0, rx_start, rx_start + len_Rx_Signal, 1);

    // Padding Signal

    int len_inp_sig = len_Rx_Signal;
  
    int len_out_sig = len_inp_sig + len_RRC_Coeff - 1; 


    // Convolution
    double complex *Rx_filter_signal = Convolution(Rx_Signal, RRC_Filter_Tx, len_inp_sig, len_RRC_Coeff);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Packet Detection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    int len_Corr_Out = 0;

    double complex *Corr_Out = Packet_Detection(Rx_Signal, len_Rx_Signal, &len_Corr_Out);

    free(Rx_Signal);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Packet Selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    int packet_idx = Packet_Selection(Corr_Out, len_Corr_Out);

    free(Corr_Out);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Down Sampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    int oversampling_rate = 2;

    int rx_frame_size = ((oversampling_rate * Data_Frame_Size + packet_idx - 1) - packet_idx) / oversampling_rate + 1;

    double complex* rx_frame = Allocate_Array_1D(rx_frame_size);

    int index = 0;

    for(int i = packet_idx; i < oversampling_rate * Data_Frame_Size + packet_idx - 1; i += oversampling_rate)
    {
        rx_frame[index] = Rx_filter_signal[i];
        index += 1;
    }

    free(Rx_filter_signal);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Coarse CFO Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    double complex* rx_frame_after_coarse = Allocate_Array_1D(rx_frame_size);

    Coarse_CFO_Estimation(rx_frame, rx_frame_after_coarse, rx_frame_size);

    free(rx_frame);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fine CFO Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    double complex* rx_frame_after_fine = Allocate_Array_1D(rx_frame_size);

    Fine_CFO_Estimation(rx_frame_after_coarse, rx_frame_after_fine, rx_frame_size);

    free(rx_frame_after_coarse);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Channel Estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    double complex* H_est = Allocate_Array_1D(64); // Size = Short Preamble Slot Length * 4 amd Short Preamble Slot Length = 16

    Channel_Estimation(rx_frame_after_fine, H_est, rx_frame_size);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% One Tap Equalizer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    double complex** Rx_Payload_Time = Allocate_Array_2D(data_frames_number, N_FFT); // Extracting Rx_Payload_Time and removing circular prefix from it in this step

    for(int i = 0; i < data_frames_number; ++i)
    {
        int lo = 320 + i * 80 + 16;
        int hi = 320 + (i + 1) * 80;
        Slice_Repeater(rx_frame_after_fine, Rx_Payload_Time[i], 0, lo, hi, 1);
    }

    free(rx_frame_after_fine);

    double complex** Rx_Payload_Frequency = Allocate_Array_2D(data_frames_number, N_FFT);

    for(int i = 0; i < data_frames_number; ++i)
    {
        fft(Rx_Payload_Time[i], Rx_Payload_Frequency[i], N_FFT);
    }

    Deallocate_Array_2D(Rx_Payload_Time, data_frames_number);

    double complex** Rx_Payload_Frequency_Equalizer = Allocate_Array_2D(data_frames_number, N_FFT);

    for(int i = 0; i < data_frames_number; ++i)
    {
        for(int j = 0; j < N_FFT; ++j)
        {
            Rx_Payload_Frequency_Equalizer[i][j] = Rx_Payload_Frequency[i][j] / H_est[j];
        }
    }

    free(H_est);
    Deallocate_Array_2D(Rx_Payload_Frequency, data_frames_number);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% De Mapping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    double complex** Rx_Payload_No_Pilot = Allocate_Array_2D(data_frames_number, 48);

    for(int i = 0; i < data_frames_number; ++i)
    {
        Slice_Repeater(Rx_Payload_Frequency_Equalizer[i], Rx_Payload_No_Pilot[i], 0, 6, 11, 1);
        Slice_Repeater(Rx_Payload_Frequency_Equalizer[i], Rx_Payload_No_Pilot[i], 5, 12, 25, 1);
        Slice_Repeater(Rx_Payload_Frequency_Equalizer[i], Rx_Payload_No_Pilot[i], 18, 26, 32, 1);
        Slice_Repeater(Rx_Payload_Frequency_Equalizer[i], Rx_Payload_No_Pilot[i], 24, 33, 39, 1);
        Slice_Repeater(Rx_Payload_Frequency_Equalizer[i], Rx_Payload_No_Pilot[i], 30, 40, 53, 1);
        Slice_Repeater(Rx_Payload_Frequency_Equalizer[i], Rx_Payload_No_Pilot[i], 43, 54, 59, 1);
    }

    Deallocate_Array_2D(Rx_Payload_Frequency_Equalizer, data_frames_number);

     //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AGC for RX_Data_Payload %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    double complex** Rx_Payload_Final = Allocate_Array_2D(data_frames_number, 48);

    AGC_Receiver(Rx_Payload_No_Pilot, Rx_Payload_Final);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QPSK Demodulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    double complex** Rx_Payload_Demod = Allocate_Array_2D(data_frames_number, 96);

    QPSK_Demodulator(Rx_Payload_Final, Rx_Payload_Demod, data_frames_number);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Combining all Data Frames %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    int total_bits = data_frames_number * 96;

    double complex* Data_Rx = Allocate_Array_1D(total_bits);

    index = 0;

    for(int i = 0; i < data_frames_number; ++i)
    {
        for(int j = 0; j < 96; ++j)
        {
            Data_Rx[index] = Rx_Payload_Demod[i][j];
            index++;
        }
    }

    Deallocate_Array_2D(Rx_Payload_Demod, data_frames_number);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EVM Calculation (Before AGC) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    double complex error = 0;
    
    double error_square_sum = 0, data_payload_square_sum = 0;
    
    for(int i = 0; i < data_frames_number; ++i)
    {
        for(int j = 0; j < 48; ++j)
        {
            error = Rx_Payload_No_Pilot[i][j] - Data_Payload_Mod[i][j];
            error_square_sum += pow(cabs(error), 2);
            data_payload_square_sum += pow(cabs(Data_Payload_Mod[i][j]), 2);
        }
    }

    Deallocate_Array_2D(Rx_Payload_No_Pilot, data_frames_number);

    double evm = 0, evm_dB = 0;

    evm = sqrt(error_square_sum / (data_frames_number * 48)) / sqrt(data_payload_square_sum / (data_frames_number * 48));

    evm_dB = 20 * log10(evm);   

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EVM Calculation (After AGC) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    error = 0;
    error_square_sum = 0;
    data_payload_square_sum = 0;
    
    for(int i = 0; i < data_frames_number; ++i)
    {
        for(int j = 0; j < 48; ++j)
        {
            error = Rx_Payload_Final[i][j] - Data_Payload_Mod[i][j];
            error_square_sum += pow(cabs(error), 2);
            data_payload_square_sum += pow(cabs(Data_Payload_Mod[i][j]), 2);
        }
    }

    Deallocate_Array_2D(Rx_Payload_Final, data_frames_number);

    double evm_AGC = 0, evm_AGC_dB = 0;

    evm_AGC = sqrt(error_square_sum / (data_frames_number * 48)) / sqrt(data_payload_square_sum / (data_frames_number * 48));

    evm_AGC_dB = 20 * log10(evm_AGC);   

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BER Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    double sum = 0, ber = 0;

    for(int i = 0; i < total_bits; ++i)
    {
        sum += abs(creal(Data[i]) - creal(Data_Rx[i]));
    }

    ber = sum / total_bits;
    
    Res[0] = evm_dB;
    Res[1] = evm_AGC_dB;
    Res[2] = ber;

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Converting Bits to Message %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    int msg_Rx_size = total_bits / 8;

    char msg_RX[msg_Rx_size];

    Message_Generator(Data_Rx, msg_RX, total_bits);

    free(Data_Rx);

    printf("\nReceived Message: \n");
    for(int i = 0; i < msg_Rx_size; ++i)
    {
        printf("%c", msg_RX[i]);
    }
    printf("\n");


}

int main() 
{
    srand(time(NULL));

    double complex* TX_signal_repeated = Transmitter();

    double SNR[num_snr];

    for(int i = 0; i < num_snr; ++i)
    {
        SNR[i] = 1 + i;
    }

    double EVM_dB[num_snr], EVM_AGC_dB[num_snr], BER[num_snr];

    for(int i = 0; i < num_snr; ++i)
    {
        printf("\n\nFor SNR = %lf \n", SNR[i]);

        double complex* Tx_OTA_signal = Allocate_Array_1D(len_Tx_Signal_repeated);  

        Transmission_Over_Air(TX_signal_repeated, Tx_OTA_signal, SNR[i], len_Tx_Signal_repeated);
        
        double Res[3]; // EVM_dB, EVM_AGC_DB, BER
        Receiver(Tx_OTA_signal, len_Tx_Signal_repeated, data_frames_number, Res);

        free(Tx_OTA_signal);

        EVM_dB[i] = Res[0];
        EVM_AGC_dB[i] = Res[1];
        BER[i] = Res[2];   
        
        printf("\nEVM dB     = %lf", Res[0]);
        printf("\nEVM_ACG dB = %lf", Res[1]);
        printf("\nBER        = %lf", Res[2]);
    }

    free(TX_signal_repeated);
    free(Data);
    Deallocate_Array_2D(Data_Payload_Mod, data_frames_number);

    write_double_array_to_file(SNR, num_snr, "Output_SNR.txt");
    write_double_array_to_file(EVM_dB, num_snr, "Output_EVM_AGC.txt");
    write_double_array_to_file(EVM_AGC_dB, num_snr, "Output_EVM_AGC_DB.txt");
    write_double_array_to_file(BER, num_snr, "Output_BER.txt");

    printf("\nCode Run Successful!");

    return 0;    
}