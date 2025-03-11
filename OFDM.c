#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>
#include <complex.h>  // Required for complex numbers
#include <time.h>

#define N_FFT 64

#define PI 3.14159265358979323846

int len_Tx_Signal = 0, data_frames_number = 0, Data_Frame_Size = 0; // Used for Up Sampling and Down Sampling

// data_frames_number = number of payloads

complex double *Data = NULL; // This is where all the binary data is stored

double RRC_Filter_Tx[21] = {-0.000454720514876223, 0.00353689555574986, -0.00714560809091226, 0.00757906190517828, 0.00214368242727367, -0.0106106866672496, 0.0300115539818315, -0.0530534333362480, -0.0750288849545787, 0.409168714634052, 0.803738600397980, 0.409168714634052, -0.0750288849545787, -0.0530534333362480, 0.0300115539818315, -0.0106106866672496, 0.00214368242727367, 0.00757906190517828, -0.00714560809091226, 0.00353689555574986, -0.000454720514876223};

int len_RRC_Coeff = 21;
int len_RRC_rx = 10;

// Cyclically Shifts frequency domain array before discrete ifft
void ifft_shift(complex double *X) 
{
    int mid = N_FFT / 2;

    for(int i = 0; i < mid; ++i)
    {
        complex double temp = X[i];
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

void Convolution (complex double *Out, double *impulse, complex double *In, int lenOut, int lenImp)
{
    for(int i = 0; i < lenOut; ++i)
    {
        for(int j = 0; j < lenImp; ++j)
        {
            Out[i] += impulse[j] * In[i + j - 1];
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

    Slice_Repeater(virtual_subcarrier, preamble_freq, 0, 0, 6, 1);
    Slice_Repeater(P_k, preamble_freq, 6, 0, 54, 1);
    Slice_Repeater(virtual_subcarrier, preamble_freq, 59, 6, 11, 1);
    

    // printf("Preamble (Frequency Domain):\n");
    // Display(preamble_freq, N_FFT);

    ifft(preamble_freq, preamble_time);

    // printf("Preamble Shift (Time Domain):\n");
    // Display(preamble_time, N_FFT);

    if(type == 0) // Short Preamble
        Slice_Repeater(preamble_time, Preamble, 0, 0, 16, 10);
    else//long preamble
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
complex double** Allocate_Array_2D(int data_frames_number, int cols) {
    complex double **array = (complex double **)malloc(data_frames_number * sizeof(complex double *));
    if (array == NULL) {
        printf("Memory allocation failed for 2D array (row pointers)!\n");
        exit(1);
    }

    for (int i = 0; i < data_frames_number; i++) {
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

void QPSK_Modulator(complex double **Data_Payload, complex double **Data_Payload_Mod, int data_frames_number)
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
    unsigned char message[] = "AAAAAAAAAAAAAAAAAAAAAAAA";
    //"The Supreme Lord Shree Krishna said: I taught this eternal science of Yog to the Sun God, Vivasvan, who passed it on to Manu; and Manu, in turn, instructed it to Ikshvaku.";

    int message_size = sizeof(message)/sizeof(message[0]) - 1;

    data_frames_number = ceil(8 * message_size / 96.0); 

    int total_bits = data_frames_number * 96;

    Data = Allocate_Array_1D(total_bits);

    if (Data == NULL) 
    {
        printf("Memory allocation failed!\n");
        return;
    }

    // Encoder 
    int msg_bits[8] = {0};

    int index = 0;

    for(int i = 0; i < message_size; ++i) //convert each character to 8bit binary and store in data array
    {
        Decimal_To_Binary(message[i], msg_bits);

        for(int j = 0; j < 8; ++j)
        {
            Data[index] = msg_bits[j];
            index++;
        }
    }    
}

complex double *Transmitter() // Functions Used: Data_Generator(), Preamble_Generator(), IFFT(), QPSK_Modulator(), Slice_Repeater(), Convoulation()
{
    Data_Generator(); // This function initalizes Data array, data frames number and data frames size 

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

    // Preparing Payload

    int M = 4; // QPSK modulation
    int bits_per_symb = log2(M); // Bits per symbol =2 (have to make pairs later)
    int n_bits = 96; // Number of bits in one frame

    

    // Creating Data Payloads

    complex double **Data_Payload = Allocate_Array_2D(data_frames_number, 96);

    for(int i = 0; i < data_frames_number; ++i)
    {
        int row_start = i * 96;
        for(int j = 0; j < 96; ++j)
        {
            Data_Payload[i][j] = Data[row_start + j];
        }
    }

    // Combining consecutive bits for QPSK modulation: Creating 1 X 48 complex numbers from 1 * 2 real numbers

    complex double **Data_Payload_Mod = Allocate_Array_2D(data_frames_number, 48);

    QPSK_Modulator(Data_Payload, Data_Payload_Mod, data_frames_number);

    // Data Frames Generation

    int pilot[] = {1,1,1,-1};
    complex double **Data_Frames = Allocate_Array_2D(data_frames_number, 64);

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

    complex double **Data_Frames_IFFT = Allocate_Array_2D(data_frames_number, 64);

    for(int i = 0; i < data_frames_number; ++i)
    {
        ifft(Data_Frames[i], Data_Frames_IFFT[i]);
    }

    complex double **Data_Frames_Time_TX = Allocate_Array_2D(data_frames_number, 80);

    for(int i = 0; i < data_frames_number; ++i)
    {
        Slice_Repeater(Data_Frames_IFFT[i], Data_Frames_Time_TX[i], 0, 48, 64, 1);
        Slice_Repeater(Data_Frames_IFFT[i], Data_Frames_Time_TX[i], 16, 0, 64, 1);
    }

    Data_Frame_Size = data_frames_number * 80 + 160 + 160;
    complex double *Data_Frame_TX = Allocate_Array_1D(Data_Frame_Size);

    Slice_Repeater(Short_Preamble, Data_Frame_TX, 0, 0, 160, 1);   // Filled = 160
    Slice_Repeater(Long_Preamble, Data_Frame_TX, 160, 0, 160, 1);  // Filled = 320

    int start = 320;

    for(int i = 0; i < data_frames_number; ++i)
    {
        Slice_Repeater(Data_Frames_Time_TX[i], Data_Frame_TX, start, 0, 80, 1);
        start += 80;
    }

    // Oversampling 

    int oversampling_rate_tx = 2; 

    complex double *Data_Frame_TX_Oversamp = Allocate_Array_1D(oversampling_rate_tx * Data_Frame_Size);

    for(int i = 0; i < Data_Frame_Size; ++i)
    {
        Data_Frame_TX_Oversamp[2*i] = Data_Frame_TX[i];
        Data_Frame_TX_Oversamp[2*i + 1] = 0;
    }  
    
    // Convoulation of the Tx Frame with RRC Filter

    int len_inp_sig = 2*Data_Frame_Size;
  
    int len_out_sig = len_inp_sig + len_RRC_Coeff - 1; 

    complex double *inp_signal_padded = Allocate_Array_1D(len_inp_sig + 2 * (len_RRC_Coeff - 1));

    Slice_Repeater(Data_Frame_TX_Oversamp, inp_signal_padded, len_RRC_Coeff - 1, 0, len_inp_sig, 1);

    complex double *Tx_signal = Allocate_Array_1D(len_out_sig);

    Convolution(Tx_signal,RRC_Filter_Tx,inp_signal_padded, len_out_sig, len_RRC_Coeff);

    Display(Tx_signal, len_out_sig);

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

    z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * 3.14159 * u2);

    return mean + sqrt(variance) * z0;
}


complex double* Transmission_Over_Air(complex double* TX_signal, double snr, int len_Tx_Signal)
{
    complex double* TX_OTA_signal = Allocate_Array_1D(len_Tx_Signal); 

    double Tx_signal_power = 0.0;
    for (int i = 0; i < len_Tx_Signal; i++) 
    {
        Tx_signal_power += (cabs(TX_signal[i])*cabs(TX_signal[i]));
    }
    Tx_signal_power /= len_Tx_Signal; 

    double snr_linear = pow(10, Tx_signal_power / 10);

    double noise_power = Tx_signal_power / snr_linear;

    for(int i = 0 ; i < len_Tx_Signal; ++i)
    {
        TX_OTA_signal[i] += sqrt(noise_power) * (gaussian_noise(0, 1) + I*gaussian_noise(0, 1));
    }

    return TX_OTA_signal;
}

complex double* Packet_Detection(complex double *Rx_Signal, int len_RX_Signal, int* len_Corr_Out)
{

    int delay_param= 16;
    int window_length =32;
    int lenCorr_arr = len_RX_Signal-delay_param+1-window_length;

    complex double* corr_out = Allocate_Array_1D(lenCorr_arr);

    complex double corr_temp_arr = 0;
    complex double peak_temp_arr = 0;

    for (int i=0;i<lenCorr_arr;++i)
    {
        for (int k=0;k<window_length;++k)
        {
            corr_temp_arr += Rx_Signal[i+k-1]*Rx_Signal[i+k-1+delay_param];
            peak_temp_arr += abs(Rx_Signal[i+k-1+delay_param]*Rx_Signal[i+k-1+delay_param]);
        }  
        corr_out[i] = ( abs(corr_temp_arr) * abs(corr_temp_arr) ) / ( abs(peak_temp_arr) * abs(peak_temp_arr) );
    }

    *len_Corr_Out = lenCorr_arr;
    return corr_out;   
}

int Packet_Selection(complex double* Corr_Out, int len_Corr_Out)
{
    double packet_threshold = 0.75;
    complex double *packet_idx_arr = Allocate_Array_1D(len_Corr_Out);
    int idx_count = 0;

    for (int i=0;i<len_Corr_Out;++i)
    {
        if(abs(Corr_Out[i])>packet_threshold)
        {
            idx_count+=1;
            packet_idx_arr[idx_count]=i;
        }

    }
    Slice_Repeater(packet_idx_arr,packet_idx_arr,0,0,idx_count,1);
    complex double *temp1 = Allocate_Array_1D(idx_count+1);
    complex double *temp2 = Allocate_Array_1D(idx_count+1);
    complex double *temp3 = Allocate_Array_1D(idx_count+1);

    Slice_Repeater(packet_idx_arr,temp1,0,0,idx_count,1);
    Slice_Repeater(packet_idx_arr,temp2,1,0,idx_count,1);

    for(int j=0;j<idx_count+1;++j)
    {
        temp3[j] = temp1[j] - temp2[j];
    }

    complex double *packet_front = Allocate_Array_1D(idx_count+1);
    int packet_front_count =0;

    for(int j=0;j<idx_count+1;++j)
    {
        if(abs(temp3[j])>300)
        {
            packet_front_count+=1;
            packet_front[packet_front_count] =j;


        }
    }

    Slice_Repeater(packet_front,packet_front,0,0,packet_front_count,1);
    complex double *packet_front_idx = Allocate_Array_1D(packet_front_count);

    for (int i =0;i<packet_front_count;++i)
    {
        packet_front_idx[i] = packet_idx_arr[abs(packet_front[i])];

    }


    int threshold_length=230, packet_idx = 0;

    for (int x=0;x<packet_front_count-1;++x)
    {
        if(abs(Corr_Out[abs(packet_front_idx[x])+230])>threshold_length)
        {
            packet_idx = packet_front_idx[x] + len_RRC_rx+1;
            break;
        }
    }

    return packet_idx;
}

void Receiver(complex double* Tx_OTA_signal, int len_Tx_Signal, int data_frames_number)
{
    // Capturing Packets

    int len_Rx_Signal = len_Tx_Signal * 0.307; // Number of packets Captured

    int rx_start = rand() % (len_Tx_Signal - len_Rx_Signal);

    complex double* Rx_Signal = Allocate_Array_1D(len_Rx_Signal);

    Slice_Repeater(Tx_OTA_signal, Rx_Signal, 0, 0, len_Rx_Signal, 1);

    // Padding Signal

    int len_inp_sig = len_Rx_Signal;
  
    int len_out_sig = len_inp_sig + len_RRC_Coeff - 1; 

    complex double *inp_signal_padded = Allocate_Array_1D(len_inp_sig + 2 * (len_RRC_Coeff - 1));

    Slice_Repeater(Rx_Signal, inp_signal_padded, len_RRC_Coeff - 1, 0, len_inp_sig, 1);

    complex double *Rx_filter_signal = Allocate_Array_1D(len_out_sig);

    // Convoulation

    Convolution(Rx_filter_signal,RRC_Filter_Tx,inp_signal_padded, len_out_sig, len_RRC_Coeff);

    // Packet Detection

    int len_Corr_Out = 0;

    complex double *Corr_Out = Packet_Detection(Rx_Signal, len_out_sig, &len_Corr_Out);

    // Packet Selection

    int packet_idx = Packet_Selection(Corr_Out, len_Corr_Out);


    // Down Sampling

    int oversampling_rate = 2;

    complex double* rx_frame = Allocate_Array_1D(oversampling_rate * Data_Frame_Size + packet_idx - 1);

    int index = 0;

    for(int i = packet_idx; i < oversampling_rate * Data_Frame_Size + packet_idx - 1; i += oversampling_rate)
    {
        rx_frame[index] = Rx_filter_signal[i];
    } 
    
}


int main() 
{
    srand(time(NULL));

    complex double* TX_signal = Transmitter();

    // Display(TX_signal, len_Tx_Signal);

    double SNR[5] = {100, 101, 102, 103, 104};

    printf("\n\nCode Run Successful!");

    for(int i = 0; i < 1; ++i)
    {
        complex double* Tx_OTA_signal = Transmission_Over_Air(TX_signal, SNR[i], len_Tx_Signal);

        Receiver(Tx_OTA_signal, len_Tx_Signal, data_frames_number);
    }
    return 0;    
}