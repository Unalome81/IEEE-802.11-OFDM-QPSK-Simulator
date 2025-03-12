Hello, this is the implementation of IEEE 802.11 Wireless Communication System as part of our AELD Project.

Transmitter has been completed. Receiver has been partially completed.

Line 1-250 are complete and debug has been done.

The project members are:

Vivaswan Nawani
Krishna Ayyagari
Nikhil Kumar
Nishita Bharat Lohana

File Description:

1) OFDM.c : This file contains the main implementation of   IEEE 802.11 protocol in C
2) Code_Output.txt: This file contains the output written from OFDM.c to this file using function write_complex_array_to_file()
3) IEEE_802_11_a_Code: This is the implementation of IEEE 802.11 in MATLAB
4) IEEE_802_11_a_Code_Tester: This is used to debug OFDM.c. Instead of a random or a very large input, it contains 

    Message = "AAAAAAAAAAAAAAAAAAAAAAA"
    Noise = 0
    starting index for capturing 3000 packets in receiver = 0

5) Matlab_Output.txt : This file contains output of a complex double or a double array from IEEE_802_11_a_Code_Tester to compare with corresponding OFDM.c values.
6) compare_complex.py : This class is used to find mismatches and for error analysis when there are complex double values in Code_Output.txt and Matlab_Output.txt
7) compare_double.py : This class is used to find mismatches and for error analysis when there are double values in Code_Output.txt and Matlab_Output.txt

Also, right now we are checking functionality for SNR = 100 only, once that works well, we will check it in the range 8:15 with the final goal of transmitting a reasonably long text message (200-250 characters) and receiving it.