/* 
 * File:   L4_DelayDiversity.c
 * Author: kopic
 *
 */

/*
 * Course: Signal Processing for MIMO Communications
 * Lab4: Delay Diversity
 */

/*
  COMPILE with:
    gcc -o L4_DelayDiversity.exe L4_DelayDiversity.c -lm -Wall
  Usage example:
    ./L4_DelayDiversity.exe 1000000 2
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265359
#define SQRT_OF_2 1.41421356237
      
#define SCALE_QPSK 0.707106781186547 
#define QPSK 2
     
#define FILE_NAME_SIZE 20
#define FILE_BER "BER_DD"
    
#define fDT 4e-5 // 43.2 km/h -> 12 m/s and T = 1us
#define BEAMS 16
     
void usage(char* progName);
double rand_ra();
void generate_info_bits(int len, unsigned char* bits);
void print_info_bits(unsigned char* bits, int len);

void generate_symbol_qpsk(unsigned char* bits, double* symbol_I, double* symbol_Q);
void print_symbols(double* symbols, int len);

void ddEncoder(unsigned char* txBits, double* txSymI, double* txSymQ,
        int bitsPerSymbol, double* sqrtSNR, int j);

void box_muller(double* real, double* imag);
void add_noise(double* sI, double* sQ, int len);

void rayleigh_fading_init(double* alpha, double* phi, int rxAntennas, int j,
        int howManyBeams, double fadingFactor);
void rayleigh_fading(double* alpha, double* phi, int antIndex, int howManyBeams, 
        int symbolIndex, double* hi, double* hq);

void received_signal(double txSymI, double txSymQ, double* rxSymI, double* rxSymQ, 
        double* alpha, double* phi, double* hi, double* hq, int rxAntennas, int j, double current_sqrt_EsN0);
void decode_qpsk(double recSymI, double recSymQ, 
        unsigned char* recBit, double sqrtSNR);
void ddDecoder(double* rxSymI, double* rxSymQ, double* rxComI, double* rxComQ, 
        double* hi, double* hq, int rxAntennas, int j);

void fix_file_name(char* fileName, char* rxAntennas, char* kernel);
int save_asignal(double* t, double* signal, int len, char* fileName);  

int main(int argc, char** argv) {

    /*
     * Initialize the random generator. Must be done once per run of main
     */
    srand((unsigned) time(NULL));

    if (argc < 3) {
        usage(argv[0]);
        return 0;
    }

    int bitsPerSymbol = QPSK;
    double EbN0[] = {0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25};
    int howManySymbols = atoi(argv[1]);
    int howManyRxAntennas = atoi(argv[2]); 

    // Reuse to generate info bits at every second iteration
    unsigned char txBits[bitsPerSymbol];

    unsigned char rxBits[bitsPerSymbol];

    // SNR per bit
    double* snr = (double*) malloc(sizeof(EbN0)); 
    double* sqrtSNR = (double*) malloc(sizeof(EbN0));
    for (int i = 0; i < sizeof(EbN0) / sizeof(double); i++) {
        snr[i] = pow( 10, EbN0[i] / 10);
        
        /* 
         * To normalize energy per bit.
         * Models "modulation loss".
         * Makes comparison of different modulations fair.
         */ 
        snr[i] *= bitsPerSymbol;
        sqrtSNR[i] = sqrt(snr[i]);
    }
    
    // Tx symbol 
    double txSymI, txSymQ;
    
    // Rx symbols before decoding
    double* rxSymI = (double*) calloc(howManyRxAntennas * 2, sizeof(double));
    double* rxSymQ = (double*) calloc(howManyRxAntennas * 2, sizeof(double));
    
    // Rx symbol after decoding
    double rxComI = 0, rxComQ = 0;
    
    double* alpha = (double*) malloc( howManyRxAntennas * 2 * BEAMS * sizeof(double));
    double* phi = (double*) malloc( howManyRxAntennas * 2 * BEAMS * sizeof(double));

    double* ber = (double*) malloc(sizeof(EbN0));
    
    double* hi = (double*) calloc(howManyRxAntennas * 2, sizeof(double));
    double* hq = (double*) calloc(howManyRxAntennas * 2, sizeof(double));
    
    for (int i = 0; i < sizeof(EbN0) / sizeof(double); i++) {
        
        ber[i] = 0;
        rayleigh_fading_init(alpha, phi, howManyRxAntennas, 0, BEAMS, fDT);
        rayleigh_fading_init(alpha, phi, howManyRxAntennas, 1, BEAMS, fDT);
        for (int j = 0; j < howManySymbols * 2; j++) {
            
            // Delay encoder
            ddEncoder(txBits, &txSymI, &txSymQ, bitsPerSymbol, sqrtSNR + i, j); 

            // Channel
            received_signal(txSymI, txSymQ, rxSymI, rxSymQ, alpha, phi, hi, hq, 
                    howManyRxAntennas, j, sqrtSNR[i]);
            
            // Delay decoder
            ddDecoder(rxSymI, rxSymQ, &rxComI, &rxComQ, hi, hq, howManyRxAntennas, j);   
            

            if (j % 2 == 1) {       
                      
                decode_qpsk(rxComI, rxComQ, rxBits, sqrtSNR[i]);              
                
                for (int k = 0; k < bitsPerSymbol; k++) {
                    if (txBits[k] != rxBits[k]) ber[i]++;    
                }
                rxComI = 0, rxComQ = 0;
            }
        }

        ber[i] /= (howManySymbols * bitsPerSymbol);
        printf("%f\t%f\n", EbN0[i], ber[i]);
    }
    
    char fileNameBER[FILE_NAME_SIZE];  
    fix_file_name(fileNameBER, argv[2], FILE_BER);
    
    save_asignal(EbN0, ber, sizeof(EbN0)/sizeof(double), fileNameBER);
        
    free(rxSymI); free(rxSymQ);
    free(snr); free(sqrtSNR);
    free(alpha); free(phi);
    free(hi); free(hq);
    free(ber);

    return (EXIT_SUCCESS);
}


void usage(char* progName) {
    printf("Usage: %s <number of symbols> <number of receive antennas>\n", progName);
}

/*
 * rand_ra (): Returns a random number from 0 to 1 
 */
double rand_ra() {
    return ((double) rand()) / RAND_MAX;
}

void generate_info_bits(int len, unsigned char* bits) {
    double tmp;
    for (int k = 0; k < len; k++) {
        tmp = rand_ra();
        bits[k] = (tmp > 0.5) ? 1 : 0;
    }
}

void print_info_bits(unsigned char* bits, int len) {
    printf("Generated info bits:\n");
    for (int j = 0; j < len; j++) {
        printf("%u", bits[j]); 
    }
    printf("\n");    
}

void generate_symbol_qpsk(unsigned char* bits, double* symbol_I, double* symbol_Q) {
    
        *symbol_I = (bits[0]) ? 1 : -1;
        *symbol_I *= SCALE_QPSK;

        *symbol_Q = (bits[1]) ? 1 : -1;
        *symbol_Q *= SCALE_QPSK;
}

void print_symbols(double* symbols, int len) {
    for (int i = 0; i < len; i++) {
        printf("%3d", (int) symbols[i]);
    }
    printf("\n");
}

void ddEncoder(unsigned char* txBits, double* txSymI, double* txSymQ,
        int bitsPerSymbol, double* sqrtSNR, int j) {
    /*
    * Fill in the function
    * Output: a symbol to be transmitted from transmit antenna
    */ 

     static double stored_sI = 0.0, stored_sQ = 0.0; 
    
    if (j%2 == 0) {
        // Even slot: generate new bits and symbol 
        generate_info_bits(bitsPerSymbol, txBits); 
        // Copy generated bits to prevBits 
     
        generate_symbol_qpsk(txBits, &stored_sI, &stored_sQ);

        *txSymI = stored_sI; /// (SQRT_OF_2);
        *txSymQ = stored_sQ; // / (SQRT_OF_2);
        // printf("Transmitted symbol from ANT 1: I = %.3f, Q = %.3f\n", 
        //         *txSymI, *txSymQ);
      
    } else {
        // Odd slot: use previous bits and symbol
        *txSymI = stored_sI; /// (SQRT_OF_2); 
        *txSymQ = stored_sQ; // (SQRT_OF_2); 
        // printf("Transmitted symbol from ANT 2: I = %.3f, Q = %.3f\n", 
        //        *txSymI, *txSymQ);
    }
}

void box_muller(double* real, double* imag) {
    double uniform1, uniform2;

    do {
        uniform1 = rand_ra();
    } while (uniform1 == 0);

    uniform2 = rand_ra();

    *real = sqrt(-2 * log(uniform1)) * cos(2 * PI * uniform2);
    *imag = sqrt(-2 * log(uniform1)) * sin(2 * PI * uniform2);
}

void add_noise(double* sI, double* sQ, int len) {
    double noiseI, noiseQ;

    for (int i = 0; i < len; i++) {
        box_muller(&noiseI, &noiseQ);

        sI[i] += noiseI;
        sQ[i] += noiseQ;
    }
}

void rayleigh_fading_init(double* alpha, double* phi, int rxAntennas, int j,
        int howManyBeams, double fadingFactor) {

    double x, y;
    for (int i = 0; i < rxAntennas; i++) {
        x = rand_ra();
        for (int k = 0; k < howManyBeams; k++) {
            int baseIdx = (i * 2 * howManyBeams) + (j * howManyBeams);
            alpha[baseIdx + k] = 2 * PI 
                    * fadingFactor * cos(2 * PI * (k - x) / howManyBeams);
            y = rand_ra();
            phi[baseIdx + k] = 2 * PI * y;
        }
    }
}

void rayleigh_fading(double* alpha, double* phi, int antIndex, int howManyBeams, 
        int symbolIndex, double* hi, double* hq) {
    /*
    * Lab 2
    */
    double current_hi_sum = 0.0; 
    double current_hq_sum = 0.0;
   
    for (int j = 0; j < howManyBeams; j++) {
        double phase = alpha[j] * symbolIndex + phi[j];
        current_hi_sum += cos(phase);
        current_hq_sum += sin(phase);
    }
    double norm = 1.0 / sqrt((double)howManyBeams);

    *hi = current_hi_sum * norm;
    *hq = current_hq_sum * norm;          
}

void received_signal(double txSymI, double txSymQ, double* rxSymI, double* rxSymQ, 
        double* alpha, double* phi, double* hi, double* hq, int rxAntennas, int j, double current_sqrt_EsN0) {
    /*
    * Fill in the function.
    * Output: Rx symbols after fading the signal and adding noise
    */

    int active_tx_slot_idx = j % 2; // 0 if Tx1 is active (for s_t), 1 if Tx2 is active (for s_{t-T})

    for (int r = 0; r < rxAntennas; r++) {
            // Channel from TxAnt 1 to RxAnt r (h_r0)
        //pointer to the start of alpha parameters for the link RxAnt r, TxAnt 0; 
        double* current_link_alpha_tx0 = alpha + (r * 2 * BEAMS) + (active_tx_slot_idx * BEAMS);
        //pointer to the start of phi parameters for the link RxAnt r, TxAnt 0;
        double* current_link_phi_tx0 = phi + (r * 2 * BEAMS) + (active_tx_slot_idx * BEAMS);

        // Fading for the first transmit antenna
        rayleigh_fading(current_link_alpha_tx0, current_link_phi_tx0, r, BEAMS, 
                j, &hi[r * 2 + 0], &hq[r * 2 + 0]);  
    }

   for (int r=0; r < rxAntennas; r++) {
        double h_active_I = hi[r*2 + active_tx_slot_idx];
        double h_active_Q = hq[r*2 + active_tx_slot_idx];
        rxSymI[2*r + active_tx_slot_idx] = h_active_I * txSymI - h_active_Q * txSymQ;
        rxSymQ[2*r + active_tx_slot_idx] = h_active_I * txSymQ + h_active_Q * txSymI;
        
    }

    if (current_sqrt_EsN0 < 1e-9) {
        printf("Warning: current_sqrt_EsN0 is too small, setting sigma_n to a large value.\n");
        current_sqrt_EsN0 = 1e-9; // Avoid division by zero
    }
    double sigma_n = 1.0 / current_sqrt_EsN0; 
    // Add noise (AWGN)

    for (int r = 0; r < rxAntennas; r++) {

        double noiseI_sample, noiseQ_sample;
        box_muller(&noiseI_sample, &noiseQ_sample); // Generate fresh noise for each component

        // Add scaled noise to the received signal for the current active path
        rxSymI[2*r + active_tx_slot_idx] += noiseI_sample * sigma_n;
        rxSymQ[2*r + active_tx_slot_idx] += noiseQ_sample * sigma_n;
    }



}

void decode_qpsk(double recSymI, double recSymQ, 
        unsigned char* recBits, double sqrtSNR) {
    
    /*
     * Does not really matter since compared to 0, 
     * but keeps everything consistent 
     */
    /*
     * Noise has power 2\sigma^2, so the signal is scaled with sqrt(2)
     * We revert it here
     */ 
    // recSymI /= (sqrtSNR * SQRT_OF_2);
    // recSymQ /= (sqrtSNR * SQRT_OF_2);
    
    recBits[0] = (recSymI > 0) ? 1 : 0;
    recBits[1] = (recSymQ > 0) ? 1 : 0;
}

void ddDecoder(double* rxSymI, double* rxSymQ, double* rxComI, double* rxComQ, 
        double* hi, double* hq, int rxAntennas, int j) {
    /*
    * Fill in the function.
    * Output: Rx symbol after decoding
    */

       
    double s_hat_I_accumulator = 0.0;
    double s_hat_Q_accumulator = 0.0;
    for (int r = 0; r < rxAntennas; r++) {
        // h1: channel for Tx1's transmission of s (from previous even slot, j-1)
        double h1_I = hi[r*2 + 0]; 
        double h1_Q = hq[r*2 + 0]; 

        // h2: channel for Tx2's transmission of s (from current odd slot, j)
        double h2_I = hi[r*2 + 1]; 
        double h2_Q = hq[r*2 + 1]; 

        // y1: received signal from Tx1's transmission (from previous even slot, j-1)
        double y1_I  = rxSymI[r*2 + 0]; 
        double y1_Q  = rxSymQ[r*2 + 0];

        // y2: received signal from Tx2's transmission (from current odd slot, j)
        double y2_I = rxSymI[r*2 + 1];
        double y2_Q = rxSymQ[r*2 + 1];

        // MRC: z_r = y1 * h1_conj + y2 * h2_conj for receive antenna r
        // Re(y1*h1_conj) = y1_I*h1_I + y1_Q*h1_Q
        // Im(y1*h1_conj) = y1_Q*h1_I - y1_I*h1_Q
        
        // Re(y2*h2_conj) = y2_I*h2_I + y2_Q*h2_Q
        // Im(y2*h2_conj) = y2_Q*h2_I - y2_I*h2_Q

        s_hat_I_accumulator += (y1_I * h1_I + y1_Q * h1_Q) + (y2_I * h2_I + y2_Q * h2_Q);
        s_hat_Q_accumulator += (y1_Q * h1_I - y1_I * h1_Q) + (y2_Q * h2_I - y2_I * h2_Q);
    }
    *rxComI = s_hat_I_accumulator;
    *rxComQ = s_hat_Q_accumulator;
    // if (tx_idx == 0) {
    //     for (int r = 0; r < rxAntennas; r++) {
    //         double h_r_I_T = hi[r*2 + 0]; 
    //         double h_r_Q_T = hq[r*2 + 0]; 

    //         double h_r_I_T1 = hi[r*2 + 1]; 
    //         double h_r_Q_T1 = hq[r*2 + 1]; 

    //         double y_r_t_I_T  = rxSymI[r*2 + 0]; 
    //         double y_r_t_Q_T  = rxSymQ[r*2 + 0];

    //         double y_r_t_I_T1 = rxSymI[r*2 + 1];
    //         double y_r_t_Q_T1 = rxSymQ[r*2 + 1];


    //         s_hat_I_accumulator += (h_r_I_T * y_r_t_I_T) + (h_r_Q_T * y_r_t_Q_T) + (h_r_I_T1 * y_r_t_I_T1) + (h_r_Q_T1 * y_r_t_Q_T1);
    //         s_hat_Q_accumulator += (h_r_I_T * y_r_t_Q_T) - (h_r_Q_T * y_r_t_I_T) + (h_r_I_T1 * y_r_t_Q_T1) - (h_r_Q_T1 * y_r_t_I_T1);
    //     }
    //     *rxComQ = s_hat_Q_accumulator;
    //     *rxComI = s_hat_I_accumulator;
    // }
         
}

void fix_file_name(char* fileName, char* rxAntennas, char* kernel) {

    char tmp[FILE_NAME_SIZE];
    strcpy(tmp, kernel);   

    strcat(tmp, "_Nr");
    strcat(tmp, rxAntennas);
    strcat(tmp, ".txt");

    strcpy(fileName, tmp);
}

int save_asignal(double* t, double* signal, int len, char* fileName) {
    FILE* fptr = fopen(fileName, "w");

    for (int k = 0; k < len; k++) {
        fprintf(fptr, "%.9f\t%.9f\t\n", t[k], signal[k]);
    }
    fclose(fptr);

    return 0;
}
