/* 
 * File:   L4_PhaseRollDiversity.c
 * Author: kopic
 *
 */

/*
 * Course: Signal Processing for MIMO Communications
 * Lab4: Phase-roll Diversity
 */

/*
  COMPILE with:
    gcc -o L4_PhaseRollDiversity.exe L4_PhaseRollDiversity.c -lm -Wall
  Usage example:
    ./L4_PhaseRollDiversity.exe 1000000 2
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define PI 3.14159265359
#define SQRT_OF_2 1.41421356237
#define THETA 1/4
      
#define SCALE_QPSK 0.707106781186547 
#define QPSK 2
     
#define FILE_NAME_SIZE 20
#define FILE_BER "BER_PR"
    
#define fDT 4e-5 // 43.2 km/h -> 12 m/s and T = 1us
#define BEAMS 16
     
void usage(char* progName);
double rand_ra();
void generate_info_bits(int len, unsigned char* bits);
void print_info_bits(unsigned char* bits, int len);

void generate_symbol_qpsk(unsigned char* bits, double* symbol_I, double* symbol_Q);
void print_symbols(double* symbols, int len);

void prEncoder(unsigned char* txBits, double* txSymI1, 
        double* txSymQ1, double* txSymI2, double* txSymQ2,
        int bitsPerSymbol, double* sqrtSNR, int j);

void box_muller(double* real, double* imag);
void add_noise(double* sI, double* sQ, int len);

void rayleigh_fading_init(double* alpha, double* phi, int rxAntennas,
        int howManyBeams, double fadingFactor);
void rayleigh_fading(double* alpha, double* phi, int antIndex, int howManyBeams, 
        int symbolIndex, double* hi, double* hq);

void received_signal(double txSymI1, double txSymQ1, double txSymI2, double txSymQ2, 
        double* rxSymI, double* rxSymQ, double* alpha, double* phi, 
        double* hi, double* hq, int rxAntennas, int j);
void decode_qpsk(double recSymI, double recSymQ, 
        unsigned char* recBit, double sqrtSNR);
void prDecoder(double* rxSymI, double* rxSymQ, double* rxComI, double* rxComQ, 
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

    // To generate info bits
    unsigned char txBits[bitsPerSymbol];

    // To decode Rx symbol
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
    
    // Tx symbols (there are 2 transmit antennas)
    double txSymI1, txSymQ1, txSymI2, txSymQ2;
    
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
        // Initialize Rayleigh fading
        rayleigh_fading_init(alpha, phi, howManyRxAntennas, BEAMS, fDT);
        for (int j = 0; j < howManySymbols * 2; j++) {
            
            // Phase-roll encoder
            prEncoder(txBits, &txSymI1, &txSymQ1, &txSymI2, &txSymQ2, 
                    bitsPerSymbol, sqrtSNR + i, j);

            // Channel
            received_signal(txSymI1, txSymQ1, txSymI2, txSymQ2, rxSymI, rxSymQ, 
                    alpha, phi, hi, hq, howManyRxAntennas, j);
            
            // Phase-roll decoder
            prDecoder(rxSymI, rxSymQ, &rxComI, &rxComQ, hi, hq, howManyRxAntennas, j);

              
            decode_qpsk(rxComI, rxComQ, rxBits, sqrtSNR[i]);              
            
            for (int k = 0; k < bitsPerSymbol; k++) {
                if (txBits[k] != rxBits[k]) ber[i]++;    
            }
            rxComI = 0, rxComQ = 0;
        
        }

        ber[i] /= (2*howManySymbols * bitsPerSymbol);
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

void prEncoder(unsigned char* txBits, double* txSymI1, double* txSymQ1, 
        double* txSymI2, double* txSymQ2,
        int bitsPerSymbol, double* sqrtSNR, int j) {

    /*
    * Fill in the function
    * Output: symbols to be transmitted from transmit antennas
    */ 

    double s1_I_base, s1_Q_base;
    generate_info_bits(bitsPerSymbol, txBits); 
    generate_symbol_qpsk(txBits, &s1_I_base, &s1_Q_base);

    *txSymI1 = s1_I_base;
    *txSymQ1 = s1_Q_base;
    double phase_roll = 2.0 * PI * THETA * j;
    double cos_phase = cos(phase_roll);
    double sin_phase = sin(phase_roll);

    *txSymI2 = (*txSymI1) * cos_phase - (*txSymQ1) * sin_phase;
    *txSymQ2 = (*txSymI1) * sin_phase + (*txSymQ1) * cos_phase;

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

void rayleigh_fading_init(double* alpha, double* phi, int rxAntennas,
        int howManyBeams, double fadingFactor) {

    double x, y;
    
    for (int i = 0; i < rxAntennas; i++) {
        for (int j = 0; j < 2; j++) {
            x = rand_ra();
            for (int k = 0; k < howManyBeams; k++) {
                alpha[i * 2 * howManyBeams + j * howManyBeams + k] = 2 * PI 
                        * fadingFactor * cos(2 * PI * (k - x) / howManyBeams);
                y = rand_ra();
                phi[i * 2 * howManyBeams + j * howManyBeams + k] = 2 * PI * y;
            }
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

void received_signal(double txSymI1, double txSymQ1, double txSymI2, double txSymQ2, 
        double* rxSymI, double* rxSymQ, 
        double* alpha, double* phi, double* hi, double* hq, int rxAntennas, int j) {

    /*
    * Fill in the function.
    * Output: Rx symbols after fading the signal and adding noise
    */

    for (int r=0; r < rxAntennas; r++) {

        double* current_link_alpha_tx0 = alpha + (r * 2 * BEAMS);
        double* current_link_phi_tx0 = phi + (r * 2 * BEAMS);

        rayleigh_fading(current_link_alpha_tx0, current_link_phi_tx0, r, BEAMS, j, &hi[2*r], &hq[2*r]);

        double* current_link_alpha_tx1 = alpha + (r * 2 * BEAMS) + (1 * BEAMS);
        double* current_link_phi_tx1 = phi + (r * 2 * BEAMS) + (1 * BEAMS);

        rayleigh_fading(current_link_alpha_tx1, current_link_phi_tx1, r, BEAMS, j, &hi[2*r + 1], &hq[2*r + 1]);


        
        rxSymI[r] = (hi[2*r] * txSymI1 - hq[2*r] * txSymQ1) + 
            (hi[2*r + 1] * txSymI2 - hq[2*r + 1] * txSymQ2);
        rxSymQ[r] = (hi[2*r] * txSymQ1 + hq[2*r] * txSymI1) + 
                    (hi[2*r + 1] * txSymQ2 + hq[2*r + 1] * txSymI2);

        
    }
    add_noise(rxSymI, rxSymQ, rxAntennas);

    
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
    recSymI /= (sqrtSNR * SQRT_OF_2);
    recSymQ /= (sqrtSNR * SQRT_OF_2);
    
    recBits[0] = (recSymI > 0) ? 1 : 0;
    recBits[1] = (recSymQ > 0) ? 1 : 0;
}

void prDecoder(double* rxSymI, double* rxSymQ, double* rxComI, double* rxComQ, 
        double* hi, double* hq, int rxAntennas, int j) {
    
    /*
    * Fill in the function.
    * Output: Rx symbol after decoding
    */  

    *rxComI = 0.0;
    *rxComQ = 0.0;

    double conjugate_phase_roll = -2.0 * PI * j * THETA;
    double p_star_I = cos(conjugate_phase_roll);
    double p_star_Q = sin(conjugate_phase_roll);

    for (int m = 0; m < rxAntennas; m++) {
        // Get channel coefficients for Rx antenna m (computed in received_signal)
        // h_1m = (hi[2*m], hq[2*m])
        // h_2m = (hi[2*m+1], hq[2*m+1])
        double h1_i = hi[2 * m];
        double h1_q = hq[2 * m];
        double h2_i = hi[2 * m + 1];
        double h2_q = hq[2 * m + 1];

        // Calculate conjugates: h_1m* and h_2m*
        double h1_star_i = h1_i;
        double h1_star_q = -h1_q;
        double h2_star_i = h2_i;
        double h2_star_q = -h2_q;

        // Calculate Term2 = h_2m* * P*
        double term2_I = (h2_star_i * p_star_I) - (h2_star_q * p_star_Q);
        double term2_Q = (h2_star_i * p_star_Q) + (h2_star_q * p_star_I);

        // Calculate C_m = h_1m* + Term2
        // C_m = h_1m* + h_2m* * e^(-j2πkθ)
        double c_m_I = h1_star_i + term2_I;
        double c_m_Q = h1_star_q + term2_Q;

        // Get received signal y_m for current antenna m
        double y_m_I = rxSymI[m];
        double y_m_Q = rxSymQ[m];

        // Calculate z_m = y_m * C_m
        double z_m_I = (y_m_I * c_m_I) - (y_m_Q * c_m_Q);
        double z_m_Q = (y_m_I * c_m_Q) + (y_m_Q * c_m_I);

        // Accumulate z_m for Maximal Ratio Combining
        *rxComI += z_m_I;
        *rxComQ += z_m_Q;
    }


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

