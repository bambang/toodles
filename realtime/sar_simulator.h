#define PI 3.14159265
#define C 300000000
//#define MEMORY_SIZE 104857600
#define MEMORY_SIZE 209715200
//#define MEMORY_SIZE 429496730
//#define MEMORY_SIZE 536870912
//#define MEMORY_SIZE 1073741824
//#define MEMORY_SIZE 2147483648
//
//

#include <complex.h>

typedef struct{
  char real_or_complex;
  unsigned int rows;
  unsigned int cols;
}radar_metadata;

void chirp_generator();
void chirp_matched_generator();
void insert_waveform_in_scene();
void radar_imager();
void pulse_compress_image();
void gbp();
void gbp_fft();
void pulse_compress_signal();
void fft_waveform(unsigned int kernel_length, double complex* kernel, double complex* output);
void filter_dc();
int write_complex_data();
int write_real_data();
int simulate();
void process_data();
int read_radar_file();
