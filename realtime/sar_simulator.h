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
  double* chirp_time_vector;
  double* matched_time_vector;
  double complex* chirp_fft;
  double complex* matched_fft;
  double complex* pulse_compressed_waveform;
  double complex* chirp_signal;
  double complex* matched_chirp;
  double complex* scene_with_waveform;
  double complex* radar_image;
  double complex* pulse_compressed_radar_image;
  double complex* sar_image;
  double complex* sar_fft;
  double complex* sar_img_shifted;
}data_arrays;

typedef struct{
  char real_or_complex;
  unsigned int rows;
  unsigned int cols;
}radar_metadata;

typedef struct{
  long unsigned int start_frequency;
  long unsigned int bandwidth;
  unsigned int chirp_length;
  unsigned int nrows;
  unsigned int ncols;
  float btproduct;
  int altitude;
  float beamwidth;
  double signal_distance;
  char mode;
  char radar_data_filename[255];
  char real_or_complex_simulation[2];
}radar_variables;

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
float calculate_compressed_pulse_resolution();
void normalize_radar_data();
void free_memory();
