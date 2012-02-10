#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fftw3.h>
#include <sys/time.h>

void chirp_generator(unsigned int start_frequency, unsigned int bandwidth, unsigned int* chirp_length, double complex** chirp_signal, double* signal_distance);
void chirp_matched_generator(unsigned int start_frequency, unsigned int bandwidth, unsigned int* chirp_length, double complex** chirp_signal);
double* scene_generator(unsigned int nrows, unsigned int ncols);
void insert_waveform_in_scene(unsigned int waveform_length, double complex* waveform, unsigned int nrows, unsigned int ncols, double* scene, unsigned int* nnrows, unsigned int* nncols, double complex** scene_with_waveform);
double complex* radar_imager(unsigned int nnrows, unsigned int nncols, double complex* scene, unsigned int altitude, float beamwidth);
void gbp(unsigned int nrows, unsigned int ncols, double complex* raw_data, double complex** sar_image);
double complex* pulse_compressed(unsigned int kernel_length, double complex* kernel, unsigned int nrows, unsigned int ncols, double complex* image);
void gbp_fft(unsigned int nrows, unsigned int ncols, double complex* sar_data, double complex** sar_fft);

#define PI 3.14159265
// 10MB
#define MEMORY_SIZE 10485760

double scene[MEMORY_SIZE/sizeof(double)];
double complex chirp_signal[MEMORY_SIZE/sizeof(double complex)];
double complex matched_chirp[MEMORY_SIZE/sizeof(double complex)];
double complex scene_with_waveform[MEMORY_SIZE/sizeof(double complex)];
double complex radar_image[MEMORY_SIZE/sizeof(double complex)];
double complex pulse_compressed_radar_image[MEMORY_SIZE/sizeof(double complex)];
double complex sar_image[MEMORY_SIZE/sizeof(double complex)];
double complex sar_fft[MEMORY_SIZE/sizeof(double complex)];

int main(int argc, char** argv){
  double complex* chirp_signal;
  double complex* matched_chirp;
  double complex* scene_with_waveform;
  double complex* radar_image;
  double complex* pulse_compressed_radar_image;
  double complex* sar_image;
  double complex* sar_fft;
  double* scene;

  unsigned int chirp_length;
  unsigned int nrows = 0;
  unsigned int ncols = 0;
  unsigned int nnrows, nncols;
  int altitude = 0;
  float beamwidth = 0;
  unsigned int start_frequency = 0;
  unsigned int bandwidth = 0;
  double signal_distance = 0;

  struct timeval otime, ntime;

  if(argc != 7){
    printf("Input: scene_cols scene_rows altitude beamwidth start_frequency bandwidth\n");
    return;
  }
  ncols = atoi(argv[1]);
  nrows = atoi(argv[2]);
  altitude = atoi(argv[3]);
  beamwidth = atoi(argv[4])*PI/180;
  start_frequency = atoi(argv[5]);
  bandwidth = atoi(argv[6]);

  gettimeofday(&otime, NULL);

  scene = scene_generator(nrows, ncols);
  if(scene == 0)
    return;
  
  gettimeofday(&ntime, NULL);
  printf("Scene generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  gettimeofday(&otime, NULL);

  chirp_generator(start_frequency, bandwidth, &chirp_length, &chirp_signal, &signal_distance);
  if(chirp_signal == 0){
    free(scene);
    return;
  }

  gettimeofday(&ntime, NULL);
  printf("Chirp generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  gettimeofday(&otime, NULL);

  chirp_matched_generator(start_frequency, bandwidth, &chirp_length, &matched_chirp);
  if(matched_chirp == 0){
    free(scene);
    free(chirp_signal);
    return;
  }

  gettimeofday(&ntime, NULL);
  printf("Matched chirp generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  gettimeofday(&otime, NULL);

  insert_waveform_in_scene(chirp_length, chirp_signal, nrows, ncols, scene, &nnrows, &nncols, &scene_with_waveform);
  if(scene_with_waveform == 0){
    free(scene);
    free(chirp_signal);
    free(matched_chirp);
    return;
  }

  gettimeofday(&ntime, NULL);
  printf("Scene with waveform generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  gettimeofday(&otime, NULL);

  radar_image = radar_imager(nnrows, nncols, scene_with_waveform, altitude, beamwidth);
  if(radar_image == 0){
    free(scene);
    free(chirp_signal);
    free(matched_chirp);
    free(scene_with_waveform);
    return;
  }

  gettimeofday(&ntime, NULL);
  printf("Radar image generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  gettimeofday(&otime, NULL);

  pulse_compressed_radar_image = pulse_compressed(chirp_length, matched_chirp, nnrows, nncols, radar_image);
  if(pulse_compressed_radar_image == 0){
    free(scene);
    free(chirp_signal);
    free(matched_chirp);
    free(scene_with_waveform);
    free(radar_image);
    return;
  }

  gettimeofday(&ntime, NULL);
  printf("Pulse compression of radar image took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  gettimeofday(&otime, NULL);

  gbp(nnrows, nncols, pulse_compressed_radar_image, &sar_image);
  if(sar_image == 0){
    free(scene);
    free(chirp_signal);
    free(matched_chirp);
    free(scene_with_waveform);
    free(pulse_compressed_radar_image);
    return;
  }

  gettimeofday(&ntime, NULL);
  printf("GBP of pulse-compressed radar image took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  gettimeofday(&otime, NULL);


  gbp_fft(nnrows, nncols, sar_image, &sar_fft);
  if(sar_fft == 0){
    free(scene);
    free(chirp_signal);
    free(matched_chirp);
    free(scene_with_waveform);
    free(pulse_compressed_radar_image);
    free(sar_image);
    return;
  }

  gettimeofday(&ntime, NULL);
  printf("FFT generation of GBP image took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));

  printf("Chirp length: %u\n", chirp_length);
  printf("Number of rows in scene: %i\n", nrows);
  printf("Number of columns in scene: %i\n", ncols);
  printf("Number of rows in final scene: %i\n", nnrows);
  printf("Number of columns in final scene: %i\n", nncols);
  printf("Platform altitude: %i\n", altitude);
  printf("Antenna beamwidth: %f\n", beamwidth);
  printf("Start frequency: %i\n", start_frequency);
  printf("Bandwidth: %i\n", bandwidth);

  FILE* dimensions = fopen("dimensions.dat", "w");
  FILE* scenef = fopen("scene.dat", "w");
  FILE* chirpf = fopen("chirp.dat", "w");
  FILE* matchedf = fopen("matched.dat", "w");
  FILE* scene_with_waveformf = fopen("scene_with_waveform.dat", "w");
  FILE* radar_imagef = fopen("radar_image.dat", "w");
  FILE* pulse_compressedf = fopen("pulse_compressed_image.dat", "w");
  FILE* sar_imagef = fopen("sar_image.dat", "w");
  FILE* sar_fftf = fopen("sar_fft.dat", "w");

  fprintf(dimensions, "%u\n%u\n%u\n%u\n%u\n%f\n", chirp_length, nrows, ncols, nnrows, nncols, signal_distance);

  int i,j;
  for(i = 0; i < ncols; i++){
    for(j = 0; j < nrows; j++){
      fprintf(scenef, "%g\t",scene[i*nrows+j]);
    }
    fprintf(scenef, "\n");
  }

  for(i = 0; i < chirp_length; i++){
    fprintf(chirpf, "%f\t", creal(chirp_signal[i]));
    fprintf(chirpf, "%f\n", cimag(chirp_signal[i]));
  }

  for(i = 0; i < chirp_length; i++){
    fprintf(matchedf, "%f\t", creal(matched_chirp[i]));
    fprintf(matchedf, "%f\n", cimag(matched_chirp[i]));
  }

  for(i = 0; i < nncols; i++){
    for(j = 0; j < nnrows; j++){
      fprintf(scene_with_waveformf, "%g\t", creal(scene_with_waveform[i*nnrows+j]));
      fprintf(scene_with_waveformf, "%g\t", cimag(scene_with_waveform[i*nnrows+j]));
    }
    fprintf(scene_with_waveformf, "\n");
  }

  for(i = 0; i < nncols; i++){
    for(j = 0; j < nnrows; j++){
      fprintf(radar_imagef, "%f\t", creal(radar_image[i*nnrows+j]));
      fprintf(radar_imagef, "%f\t", cimag(radar_image[i*nnrows+j]));
    }
    fprintf(radar_imagef, "\n");
  }

  for(i = 0; i < nncols; i++){
    for(j = 0; j < nnrows; j++){
      fprintf(pulse_compressedf, "%f\t", creal(pulse_compressed_radar_image[i*nnrows+j]));
      fprintf(pulse_compressedf, "%f\t", cimag(pulse_compressed_radar_image[i*nnrows+j]));
    }
    fprintf(pulse_compressedf, "\n");
  }

  for(i = 0; i < nncols; i++){
    for(j = 0; j < nnrows; j++){
      fprintf(sar_imagef, "%f\t", creal(sar_image[i*nnrows+j]));
      fprintf(sar_imagef, "%f\t", cimag(sar_image[i*nnrows+j]));
    }
    fprintf(sar_imagef, "\n");
  }

  for(i = 0; i < nncols; i++){
    for(j = 0; j < nnrows; j++){
      fprintf(sar_fftf, "%f\t", creal(sar_fft[i*nnrows+j]));
      fprintf(sar_fftf, "%f\t", cimag(sar_fft[i*nnrows+j]));
    }
    fprintf(sar_fftf, "\n");
  }

  fclose(sar_fftf);
  fclose(pulse_compressedf);
  fclose(dimensions);
  fclose(scenef);
  fclose(chirpf);
  fclose(matchedf);
  fclose(scene_with_waveformf);
  fclose(radar_imagef);
  fclose(sar_imagef);

  free(sar_fft);
  free(sar_image);
  free(radar_image);
  free(pulse_compressed_radar_image);
  free(scene_with_waveform);
  free(scene);
  free(chirp_signal);
  free(matched_chirp);
}

void chirp_generator(unsigned int start_frequency, unsigned int bandwidth, unsigned int *chirp_length, double complex** chirp_signal, double* signal_distance){
  float chirp_rate = 100;
  float end_time = bandwidth/chirp_rate;
  unsigned int sample_frequency = 5*bandwidth;
  unsigned int time_steps = end_time*sample_frequency;
  *chirp_length = time_steps;

  *signal_distance = end_time*300000000;

  

  double* chirp_time_vector = malloc(end_time*sample_frequency*sizeof(double));
  if(chirp_time_vector == 0){
    printf("Could not allocate memory for time vector.\n");
    return;
  }

  int i;
  float last_time = 0;
  for(i = 0; i < end_time*sample_frequency; i++){
    chirp_time_vector[i] = last_time;
    last_time += (float)1/sample_frequency;
  }

  *chirp_signal = malloc(end_time*sample_frequency*sizeof(double complex));
  if(chirp_signal == 0){
    printf("Could not allocate memory for chirp signal.\n");
    return;
  }

  double complex* signal = *chirp_signal;
  for(i = 0; i < end_time*sample_frequency; i++){
    signal[i] = cos( 2*PI*(start_frequency*chirp_time_vector[i]+chirp_rate*chirp_time_vector[i]*chirp_time_vector[i] ))+_Complex_I*sin( 2*PI*(start_frequency*chirp_time_vector[i]+chirp_rate*chirp_time_vector[i]*chirp_time_vector[i] ));
  }
}

void chirp_matched_generator(unsigned int start_frequency, unsigned int bandwidth, unsigned int *chirp_length, double complex** chirp_signal){
  float chirp_rate = 100;
  float end_time = bandwidth/chirp_rate;
  unsigned int sample_frequency = 5*bandwidth;
  unsigned int time_steps = end_time*sample_frequency;
  *chirp_length = time_steps;

  double* chirp_time_vector = malloc(end_time*sample_frequency*sizeof(double));
  if(chirp_time_vector == 0){
    printf("Could not allocate memory for matched time vector.\n");
    return;
  }
  int i;
  float last_time = end_time;
  for(i = 0; i < end_time*sample_frequency; i++){
    chirp_time_vector[i] = last_time;
    last_time -= (float)1/sample_frequency;
  }

  *chirp_signal = malloc(end_time*sample_frequency*sizeof(double complex));
  if(chirp_signal == 0){
    printf("Could not allocate memory for matched chirp signal.\n");
    return;
  }
  double complex* signal = *chirp_signal;
  for(i = 0; i < end_time*sample_frequency; i++){
    signal[i] = cos( 2*PI*(start_frequency*chirp_time_vector[i]+chirp_rate*chirp_time_vector[i]*chirp_time_vector[i] ))-_Complex_I*sin( 2*PI*(start_frequency*chirp_time_vector[i]+chirp_rate*chirp_time_vector[i]*chirp_time_vector[i] ));
  }
}

double complex* pulse_compressed(unsigned int kernel_length, double complex* kernel, unsigned int nrows, unsigned int ncols, double complex* image){
  double complex* compressed = malloc(nrows*ncols*sizeof(double complex));
  if(compressed == 0){
    printf("Could not allocate memory for pulse-compressed image.\n");
    return 0;
  }

  int i,j,k;
  int index;
  double complex* image_column;
  double complex* compressed_column;
  for(k = 0; k < ncols; k++){
    image_column = &image[k*nrows];
    compressed_column = &compressed[k*nrows];
    for(i = 0; i < nrows; i++){
      for(j = 0; j < kernel_length; j++){
	index = i-j;
	if(index >= 0 && index < nrows){
	  compressed_column[i] += image_column[index]*kernel[j];
	}
      }
    }
  }
  return compressed;
}

double* scene_generator(unsigned int nrows, unsigned int ncols){
  double* scene = malloc(nrows*ncols*sizeof(double));
  if(scene == 0){
    printf("Could not allocate memory for scene.\n");
    return 0;
  }
  memset(scene, 0, nrows*ncols*sizeof(double));
  scene[nrows*ncols/2+nrows/2] = 1;
  return scene;
}

void insert_waveform_in_scene(unsigned int waveform_length, double complex* waveform, unsigned int nrows, unsigned int ncols, double* scene, unsigned int* nnrows, unsigned int* nncols, double complex** scene_with_waveform){
	*nncols = ncols;
	*nnrows = nrows*waveform_length;

	*scene_with_waveform = malloc((*nncols)*(*nnrows)*sizeof(double complex));
	if(scene_with_waveform == 0){
	  printf("Could not allocate memory for scene with waveform.\n");
	  return;
	}

	double complex* scenew = *scene_with_waveform;
	int i;
	for(i = 0; i < (*nncols)*(*nnrows); i++){
	  scenew[i] = 0;
	}

	for(i = 0; i < waveform_length; i++){
	  scenew[(*nncols/2)*(*nnrows)+(*nnrows/2)+i] = (*nnrows)*waveform[i];
	}
}

double complex* radar_imager(unsigned int nnrows, unsigned int nncols, double complex* scene, unsigned int altitude, float beamwidth){
  unsigned int beamcrossrange = round(altitude*tan(beamwidth));

  double complex* radar_image = malloc(nnrows*nncols*sizeof(double complex));
  if(radar_image == 0){
    printf("Could not allocate memory for radar image.\n");
    return 0;
  }

  memset(radar_image, 0, nnrows*nncols*sizeof(double complex));

  int i,j,k;
  int beam_value;
  for(i = 0; i < nncols; i++){
    for(j = 0; j < 2*beamcrossrange; j++){
      beam_value = j-beamcrossrange;
      for(k = 0; k < nnrows; k++){
	if(i+beam_value <= nncols){
	  int dist = sqrt(beam_value*beam_value + k*k);
	  if(dist <= nnrows){
	    if(i+beam_value >= 1){
	      radar_image[i*nnrows+dist] += (nnrows*nnrows)*scene[(i+beam_value)*nnrows+k]/(dist*dist);
	    }
	  }
	}
      }
    }
  }
 return radar_image;
}

void gbp(unsigned int nrows, unsigned int ncols, double complex* raw_data, double complex** sar_image){
  *sar_image = malloc(nrows*ncols*sizeof(double complex));
  if(sar_image == 0){
    printf("Could not allocate memory for SAR image.\n");
    return;
  }

  double complex* image = *sar_image;

  int i;
  for(i = 0; i < nrows*ncols; i++){
    image[i] = 0;
  }

  int j,k,l;
  for(j = 0; j < ncols; j++){
    for(k = 0; k < nrows; k++){
      for(l = 0; l < ncols; l++){
	unsigned int range_index = sqrt((l-j)*(l-j)+k*k);
	if(range_index < nrows){
	  image[j*nrows+k] += raw_data[l*nrows+range_index]/nrows;
	}
      }
    }
  }

}

void gbp_fft(unsigned int nrows, unsigned int ncols, double complex* sar_data, double complex** sar_fft){
  *sar_fft = fftw_malloc(nrows*ncols*sizeof(double complex));

  if(sar_fft == 0){
    printf("Could not allocate memory for SAR FFT image.\n");
    return;
  }

  fftw_plan fft = fftw_plan_dft_2d(nrows, ncols, sar_data, *sar_fft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(fft);
  fftw_destroy_plan(fft);
}
