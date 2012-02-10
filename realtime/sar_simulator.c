#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fftw3.h>
#include <sys/time.h>

void chirp_generator(unsigned int start_frequency, unsigned int bandwidth, unsigned int* chirp_length,  double* signal_distance);
void chirp_matched_generator(unsigned int start_frequency, unsigned int bandwidth, unsigned int* chirp_length);
void scene_generator(unsigned int nrows, unsigned int ncols);
void insert_waveform_in_scene(unsigned int waveform_length, unsigned int nrows, unsigned int ncols, unsigned int* nnrows, unsigned int* nncols);
void radar_imager(unsigned int nnrows, unsigned int nncols, unsigned int altitude, float beamwidth);
void gbp(unsigned int nrows, unsigned int ncols);
void pulse_compressed(unsigned int kernel_length, unsigned int nrows, unsigned int ncols);
void gbp_fft(unsigned int nrows, unsigned int ncols);
void pulse_compressed_signal(unsigned int kernel_length);
void fft_waveform(unsigned int kernel_length, double complex* kernel, double complex* output);

#define PI 3.14159265
// 10MB
#define MEMORY_SIZE 104857600

double scene[MEMORY_SIZE/sizeof(double)];
double chirp_time_vector[MEMORY_SIZE/sizeof(double)];
double matched_time_vector[MEMORY_SIZE/sizeof(double)];
double complex chirp_fft[MEMORY_SIZE/sizeof(double)];
double complex matched_fft[MEMORY_SIZE/sizeof(double)];
double complex pulse_compressed_waveform[2*MEMORY_SIZE/sizeof(double complex)];
double complex chirp_signal[MEMORY_SIZE/sizeof(double complex)];
double complex matched_chirp[MEMORY_SIZE/sizeof(double complex)];
double complex scene_with_waveform[MEMORY_SIZE/sizeof(double complex)];
double complex radar_image[MEMORY_SIZE/sizeof(double complex)];
double complex pulse_compressed_radar_image[MEMORY_SIZE/sizeof(double complex)];
double complex sar_image[MEMORY_SIZE/sizeof(double complex)];
double complex sar_fft[MEMORY_SIZE/sizeof(double complex)];
double complex sar_img_shifted[MEMORY_SIZE/sizeof(double complex)];

int main(int argc, char** argv){
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

  memset(scene, 0, MEMORY_SIZE);
  memset(chirp_time_vector, 0, MEMORY_SIZE);
  memset(matched_time_vector, 0, MEMORY_SIZE);
  memset(chirp_signal, 0, MEMORY_SIZE);
  memset(matched_chirp, 0, MEMORY_SIZE);
  memset(scene_with_waveform, 0, MEMORY_SIZE);
  memset(radar_image, 0, MEMORY_SIZE);
  memset(pulse_compressed_radar_image, 0, MEMORY_SIZE);
  memset(sar_image, 0, MEMORY_SIZE);
  memset(sar_fft, 0, MEMORY_SIZE);
  memset(pulse_compressed_waveform, 0, MEMORY_SIZE);

  scene_generator(nrows, ncols);
  
  gettimeofday(&ntime, NULL);
  printf("Scene generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  gettimeofday(&otime, NULL);

  chirp_generator(start_frequency, bandwidth, &chirp_length, &signal_distance);

  gettimeofday(&ntime, NULL);
  printf("Chirp generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  gettimeofday(&otime, NULL);

  fft_waveform(chirp_length, chirp_signal, chirp_fft);

  gettimeofday(&ntime, NULL);
  printf("Chirp FFT generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  gettimeofday(&otime, NULL);

  chirp_matched_generator(start_frequency, bandwidth, &chirp_length);

  gettimeofday(&ntime, NULL);
  printf("Matched chirp generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  gettimeofday(&otime, NULL);
  
  fft_waveform(chirp_length, matched_chirp, matched_fft);

  gettimeofday(&ntime, NULL);
  printf("Matched chirp FFT generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  gettimeofday(&otime, NULL);

  pulse_compressed_signal(chirp_length);

  gettimeofday(&ntime, NULL);
  printf("Single pulse compression took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  gettimeofday(&otime, NULL);

  insert_waveform_in_scene(chirp_length, nrows, ncols, &nnrows, &nncols);

  gettimeofday(&ntime, NULL);
  printf("Scene with waveform generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  gettimeofday(&otime, NULL);

  radar_imager(nnrows, nncols, altitude, beamwidth);

  gettimeofday(&ntime, NULL);
  printf("Radar image generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  gettimeofday(&otime, NULL);

  pulse_compressed(chirp_length, nnrows, nncols);

  gettimeofday(&ntime, NULL);
  printf("Pulse compression of radar image took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  gettimeofday(&otime, NULL);

  gbp(nnrows, nncols);

  gettimeofday(&ntime, NULL);
  printf("GBP of pulse-compressed radar image took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  gettimeofday(&otime, NULL);

  gbp_fft(nnrows, nncols);

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
  printf("Number of complex points: %i\n",nncols*nnrows);

  FILE* dimensions = fopen("dimensions.dat", "w");
  FILE* scenef = fopen("scene.dat", "wb");
  FILE* chirpf = fopen("chirp.dat", "wb");
  FILE* matchedf = fopen("matched.dat", "wb");
  FILE* chirpfftf = fopen("chirpfft.dat", "wb");
  FILE* matchedfftf = fopen("matchedfft.dat", "wb");
  FILE* compressedf = fopen("compressed.dat", "wb");
  FILE* scene_with_waveformf = fopen("scene_with_waveform.dat", "wb");
  FILE* radar_imagef = fopen("radar_image.dat", "wb");
  FILE* pulse_compressedf = fopen("pulse_compressed_image.dat", "wb");
  FILE* sar_imagef = fopen("sar_image.dat", "wb");
  FILE* sar_fftf = fopen("sar_fft.dat", "wb");

  fprintf(dimensions, "%u\n%u\n%u\n%u\n%u\n%f\n", chirp_length, nrows, ncols, nnrows, nncols, signal_distance);
  /*
  fwrite(scene, 1, nrows*ncols*sizeof(double), scenef);
  fwrite(chirp_signal, 1, chirp_length*sizeof(complex double), chirpf);
  fwrite(matched_chirp, 1, chirp_length*sizeof(complex double), matchedf);
  fwrite(pulse_compressed_waveform, 1, chirp_length*sizeof(complex double), compressedf);
  fwrite(scene_with_waveform, 1, nnrows*nncols*sizeof(complex double), scene_with_waveformf);
  fwrite(radar_image, 1, nnrows*nncols*sizeof(complex double), radar_imagef);
  fwrite(pulse_compressed_radar_image, 1, nnrows*nncols*sizeof(complex double), pulse_compressedf);
  fwrite(sar_image, 1, nnrows*nncols*sizeof(complex double), sar_imagef);
  fwrite(sar_fft, 1, nnrows*nncols*sizeof(double complex), sar_fftf);
  */
  
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

  for(i = 0; i < chirp_length; i++){
    fprintf(compressedf, "%f\t", creal(pulse_compressed_waveform[i]));
    fprintf(compressedf, "%f\n", cimag(pulse_compressed_waveform[i]));
  }

  for(i = 0; i < chirp_length; i++){
    fprintf(chirpfftf, "%f\t", creal(chirp_fft[i]));
    fprintf(chirpfftf, "%f\n", cimag(chirp_fft[i]));
  }

  for(i = 0; i < chirp_length; i++){
    fprintf(matchedfftf, "%f\t", creal(matched_fft[i]));
    fprintf(matchedfftf, "%f\n", cimag(matched_fft[i]));
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
  fclose(chirpfftf);
  fclose(matchedfftf);
  fclose(compressedf);
  fclose(scene_with_waveformf);
  fclose(radar_imagef);
  fclose(sar_imagef);
}

void chirp_generator(unsigned int start_frequency, unsigned int bandwidth, unsigned int *chirp_length, double* signal_distance){
  /* 
   * The following relation holds for a linear chirp signal:
   * f(t) = f_0 + chirp_rate*t
   * If a signal of a bandwidth is required, we have:
   * f_0 + bandwidth = f_0 + chirp_rate*t
   * bandwidth = chirp_rate*t
   *
   * For a bandwidth-time product of 100, we have:
   * bandwidth*t = 100 => chirp_rate*t^2 = 100
   * t = sqrt(100/chirp_rate)
   */
  float btproduct = 100;
  float end_time = btproduct/bandwidth;
  float chirp_rate = bandwidth/end_time;

  unsigned int sample_frequency = 5*bandwidth;
  unsigned int time_steps = end_time*sample_frequency;
  *chirp_length = time_steps;

  *signal_distance = end_time*300000000;

  int i;
  float last_time = 0;
  for(i = 0; i < end_time*sample_frequency; i++){
    chirp_time_vector[i] = last_time;
    last_time += (float)1/sample_frequency;
  }

  for(i = 0; i < end_time*sample_frequency; i++){
    chirp_signal[i] = cos( 2*PI*(start_frequency*chirp_time_vector[i]+chirp_rate*chirp_time_vector[i]*chirp_time_vector[i] ))+_Complex_I*sin( 2*PI*(start_frequency*chirp_time_vector[i]+chirp_rate*chirp_time_vector[i]*chirp_time_vector[i] ));
  }
}

void chirp_matched_generator(unsigned int start_frequency, unsigned int bandwidth, unsigned int *chirp_length){
  float btproduct = 100;
  float end_time = btproduct/bandwidth;
  float chirp_rate = bandwidth/end_time;
  
  unsigned int sample_frequency = 5*bandwidth;
  unsigned int time_steps = end_time*sample_frequency;
  *chirp_length = time_steps;

  int i;
  float last_time = end_time;
  for(i = 0; i < end_time*sample_frequency; i++){
    matched_time_vector[i] = last_time;
    last_time -= (float)1/sample_frequency;
  }

  for(i = 0; i < end_time*sample_frequency; i++){
    matched_chirp[i] = cos( 2*PI*(start_frequency*matched_time_vector[i]+chirp_rate*matched_time_vector[i]*matched_time_vector[i] ))-_Complex_I*sin( 2*PI*(start_frequency*chirp_time_vector[i]+chirp_rate*chirp_time_vector[i]*chirp_time_vector[i] ));
  }
}

void pulse_compressed_signal(unsigned int kernel_length){
  unsigned int filter_length = 2*kernel_length;
  fftw_complex* padded_signal = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* padded_kernel = fftw_malloc(filter_length*sizeof(fftw_complex));
  memset(padded_signal, 0, filter_length*sizeof(fftw_complex));
  memset(padded_kernel, 0, filter_length*sizeof(fftw_complex));
  memcpy(padded_signal, chirp_signal, kernel_length*sizeof(fftw_complex));
  memcpy(padded_kernel, matched_chirp, kernel_length*sizeof(fftw_complex));

  fftw_complex* sigfft = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* matchfft = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* product = fftw_malloc(filter_length*sizeof(fftw_complex));

  fftw_plan sigp = fftw_plan_dft_1d(filter_length, padded_signal, sigfft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan matchp = fftw_plan_dft_1d(filter_length, padded_kernel, matchfft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan iff = fftw_plan_dft_1d(kernel_length, product, pulse_compressed_waveform, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftw_execute(sigp);
  fftw_execute(matchp);

  int i;
  for(i = 0; i < filter_length; i++)
    product[i] = sigfft[i]*matchfft[i];

  fftw_execute(iff);

  fftw_destroy_plan(sigp);
  fftw_destroy_plan(matchp);
  fftw_destroy_plan(iff);

  fftw_free(sigfft);
  fftw_free(matchfft);
  fftw_free(product);
}

void pulse_compressed(unsigned int kernel_length, unsigned int nrows, unsigned int ncols){
  // Make sure input has valid values.
  int z;
  for(z = 0; z < kernel_length; z++){
    if(isnan(matched_chirp[z]))
      matched_chirp[z] = 0;
  }

  unsigned int filter_length = nrows + kernel_length;
  fftw_complex* padded_kernel = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* kernel_fft = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* product = fftw_malloc(filter_length*sizeof(fftw_complex));
  memset(padded_kernel, 0, filter_length*sizeof(fftw_complex));
  memset(kernel_fft, 0, filter_length*sizeof(fftw_complex));
  memcpy(padded_kernel, matched_chirp, kernel_length*sizeof(fftw_complex));

  // Compute fft of filter kernel.
  fftw_plan kernelfft = fftw_plan_dft_1d(filter_length, padded_kernel, kernel_fft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(kernelfft);

  fftw_complex* column;
  fftw_complex* output_column;
  fftw_complex* padded_column = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* padded_column_fft = fftw_malloc(filter_length*sizeof(fftw_complex));
  int i,j,k;
  for(i = 0; i < ncols; i++){
    column = &radar_image[i*nrows];

    // Make sure we have valid values.
    for(k = 0; k < nrows; k++){
	if(isnan(column[k]))
	  column[k] = 0;
    }

    output_column = &pulse_compressed_radar_image[i*nrows];
    memset(padded_column, 0, filter_length*sizeof(fftw_complex));
    memcpy(padded_column, column, nrows*sizeof(fftw_complex));
    memset(padded_column_fft, 0, filter_length*sizeof(fftw_complex));
    fftw_plan colfft = fftw_plan_dft_1d(filter_length, padded_column, padded_column_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(colfft);


    for(j = 0; j < filter_length; j++){
	product[j] = padded_column_fft[j]*kernel_fft[j];
    }

    fftw_plan colifft = fftw_plan_dft_1d(nrows, product, output_column, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(colifft);

    fftw_destroy_plan(colfft);
    fftw_destroy_plan(colifft);

  }

  fftw_free(product);
  fftw_free(padded_column_fft);
  fftw_free(padded_column);
  fftw_free(padded_kernel);

  fftw_destroy_plan(kernelfft);
}

void scene_generator(unsigned int nrows, unsigned int ncols){
  scene[nrows*ncols/2+nrows/2] = 1;
}

void insert_waveform_in_scene(unsigned int waveform_length, unsigned int nrows, unsigned int ncols, unsigned int* nnrows, unsigned int* nncols){
	*nncols = ncols;
	*nnrows = nrows*waveform_length;

	int i;
	for(i = 0; i < waveform_length; i++){
	  scene_with_waveform[(*nncols/2)*(*nnrows)+(*nnrows/2)+i] = (*nnrows)*chirp_signal[i];
	}
}

void radar_imager(unsigned int nnrows, unsigned int nncols, unsigned int altitude, float beamwidth){
  unsigned int beamcrossrange = round(altitude*tan(beamwidth));

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
	      radar_image[i*nnrows+dist] += (nnrows*nnrows)*scene_with_waveform[(i+beam_value)*nnrows+k]/(dist*dist);
	    }
	  }
	}
      }
    }
  }
}

void gbp(unsigned int nrows, unsigned int ncols){
  int j,k,l;
  for(j = 0; j < ncols; j++){
    for(k = 0; k < nrows; k++){
      for(l = 0; l < ncols; l++){
	unsigned int range_index = sqrt((l-j)*(l-j)+k*k);
	if(range_index < nrows){
	  sar_image[j*nrows+k] += pulse_compressed_radar_image[l*nrows+range_index]/nrows;
	}
      }
    }
  }

}

void gbp_fft(unsigned int nrows, unsigned int ncols){
  int i, j;
  for(i = 0; i < ncols; i++){
    for(j = 0; j < nrows; j++){
	sar_img_shifted[nrows*i+j] = sar_image[nrows*i+j]*pow(-1,i+j);
    }
  }

  fftw_plan fft = fftw_plan_dft_2d(ncols, nrows, sar_img_shifted, sar_fft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(fft);
  fftw_destroy_plan(fft);
}

void fft_waveform(unsigned int kernel_length, double complex* kernel, double complex* output){
  fftw_plan fft = fftw_plan_dft_1d(kernel_length, kernel, output, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(fft);
  fftw_destroy_plan(fft);
}
