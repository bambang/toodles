/*
 * This SAR simulation and processing software is written by Alexander Rajula.
 * The software is free to use and modify.
 *
 * You may contact me at alexander@rajula.org
 */

#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fftw3.h>
#include <sys/time.h>
#include "sar_simulator.h"

double chirp_time_vector[MEMORY_SIZE/sizeof(double)];
double matched_time_vector[MEMORY_SIZE/sizeof(double)];
double complex chirp_fft[MEMORY_SIZE/sizeof(double)];
double complex matched_fft[MEMORY_SIZE/sizeof(double)];
double complex pulse_compressed_waveform[MEMORY_SIZE/sizeof(double complex)];
double complex chirp_signal[MEMORY_SIZE/sizeof(double complex)];
double complex matched_chirp[MEMORY_SIZE/sizeof(double complex)];
double complex scene_with_waveform[MEMORY_SIZE/sizeof(double complex)];
double complex radar_image[MEMORY_SIZE/sizeof(double complex)];
double complex pulse_compressed_radar_image[MEMORY_SIZE/sizeof(double complex)];
double complex sar_image[MEMORY_SIZE/sizeof(double complex)];
double complex sar_fft[MEMORY_SIZE/sizeof(double complex)];
double complex sar_img_shifted[MEMORY_SIZE/sizeof(double complex)];

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

int main(int argc, char** argv){
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

  printf("Do you wish to simulate or process radar data? (s/p): ");
  mode = getchar();
  int ret;
  if(mode == 'p'){
    printf("Please enter file name of raw data: ");
    ret = scanf("%s", radar_data_filename);
    if(ret == EOF){
	printf("Invalid input detected, closing.\n");
	return;
    }
    if(!read_radar_file()){
	printf("Failed to read radar data, closing.\n");
	return;
    }
    process_data();
  }
  else if(mode == 's'){
    printf("Simulate with real or complex values? (r/c): ");
    ret = scanf("%s", real_or_complex_simulation);
    if(*real_or_complex_simulation != 'r')
      *real_or_complex_simulation = 'c';

    printf("Antenna azimuth beamwidth in radians: ");
    ret = scanf("%f", &beamwidth);
    
    printf("Chirp start frequency: ");
    ret = scanf("%li", &start_frequency);

    if(start_frequency < 0){
	printf("Negative start frequency entered, closing.\n");
	return;
    }

    printf("Chirp bandwidth: ");
    ret = scanf("%li", &bandwidth);

    if(bandwidth < 0){
	printf("Negative bandwidth entered, closing.\n");
	return;
    }

    printf("Chirp bandwidth-time product: ");
    ret = scanf("%f", &btproduct);
    
    if(btproduct < 1){
	printf("Too small bandwidth-time product entered, closing.\n");
	return;
    }

    ret = simulate();
    if(ret == -1)
      return;
    process_data();
  }
  else{
    printf("Mode not recognized - exiting.\n");
    return;
  }

  printf("BT-product: %f\n", btproduct);
  printf("Chirp length: %u\n", chirp_length);
  printf("Number of rows in final scene: %i\n", nrows);
  printf("Number of columns in final scene: %i\n", ncols);
  printf("Platform altitude: %i\n", altitude);
  printf("Antenna beamwidth: %f\n", beamwidth);
  printf("Start frequency: %lu\n", start_frequency);
  printf("Bandwidth: %lu\n", bandwidth);
  printf("Number of complex points: %i\n",ncols*nrows);
  printf("Signal distance: %fm\n", signal_distance);

  if(*real_or_complex_simulation == 'c')
    ret = write_complex_data();
  else
    ret = write_real_data();
}

int read_radar_file(){
  FILE* fp = fopen(radar_data_filename, "r");
  if(fp == NULL)
    return -1;

  FILE* mp = fopen("radar_metadata", "r");
  int ret;
  radar_metadata meta;
  ret = fread(&meta, sizeof(radar_metadata), 1, mp);

  fclose(mp);

  int i;
  if(meta.real_or_complex == 'r'){
    for(i = 0; i < meta.rows*meta.cols; i++){
      ret = fread(radar_image + i*sizeof(complex double), sizeof(double), 1, fp);
    }
  }
  else if(meta.real_or_complex == 'c'){
    for(i = 0; i < meta.rows*meta.cols; i++){
      ret = fread(radar_image + i*sizeof(complex double), sizeof(complex double), 1, fp);
    }
  }
  else{
    printf("Invalid data mode, should be real or complex, read %c.\n", meta.real_or_complex);
    fclose(fp);
    return -1;
  }

  fclose(fp);
}

int simulate(){
  struct timeval otime, ntime;

  gettimeofday(&otime, NULL);

  chirp_generator();

  //gettimeofday(&ntime, NULL);
  //printf("Chirp generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  //gettimeofday(&otime, NULL);

  fft_waveform(chirp_length, chirp_signal, chirp_fft);

  //gettimeofday(&ntime, NULL);

  //gettimeofday(&ntime, NULL);
  //printf("Chirp FFT generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  //gettimeofday(&otime, NULL);

  chirp_matched_generator();

  //gettimeofday(&ntime, NULL);
  //printf("Matched chirp generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  //gettimeofday(&otime, NULL);
 
  fft_waveform(chirp_length, matched_chirp, matched_fft);

  //gettimeofday(&ntime, NULL);
  //printf("Matched chirp FFT generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  //gettimeofday(&otime, NULL);

  pulse_compress_signal();

  //gettimeofday(&ntime, NULL);
  //printf("Single pulse compression took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  //gettimeofday(&otime, NULL);

  printf("One chirp pulse covers %f meters.\n", signal_distance);
  printf("The target will be placed in the middle of the simulated area.\n");
  printf("Enter area azimuth length (m): ");
  float len = 0;
  int ret = 0;
  ret = scanf("%f", &len);
  ncols = len*chirp_length/signal_distance;
  if(ncols < 2){
    printf("Invalid azimuth length, exiting.\n");
    return -1;
  }
  printf("Enter area range (m): ");
  ret = scanf("%f", &len);
  nrows = len*chirp_length/signal_distance;
  if(nrows < chirp_length){
    printf("Too small range, exiting.\n");
    return -1;
  }

  insert_waveform_in_scene();

  //gettimeofday(&ntime, NULL);
  //printf("Scene with waveform generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  //gettimeofday(&otime, NULL);

  printf("Scene range: %fm\n", signal_distance*(nrows/chirp_length));
  printf("Scene azimuth length: %fm\n", signal_distance*(ncols/chirp_length));
  printf("Please input SAR platform height: ");
  ret = scanf("%d", &altitude);
  if(ret == EOF){
    printf("Invalid input detected, closing.\n");
    return;
  }

  radar_imager();

  //gettimeofday(&ntime, NULL);
  //printf("Radar image generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  //gettimeofday(&otime, NULL);
}

void process_data(){
  struct timeval otime, ntime;

  //gettimeofday(&otime, NULL);

  printf("Do you want to enable pulse compression (y/n)? ");
  char pc = 0;
  int ret;
  do{
    ret = scanf("%c", &pc);
    if(pc == 'y')
      break;
    else if(pc == 'n')
      break;
  }while(1);
  if(pc == 'y')
    pulse_compress_image();
  else if(pc == 'n'){
    memcpy(pulse_compressed_radar_image, radar_image, nrows*ncols*sizeof(complex double));
  }

  //gettimeofday(&ntime, NULL);
  //printf("Pulse compression of radar image took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  //gettimeofday(&otime, NULL);

  gbp();

  //gettimeofday(&ntime, NULL);
  //printf("GBP of pulse-compressed radar image took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  //gettimeofday(&otime, NULL);

  gbp_fft();

  //gettimeofday(&ntime, NULL);
  //printf("FFT generation of GBP image took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));

  if(mode == 's')
    filter_dc();
}

int write_complex_data(){
  char fmode = 0;
  int ret = 0;

  do{
    fmode = getchar();
  }while(fmode != '\n');

  printf("Would you like to write data in human-readable or binary format (h/b): ");
  do{
    ret = scanf("%c", &fmode);
  }while((fmode != 'h') && (fmode != 'b'));

  FILE* dimensions = fopen("dimensions.dat", "w");
  if(dimensions == NULL){
    printf("Could not open dimensions.dat for writing - exiting.\n");
    return -1;
  }
  FILE* chirpf = fopen("chirp.dat", "wb");
  if(chirpf == NULL){
    printf("Could not open chirp.dat for writing - exiting.\n");
    return -1;
  }
  FILE* matchedf = fopen("matched.dat", "wb");
  if(matchedf == NULL){
    printf("Could not open matched.dat for writing - exiting.\n");
    return -1;
  }
  FILE* chirpfftf = fopen("chirpfft.dat", "wb");
  if(chirpfftf == NULL){
    printf("Could not open chirpfft.dat for writing - exiting.\n");
    return -1;
  }
  FILE* matchedfftf = fopen("matchedfft.dat", "wb");
  if(matchedfftf == NULL){
    printf("Could not open matchedfft.dat for writing - exiting.\n");
    return -1;
  }
  FILE* compressedf = fopen("compressed.dat", "wb");
  if(compressedf == NULL){
    printf("Could not open compressed.dat for writing - exiting.\n");
    return -1;
  }
  FILE* scene_with_waveformf = fopen("scene_with_waveform.dat", "wb");
  if(scene_with_waveformf == NULL){
    printf("Could not open scene_with_waveform.dat for writing - exiting.\n");
    return -1;
  }
  FILE* radar_imagef = fopen("radar_image.dat", "wb");
  if(radar_imagef == NULL){
    printf("Could not open radar_image.dat for writing - exiting.\n");
    return -1;
  }
  FILE* pulse_compressedf = fopen("pulse_compressed_image.dat", "wb");
  if(pulse_compressedf == NULL){
   printf("Could not open pulse_compressed_image.dat for writing - exiting.\n");
   return -1;
  }
  FILE* sar_imagef = fopen("sar_image.dat", "wb");
  if(sar_imagef == NULL){
    printf("Could not open sar_image.dat for writing - exiting.\n");
    return -1;
  }
  FILE* sar_fftf = fopen("sar_fft.dat", "wb");
  if(sar_fftf == NULL){
    printf("Could not open sar_fft.dat for writing - exiting.\n");
    return -1;
  }

  fprintf(dimensions, "%u\n%u\n%u\n%f\n", chirp_length, nrows, ncols, signal_distance);
  
  if(fmode == 'b'){
	    ret = fwrite(chirp_signal, 1, chirp_length*sizeof(complex double), chirpf);
	    ret = fwrite(matched_chirp, 1, chirp_length*sizeof(complex double), matchedf);
	    ret = fwrite(chirp_fft, 1, chirp_length*sizeof(complex double), chirpfftf);
	    ret = fwrite(matched_fft, 1, chirp_length*sizeof(complex double), matchedfftf);
	    ret = fwrite(pulse_compressed_waveform, 1, chirp_length*sizeof(complex double), compressedf);
	    ret = fwrite(scene_with_waveform, 1, nrows*ncols*sizeof(complex double), scene_with_waveformf);
	    ret = fwrite(radar_image, 1, nrows*ncols*sizeof(complex double), radar_imagef);
	    ret = fwrite(pulse_compressed_radar_image, 1, nrows*ncols*sizeof(complex double), pulse_compressedf);
	    ret = fwrite(sar_image, 1, nrows*ncols*sizeof(complex double), sar_imagef);
	    ret = fwrite(sar_fft, 1, nrows*ncols*sizeof(double complex), sar_fftf);
  }
  else{
	  int i,j;
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

	  for(i = 0; i < ncols; i++){
	    for(j = 0; j < nrows; j++){
	      fprintf(scene_with_waveformf, "%g\t", creal(scene_with_waveform[i*nrows+j]));
	      fprintf(scene_with_waveformf, "%g\t", cimag(scene_with_waveform[i*nrows+j]));
	    }
	    fprintf(scene_with_waveformf, "\n");
	  }

	  for(i = 0; i < ncols; i++){
	    for(j = 0; j < nrows; j++){
	      fprintf(radar_imagef, "%f\t", creal(radar_image[i*nrows+j]));
	      fprintf(radar_imagef, "%f\t", cimag(radar_image[i*nrows+j]));
	    }
	    fprintf(radar_imagef, "\n");
	  }

	  for(i = 0; i < ncols; i++){
	    for(j = 0; j < nrows; j++){
	      fprintf(pulse_compressedf, "%f\t", creal(pulse_compressed_radar_image[i*nrows+j]));
	      fprintf(pulse_compressedf, "%f\t", cimag(pulse_compressed_radar_image[i*nrows+j]));
	    }
	    fprintf(pulse_compressedf, "\n");
	  }

	  for(i = 0; i < ncols; i++){
	    for(j = 0; j < nrows; j++){
	      fprintf(sar_imagef, "%f\t", creal(sar_image[i*nrows+j]));
	      fprintf(sar_imagef, "%f\t", cimag(sar_image[i*nrows+j]));
	    }
	    fprintf(sar_imagef, "\n");
	  }

	  for(i = 0; i < ncols; i++){
	    for(j = 0; j < nrows; j++){
	      fprintf(sar_fftf, "%f\t", creal(sar_fft[i*nrows+j]));
	      fprintf(sar_fftf, "%f\t", cimag(sar_fft[i*nrows+j]));
	    }
	    fprintf(sar_fftf, "\n");
  	  }
  }

  fclose(sar_fftf);
  fclose(pulse_compressedf);
  fclose(dimensions);
  fclose(chirpf);
  fclose(matchedf);
  fclose(chirpfftf);
  fclose(matchedfftf);
  fclose(compressedf);
  fclose(scene_with_waveformf);
  fclose(radar_imagef);
  fclose(sar_imagef);
}

int write_real_data(){
  char fmode = 0;
  int ret = 0;

  do{
    fmode = getchar();
  }while(fmode != '\n');

  printf("Would you like to write data in human-readable or binary format (h/b): ");
  do{
    ret = scanf("%c", &fmode);
  }while((fmode != 'h') && (fmode != 'b'));

  FILE* dimensions = fopen("dimensions.dat", "w");
  if(dimensions == NULL){
    printf("Could not open dimensions.dat for writing - exiting.\n");
    return -1;
  }
  FILE* chirpf = fopen("chirp.dat", "wb");
  if(chirpf == NULL){
    printf("Could not open chirp.dat for writing - exiting.\n");
    return -1;
  }
  FILE* matchedf = fopen("matched.dat", "wb");
  if(matchedf == NULL){
    printf("Could not open matched.dat for writing - exiting.\n");
    return -1;
  }
  FILE* chirpfftf = fopen("chirpfft.dat", "wb");
  if(chirpfftf == NULL){
    printf("Could not open chirpfft.dat for writing - exiting.\n");
    return -1;
  }
  FILE* matchedfftf = fopen("matchedfft.dat", "wb");
  if(matchedfftf == NULL){
    printf("Could not open matchedfft.dat for writing - exiting.\n");
    return -1;
  }
  FILE* compressedf = fopen("compressed.dat", "wb");
  if(compressedf == NULL){
    printf("Could not open compressed.dat for writing - exiting.\n");
    return -1;
  }
  FILE* scene_with_waveformf = fopen("scene_with_waveform.dat", "wb");
  if(scene_with_waveformf == NULL){
    printf("Could not open scene_with_waveform.dat for writing - exiting.\n");
    return -1;
  }
  FILE* radar_imagef = fopen("radar_image.dat", "wb");
  if(radar_imagef == NULL){
    printf("Could not open radar_image.dat for writing - exiting.\n");
    return -1;
  }
  FILE* pulse_compressedf = fopen("pulse_compressed_image.dat", "wb");
  if(pulse_compressedf == NULL){
   printf("Could not open pulse_compressed_image.dat for writing - exiting.\n");
   return -1;
  }
  FILE* sar_imagef = fopen("sar_image.dat", "wb");
  if(sar_imagef == NULL){
    printf("Could not open sar_image.dat for writing - exiting.\n");
    return -1;
  }
  FILE* sar_fftf = fopen("sar_fft.dat", "wb");
  if(sar_fftf == NULL){
    printf("Could not open sar_fft.dat for writing - exiting.\n");
    return -1;
  }

  fprintf(dimensions, "%u\n%u\n%u\n%f\n", chirp_length, nrows, ncols, signal_distance);
  
  if(fmode == 'b'){
	    ret = fwrite(chirp_signal, 1, chirp_length*sizeof(complex double), chirpf);
	    ret = fwrite(matched_chirp, 1, chirp_length*sizeof(complex double), matchedf);
	    ret = fwrite(chirp_fft, 1, chirp_length*sizeof(complex double), chirpfftf);
	    ret = fwrite(matched_fft, 1, chirp_length*sizeof(complex double), matchedfftf);
	    ret = fwrite(pulse_compressed_waveform, 1, chirp_length*sizeof(complex double), compressedf);
	    ret = fwrite(scene_with_waveform, 1, nrows*ncols*sizeof(complex double), scene_with_waveformf);
	    ret = fwrite(radar_image, 1, nrows*ncols*sizeof(complex double), radar_imagef);
	    ret = fwrite(pulse_compressed_radar_image, 1, nrows*ncols*sizeof(complex double), pulse_compressedf);
	    ret = fwrite(sar_image, 1, nrows*ncols*sizeof(complex double), sar_imagef);
	    ret = fwrite(sar_fft, 1, nrows*ncols*sizeof(double complex), sar_fftf);
  }
  else{
	  int i,j;
	  for(i = 0; i < chirp_length; i++){
	    fprintf(chirpf, "%f\n", creal(chirp_signal[i]));
	  }

	  for(i = 0; i < chirp_length; i++){
	    fprintf(matchedf, "%f\n", creal(matched_chirp[i]));
	  }

	  for(i = 0; i < chirp_length; i++){
	    fprintf(compressedf, "%f\n", creal(pulse_compressed_waveform[i]));
	  }

	  for(i = 0; i < chirp_length; i++){
	    fprintf(chirpfftf, "%f\n", creal(chirp_fft[i]));
	  }

	  for(i = 0; i < chirp_length; i++){
	    fprintf(matchedfftf, "%f\n", creal(matched_fft[i]));
	  }

	  for(i = 0; i < ncols; i++){
	    for(j = 0; j < nrows; j++){
	      fprintf(scene_with_waveformf, "%g\n", creal(scene_with_waveform[i*nrows+j]));
	    }
	    fprintf(scene_with_waveformf, "\n");
	  }

	  for(i = 0; i < ncols; i++){
	    for(j = 0; j < nrows; j++){
	      fprintf(radar_imagef, "%f\n", creal(radar_image[i*nrows+j]));
	    }
	    fprintf(radar_imagef, "\n");
	  }

	  for(i = 0; i < ncols; i++){
	    for(j = 0; j < nrows; j++){
	      fprintf(pulse_compressedf, "%f\n", creal(pulse_compressed_radar_image[i*nrows+j]));
	    }
	    fprintf(pulse_compressedf, "\n");
	  }

	  for(i = 0; i < ncols; i++){
	    for(j = 0; j < nrows; j++){
	      fprintf(sar_imagef, "%f\n", creal(sar_image[i*nrows+j]));
	    }
	    fprintf(sar_imagef, "\n");
	  }

	  for(i = 0; i < ncols; i++){
	    for(j = 0; j < nrows; j++){
	      fprintf(sar_fftf, "%f\n", creal(sar_fft[i*nrows+j]));
	    }
	    fprintf(sar_fftf, "\n");
  	  }
  }

  fclose(sar_fftf);
  fclose(pulse_compressedf);
  fclose(dimensions);
  fclose(chirpf);
  fclose(matchedf);
  fclose(chirpfftf);
  fclose(matchedfftf);
  fclose(compressedf);
  fclose(scene_with_waveformf);
  fclose(radar_imagef);
}

void chirp_generator(){
  double end_time = btproduct/bandwidth;
  double chirp_rate = bandwidth/end_time;

  unsigned long int sample_frequency = 5*bandwidth;
  unsigned long int time_steps = end_time*sample_frequency;
  chirp_length = time_steps;

  signal_distance = end_time*C;

  int i;
  double last_time = (float)10/sample_frequency;
  for(i = 0; i < end_time*sample_frequency; i++){
    chirp_time_vector[i] = last_time;
    last_time += (float)1/sample_frequency;
  }

  for(i = 0; i < time_steps; i++){
    chirp_signal[i] = cos( 2*PI*(start_frequency*chirp_time_vector[i]+chirp_rate*chirp_time_vector[i]*chirp_time_vector[i] ))+_Complex_I*sin( 2*PI*(start_frequency*chirp_time_vector[i]+chirp_rate*chirp_time_vector[i]*chirp_time_vector[i] ));
  }
}

void chirp_matched_generator(){
  float end_time = btproduct/bandwidth;
  float chirp_rate = bandwidth/end_time;
  
  unsigned long int sample_frequency = 5*bandwidth;
  unsigned long int time_steps = end_time*sample_frequency;
  chirp_length = time_steps;

  int i;
  double last_time = end_time+(float)10/sample_frequency;
  for(i = 0; i < end_time*sample_frequency; i++){
    matched_time_vector[i] = last_time;
    last_time -= (float)1/sample_frequency;
  }

  for(i = 0; i < end_time*sample_frequency; i++){
    matched_chirp[i] = cos( 2*PI*(start_frequency*matched_time_vector[i]+chirp_rate*matched_time_vector[i]*matched_time_vector[i] ))-_Complex_I*sin( 2*PI*(start_frequency*chirp_time_vector[i]+chirp_rate*chirp_time_vector[i]*chirp_time_vector[i] ));
  }
}

void pulse_compress_signal(){
  int kernel_length = chirp_length;
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

void pulse_compress_image(){
  // Make sure input has valid values.
  int kernel_length = chirp_length;
  int z;
  for(z = 0; z < kernel_length; z++){
    if(isnan(matched_chirp[z]))
      matched_chirp[z] = 0;
  }

  unsigned int filter_length = nrows + kernel_length;
  filter_length = pow(2, ceil(log(filter_length)/log(2)));

  fftw_complex* padded_kernel = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* kernel_fft = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* product = fftw_malloc(filter_length*sizeof(fftw_complex));
  memset(padded_kernel, 0, filter_length*sizeof(fftw_complex));
  memset(kernel_fft, 0, filter_length*sizeof(fftw_complex));
  memcpy(padded_kernel, matched_chirp, kernel_length*sizeof(fftw_complex));

  // Compute fft of filter kernel.
  fftw_plan kernelfft = fftw_plan_dft_1d(filter_length, padded_kernel, kernel_fft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(kernelfft);
  fftw_complex* output_column;
  fftw_complex* padded_column = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* padded_column_fft = fftw_malloc(filter_length*sizeof(fftw_complex));
  int i,j,k;
  for(i = 0; i < ncols; i++){
    double complex* column = &radar_image[i*nrows];

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

    // Re-normalize ifft output.
    for(j = 0; j < nrows; j++)
      output_column[j] /= nrows;
  }

  fftw_free(product);
  fftw_free(padded_column_fft);
  fftw_free(padded_column);
  fftw_free(padded_kernel);

  fftw_destroy_plan(kernelfft);
}

void insert_waveform_in_scene(){
	int i;
	for(i = 0; i < chirp_length; i++){
	  scene_with_waveform[ (int)(ncols/2)*nrows + (int)(nrows/2) -chirp_length/2 +i] = chirp_signal[i];
	}
}

void radar_imager(){
  unsigned int beamcrossrange = round(altitude*tan(beamwidth));

  unsigned int i,j,k;
  int beam_value;
  for(i = 0; i < ncols; i++){
    for(j = 0; j < 2*beamcrossrange; j++){
      beam_value = j-beamcrossrange;
      for(k = 0; k < nrows; k++){
	if(i+beam_value <= ncols){
	  int dist = sqrt(pow(beam_value,2)+pow(k, 2));
	  if(dist <= nrows){
	    if(i+beam_value >= 1){
	      radar_image[i*nrows+dist] += scene_with_waveform[(i+beam_value)*nrows+k]/pow(dist,2);
	    }
	  }
	}
      }
    }
  }

}

void gbp(){
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

void gbp_fft(){
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

void filter_dc(){
  int i,j;
  for(i = ncols/2 - 5; i < ncols/2+5; i++){
    for(j = nrows/2 - 5; j < nrows/2 + 5; j++){
	sar_fft[i*nrows+j] = 0;
    }
  }

  fftw_plan ifft = fftw_plan_dft_2d(ncols, nrows, sar_fft, sar_image, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(ifft);
  fftw_destroy_plan(ifft);
}
