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

radar_variables variables;
data_arrays data;

int main(int argc, char** argv){
  printf("Do you wish to simulate or process radar data? (s/p): ");
  variables.mode = getchar();
  int ret;
  if(variables.mode == 'p'){
    printf("Please enter file name of raw data: ");
    ret = scanf("%s", variables.radar_data_filename);
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
  else if(variables.mode == 's'){
    printf("Simulate with real or complex values? (r/c): ");
    ret = scanf("%s", variables.real_or_complex_simulation);
    if(*variables.real_or_complex_simulation != 'r')
      *variables.real_or_complex_simulation = 'c';

    printf("Antenna azimuth beamwidth in radians: ");
    ret = scanf("%f", &variables.beamwidth);
    
    printf("Chirp start frequency: ");
    ret = scanf("%li", &variables.start_frequency);

    if(variables.start_frequency < 0){
	printf("Negative start frequency entered, closing.\n");
	return;
    }

    printf("Chirp bandwidth: ");
    ret = scanf("%li", &variables.bandwidth);

    if(variables.bandwidth < 0){
	printf("Negative bandwidth entered, closing.\n");
	return;
    }

    printf("Chirp bandwidth-time product: ");
    ret = scanf("%f", &variables.btproduct);
    
    if(variables.btproduct < 1){
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

  printf("Chirp length: %u\n", variables.chirp_length);
  printf("Number of rows in final scene: %i\n", variables.nrows);
  printf("Number of columns in final scene: %i\n", variables.ncols);
  printf("Number of complex points: %i\n",variables.ncols*variables.nrows);
  printf("Uncompressed pulse resolution: %fm\n", variables.signal_distance);
  printf("Compressed pulse resolution: %lfm\n", calculate_compressed_pulse_resolution());

  if(*variables.real_or_complex_simulation == 'c')
    ret = write_complex_data();
  else
    ret = write_real_data();

}

void free_memory(){
  if(data.chirp_time_vector)
    free(data.chirp_time_vector);
  if(data.matched_time_vector)
    free(data.matched_time_vector);
  if(data.chirp_fft)
    free(data.chirp_fft);
  if(data.matched_fft);
    free(data.matched_fft);
  if(data.pulse_compressed_waveform)
    free(data.pulse_compressed_waveform);
  if(data.chirp_signal)
    free(data.chirp_signal);
  if(data.matched_chirp)
    free(data.matched_chirp);
  if(data.scene_with_waveform)
    free(data.scene_with_waveform);
  if(data.radar_image)
    free(data.radar_image);
  if(data.pulse_compressed_radar_image)
    free(data.pulse_compressed_radar_image);
  if(data.sar_image)
    free(data.sar_image);
  if(data.sar_fft)
    free(data.sar_fft);
  if(data.sar_img_shifted)
    free(data.sar_img_shifted);
}

float calculate_compressed_pulse_resolution(){
  double complex* waveform = data.pulse_compressed_waveform;

  // Find maximum amplitude position
  int i, max_position;
  double max_amplitude = 0;
  for(i = 0; i < variables.chirp_length; i++){
    if( sqrt( pow(creal(waveform[i]), 2) + pow(cimag(waveform[i]),2) ) > max_amplitude){
      max_amplitude = waveform[i];
      max_position = i;
    }
  }
  double half_amplitude = max_amplitude / 2;
  int low_index = 0;
  double low_value;
  for(i = max_position; i >= 0; i--){
    low_value = sqrt( pow(creal(waveform[i]), 2) + pow(cimag(waveform[i]),2));
    if( low_value <= half_amplitude){
	low_index = i;
	break;
    }
  }
  int high_index = 0;
  for(i = max_position; i < variables.chirp_length; i++){
    low_value = sqrt( pow(creal(waveform[i]), 2) + pow(cimag(waveform[i]),2));
    if( low_value <= half_amplitude){
	high_index = i;
	break;
    }
  }

  return variables.signal_distance*(high_index-low_index)/variables.chirp_length;
}

int read_radar_file(){
  FILE* fp = fopen(variables.radar_data_filename, "r");
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
      ret = fread(data.radar_image + i*sizeof(complex double), sizeof(double), 1, fp);
    }
  }
  else if(meta.real_or_complex == 'c'){
    for(i = 0; i < meta.rows*meta.cols; i++){
      ret = fread(data.radar_image + i*sizeof(complex double), sizeof(complex double), 1, fp);
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

  data.chirp_fft = malloc(variables.chirp_length*sizeof(double complex));
  fft_waveform(variables.chirp_length, data.chirp_signal, data.chirp_fft);

  //gettimeofday(&ntime, NULL);

  //gettimeofday(&ntime, NULL);
  //printf("Chirp FFT generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  //gettimeofday(&otime, NULL);

  chirp_matched_generator();

  //gettimeofday(&ntime, NULL);
  //printf("Matched chirp generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  //gettimeofday(&otime, NULL);
 
  data.matched_fft = malloc(variables.chirp_length*sizeof(double complex));
  fft_waveform(variables.chirp_length, data.matched_chirp, data.matched_fft);

  //gettimeofday(&ntime, NULL);
  //printf("Matched chirp FFT generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  //gettimeofday(&otime, NULL);

  pulse_compress_signal();

  //gettimeofday(&ntime, NULL);
  //printf("Single pulse compression took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  //gettimeofday(&otime, NULL);

  printf("The target will be placed in the middle of the simulated area.\n");
  printf("Enter area azimuth length (m): ");
  float len = 0;
  int ret = 0;
  ret = scanf("%f", &len);
  variables.ncols = len*variables.chirp_length/variables.signal_distance;
  if(variables.ncols < 2){
    printf("Invalid azimuth length, exiting.\n");
    return -1;
  }
  printf("Enter area range (m): ");
  ret = scanf("%f", &len);
  variables.nrows = len*variables.chirp_length/variables.signal_distance;
  if(variables.nrows < variables.chirp_length){
    printf("Too small range, exiting.\n");
    return -1;
  }

  insert_waveform_in_scene();

  //gettimeofday(&ntime, NULL);
  //printf("Scene with waveform generation took %lis %lfus.\n", ntime.tv_sec - otime.tv_sec, fabs(ntime.tv_usec - otime.tv_usec));
  //gettimeofday(&otime, NULL);

  printf("Scene range: %fm\n", variables.signal_distance*(variables.nrows/variables.chirp_length));
  printf("Scene azimuth length: %fm\n", variables.signal_distance*(variables.ncols/variables.chirp_length));
  printf("Please input SAR platform height: ");
  ret = scanf("%d", &variables.altitude);
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
    memcpy(data.pulse_compressed_radar_image, data.radar_image, variables.nrows*variables.ncols*sizeof(complex double));
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

  if(variables.mode == 's')
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

  fprintf(dimensions, "%u\n%u\n%u\n%f\n", variables.chirp_length, variables.nrows, variables.ncols, variables.signal_distance);
  
  if(fmode == 'b'){
	    ret = fwrite(data.chirp_signal, 1, variables.chirp_length*sizeof(complex double), chirpf);
	    ret = fwrite(data.matched_chirp, 1, variables.chirp_length*sizeof(complex double), matchedf);
	    ret = fwrite(data.chirp_fft, 1, variables.chirp_length*sizeof(complex double), chirpfftf);
	    ret = fwrite(data.matched_fft, 1, variables.chirp_length*sizeof(complex double), matchedfftf);
	    ret = fwrite(data.pulse_compressed_waveform, 1, variables.chirp_length*sizeof(complex double), compressedf);
	    ret = fwrite(data.scene_with_waveform, 1, variables.nrows*variables.ncols*sizeof(complex double), scene_with_waveformf);
	    ret = fwrite(data.radar_image, 1, variables.nrows*variables.ncols*sizeof(complex double), radar_imagef);
	    ret = fwrite(data.pulse_compressed_radar_image, 1, variables.nrows*variables.ncols*sizeof(complex double), pulse_compressedf);
	    ret = fwrite(data.sar_image, 1, variables.nrows*variables.ncols*sizeof(complex double), sar_imagef);
	    ret = fwrite(data.sar_fft, 1, variables.nrows*variables.ncols*sizeof(double complex), sar_fftf);
  }
  else{
	  int i,j;
	  for(i = 0; i < variables.chirp_length; i++){
	    fprintf(chirpf, "%f\t", creal(data.chirp_signal[i]));
	    fprintf(chirpf, "%f\n", cimag(data.chirp_signal[i]));
	  }

	  for(i = 0; i < variables.chirp_length; i++){
	    fprintf(matchedf, "%f\t", creal(data.matched_chirp[i]));
	    fprintf(matchedf, "%f\n", cimag(data.matched_chirp[i]));
	  }

	  for(i = 0; i < variables.chirp_length; i++){
	    fprintf(compressedf, "%f\t", creal(data.pulse_compressed_waveform[i]));
	    fprintf(compressedf, "%f\n", cimag(data.pulse_compressed_waveform[i]));
	  }

	  for(i = 0; i < variables.chirp_length; i++){
	    fprintf(chirpfftf, "%f\t", creal(data.chirp_fft[i]));
	    fprintf(chirpfftf, "%f\n", cimag(data.chirp_fft[i]));
	  }

	  for(i = 0; i < variables.chirp_length; i++){
	    fprintf(matchedfftf, "%f\t", creal(data.matched_fft[i]));
	    fprintf(matchedfftf, "%f\n", cimag(data.matched_fft[i]));
	  }

	  for(i = 0; i < variables.ncols; i++){
	    for(j = 0; j < variables.nrows; j++){
	      fprintf(scene_with_waveformf, "%g\t", creal(data.scene_with_waveform[i*variables.nrows+j]));
	      fprintf(scene_with_waveformf, "%g\t", cimag(data.scene_with_waveform[i*variables.nrows+j]));
	    }
	    fprintf(scene_with_waveformf, "\n");
	  }

	  for(i = 0; i < variables.ncols; i++){
	    for(j = 0; j < variables.nrows; j++){
	      fprintf(radar_imagef, "%f\t", creal(data.radar_image[i*variables.nrows+j]));
	      fprintf(radar_imagef, "%f\t", cimag(data.radar_image[i*variables.nrows+j]));
	    }
	    fprintf(radar_imagef, "\n");
	  }

	  for(i = 0; i < variables.ncols; i++){
	    for(j = 0; j < variables.nrows; j++){
	      fprintf(pulse_compressedf, "%f\t", creal(data.pulse_compressed_radar_image[i*variables.nrows+j]));
	      fprintf(pulse_compressedf, "%f\t", cimag(data.pulse_compressed_radar_image[i*variables.nrows+j]));
	    }
	    fprintf(pulse_compressedf, "\n");
	  }

	  for(i = 0; i < variables.ncols; i++){
	    for(j = 0; j < variables.nrows; j++){
	      fprintf(sar_imagef, "%f\t", creal(data.sar_image[i*variables.nrows+j]));
	      fprintf(sar_imagef, "%f\t", cimag(data.sar_image[i*variables.nrows+j]));
	    }
	    fprintf(sar_imagef, "\n");
	  }

	  for(i = 0; i < variables.ncols; i++){
	    for(j = 0; j < variables.nrows; j++){
	      fprintf(sar_fftf, "%f\t", creal(data.sar_fft[i*variables.nrows+j]));
	      fprintf(sar_fftf, "%f\t", cimag(data.sar_fft[i*variables.nrows+j]));
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

  fprintf(dimensions, "%u\n%u\n%u\n%f\n", variables.chirp_length, variables.nrows, variables.ncols, variables.signal_distance);
  
  if(fmode == 'b'){
	    ret = fwrite(data.chirp_signal, 1, variables.chirp_length*sizeof(complex double), chirpf);
	    ret = fwrite(data.matched_chirp, 1, variables.chirp_length*sizeof(complex double), matchedf);
	    ret = fwrite(data.chirp_fft, 1, variables.chirp_length*sizeof(complex double), chirpfftf);
	    ret = fwrite(data.matched_fft, 1, variables.chirp_length*sizeof(complex double), matchedfftf);
	    ret = fwrite(data.pulse_compressed_waveform, 1, variables.chirp_length*sizeof(complex double), compressedf);
	    ret = fwrite(data.scene_with_waveform, 1, variables.nrows*variables.ncols*sizeof(complex double), scene_with_waveformf);
	    ret = fwrite(data.radar_image, 1, variables.nrows*variables.ncols*sizeof(complex double), radar_imagef);
	    ret = fwrite(data.pulse_compressed_radar_image, 1, variables.nrows*variables.ncols*sizeof(complex double), pulse_compressedf);
	    ret = fwrite(data.sar_image, 1, variables.nrows*variables.ncols*sizeof(complex double), sar_imagef);
	    ret = fwrite(data.sar_fft, 1, variables.nrows*variables.ncols*sizeof(double complex), sar_fftf);
  }
  else{
	  int i,j;
	  for(i = 0; i < variables.chirp_length; i++){
	    fprintf(chirpf, "%f\n", creal(data.chirp_signal[i]));
	  }

	  for(i = 0; i < variables.chirp_length; i++){
	    fprintf(matchedf, "%f\n", creal(data.matched_chirp[i]));
	  }

	  for(i = 0; i < variables.chirp_length; i++){
	    fprintf(compressedf, "%f\n", creal(data.pulse_compressed_waveform[i]));
	  }

	  for(i = 0; i < variables.chirp_length; i++){
	    fprintf(chirpfftf, "%f\n", creal(data.chirp_fft[i]));
	  }

	  for(i = 0; i < variables.chirp_length; i++){
	    fprintf(matchedfftf, "%f\n", creal(data.matched_fft[i]));
	  }

	  for(i = 0; i < variables.ncols; i++){
	    for(j = 0; j < variables.nrows; j++){
	      fprintf(scene_with_waveformf, "%g\n", creal(data.scene_with_waveform[i*variables.nrows+j]));
	    }
	    fprintf(scene_with_waveformf, "\n");
	  }

	  for(i = 0; i < variables.ncols; i++){
	    for(j = 0; j < variables.nrows; j++){
	      fprintf(radar_imagef, "%f\n", creal(data.radar_image[i*variables.nrows+j]));
	    }
	    fprintf(radar_imagef, "\n");
	  }

	  for(i = 0; i < variables.ncols; i++){
	    for(j = 0; j < variables.nrows; j++){
	      fprintf(pulse_compressedf, "%f\n", creal(data.pulse_compressed_radar_image[i*variables.nrows+j]));
	    }
	    fprintf(pulse_compressedf, "\n");
	  }

	  for(i = 0; i < variables.ncols; i++){
	    for(j = 0; j < variables.nrows; j++){
	      fprintf(sar_imagef, "%f\n", creal(data.sar_image[i*variables.nrows+j]));
	    }
	    fprintf(sar_imagef, "\n");
	  }

	  for(i = 0; i < variables.ncols; i++){
	    for(j = 0; j < variables.nrows; j++){
	      fprintf(sar_fftf, "%f\n", creal(data.sar_fft[i*variables.nrows+j]));
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
  double end_time = variables.btproduct/variables.bandwidth;
  double chirp_rate = variables.bandwidth/end_time;

  unsigned long int sample_frequency = 5*variables.bandwidth;
  unsigned long int time_steps = end_time*sample_frequency;
  variables.chirp_length = time_steps;
  variables.signal_distance = end_time*C;

  data.chirp_time_vector = malloc(time_steps*sizeof(double));
  data.chirp_signal = malloc(time_steps*sizeof(double complex));

  int i;
  double last_time = (float)10/sample_frequency;
  for(i = 0; i < end_time*sample_frequency; i++){
    data.chirp_time_vector[i] = last_time;
    last_time += (float)1/sample_frequency;
  }

  for(i = 0; i < time_steps; i++){
    data.chirp_signal[i] = cos( 2*PI*(variables.start_frequency*data.chirp_time_vector[i]+chirp_rate*data.chirp_time_vector[i]*data.chirp_time_vector[i] ))+_Complex_I*sin( 2*PI*(variables.start_frequency*data.chirp_time_vector[i]+chirp_rate*data.chirp_time_vector[i]*data.chirp_time_vector[i] ));
  }
}

void chirp_matched_generator(){
  float end_time = variables.btproduct/variables.bandwidth;
  float chirp_rate = variables.bandwidth/end_time;
  
  unsigned long int sample_frequency = 5*variables.bandwidth;
  unsigned long int time_steps = end_time*sample_frequency;
  variables.chirp_length = time_steps;

  data.matched_time_vector = malloc(time_steps*sizeof(double));
  data.matched_chirp = malloc(time_steps*sizeof(double complex));

  int i;
  double last_time = end_time+(float)10/sample_frequency;
  for(i = 0; i < end_time*sample_frequency; i++){
    data.matched_time_vector[i] = last_time;
    last_time -= (float)1/sample_frequency;
  }

  for(i = 0; i < end_time*sample_frequency; i++){
    data.matched_chirp[i] = cos( 2*PI*(variables.start_frequency*data.matched_time_vector[i]+chirp_rate*data.matched_time_vector[i]*data.matched_time_vector[i] ))-_Complex_I*sin( 2*PI*(variables.start_frequency*data.chirp_time_vector[i]+chirp_rate*data.chirp_time_vector[i]*data.chirp_time_vector[i] ));
  }
}

void pulse_compress_signal(){
  data.pulse_compressed_waveform = malloc(variables.chirp_length*sizeof(double complex));
  int kernel_length = variables.chirp_length;
  unsigned int filter_length = 2*kernel_length;

  fftw_complex* padded_signal = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* padded_kernel = fftw_malloc(filter_length*sizeof(fftw_complex));
  memset(padded_signal, 0, filter_length*sizeof(fftw_complex));
  memset(padded_kernel, 0, filter_length*sizeof(fftw_complex));
  memcpy(padded_signal, data.chirp_signal, kernel_length*sizeof(fftw_complex));
  memcpy(padded_kernel, data.matched_chirp, kernel_length*sizeof(fftw_complex));

  fftw_complex* sigfft = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* matchfft = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* product = fftw_malloc(filter_length*sizeof(fftw_complex));

  fftw_plan sigp = fftw_plan_dft_1d(filter_length, padded_signal, sigfft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan matchp = fftw_plan_dft_1d(filter_length, padded_kernel, matchfft, FFTW_FORWARD, FFTW_ESTIMATE);
  data.pulse_compressed_waveform = malloc(kernel_length*sizeof(double complex));
  fftw_plan iff = fftw_plan_dft_1d(kernel_length, product, data.pulse_compressed_waveform, FFTW_BACKWARD, FFTW_ESTIMATE);

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
  data.pulse_compressed_radar_image = malloc(variables.nrows*variables.ncols*sizeof(double complex));
  
  // Make sure input has valid values.
  int kernel_length = variables.chirp_length;
  int z;
  for(z = 0; z < kernel_length; z++){
    if(isnan(data.matched_chirp[z]))
      data.matched_chirp[z] = 0;
  }

  unsigned int filter_length = variables.nrows + kernel_length;
  filter_length = pow(2, ceil(log(filter_length)/log(2)));

  fftw_complex* padded_kernel = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* kernel_fft = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* product = fftw_malloc(filter_length*sizeof(fftw_complex));
  memset(padded_kernel, 0, filter_length*sizeof(fftw_complex));
  memset(kernel_fft, 0, filter_length*sizeof(fftw_complex));
  memcpy(padded_kernel, data.matched_chirp, kernel_length*sizeof(fftw_complex));

  // Compute fft of filter kernel.
  fftw_plan kernelfft = fftw_plan_dft_1d(filter_length, padded_kernel, kernel_fft, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(kernelfft);
  fftw_complex* output_column;
  fftw_complex* padded_column = fftw_malloc(filter_length*sizeof(fftw_complex));
  fftw_complex* padded_column_fft = fftw_malloc(filter_length*sizeof(fftw_complex));
  int i,j,k;
  for(i = 0; i < variables.ncols; i++){
    double complex* column = &data.radar_image[i*variables.nrows];

    // Make sure we have valid values.
    for(k = 0; k < variables.nrows; k++){
	if(isnan(column[k]))
	  column[k] = 0;
    }

    output_column = &data.pulse_compressed_radar_image[i*variables.nrows];
    memset(padded_column, 0, filter_length*sizeof(fftw_complex));
    memcpy(padded_column, column, variables.nrows*sizeof(fftw_complex));
    memset(padded_column_fft, 0, filter_length*sizeof(fftw_complex));
    fftw_plan colfft = fftw_plan_dft_1d(filter_length, padded_column, padded_column_fft, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(colfft);


    for(j = 0; j < filter_length; j++){
	product[j] = padded_column_fft[j]*kernel_fft[j];
    }

    fftw_plan colifft = fftw_plan_dft_1d(variables.nrows, product, output_column, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(colifft);

    fftw_destroy_plan(colfft);
    fftw_destroy_plan(colifft);

    // Re-normalize ifft output.
    for(j = 0; j < variables.nrows; j++)
      output_column[j] /= variables.nrows;
  }

  fftw_free(product);
  fftw_free(padded_column_fft);
  fftw_free(padded_column);
  fftw_free(padded_kernel);

  fftw_destroy_plan(kernelfft);
}

void insert_waveform_in_scene(){
	data.scene_with_waveform = malloc(variables.ncols*variables.nrows*sizeof(double complex));
	int i;
	for(i = 0; i < variables.chirp_length; i++){
	  data.scene_with_waveform[ (int)(variables.ncols/2)*variables.nrows + (int)(variables.nrows/2) -variables.chirp_length/2 +i] = data.chirp_signal[i];
	}
}

void radar_imager(){
  data.radar_image = malloc(variables.nrows*variables.ncols*sizeof(double complex));
  double azimuth_coverage = round(variables.altitude*tan(0.5*variables.beamwidth));
  unsigned int beamcrossrange = variables.chirp_length*azimuth_coverage/variables.signal_distance;

  unsigned int i,j,k;
  int beam_value;
  for(i = 0; i < variables.ncols; i++){
    for(j = 0; j < 2*beamcrossrange; j++){
      beam_value = j-beamcrossrange;
      for(k = 0; k < variables.nrows; k++){
	if(i+beam_value <= variables.ncols){
	  int dist = sqrt(pow(beam_value,2)+pow(k, 2));
	  if(dist <= variables.nrows){
	    if(i+beam_value >= 1){
	      data.radar_image[i*variables.nrows+dist] += data.scene_with_waveform[(i+beam_value)*variables.nrows+k]/pow(dist,2);
	    }
	  }
	}
      }
    }
  }

}

void gbp(){
  data.sar_image = malloc(variables.nrows*variables.ncols*sizeof(double complex));
  int j,k,l;
  for(j = 0; j < variables.ncols; j++){
    for(k = 0; k < variables.nrows; k++){
      for(l = 0; l < variables.ncols; l++){
	unsigned int range_index = sqrt((l-j)*(l-j)+k*k);
	if(range_index < variables.nrows){
	  data.sar_image[j*variables.nrows+k] += data.pulse_compressed_radar_image[l*variables.nrows+range_index]/variables.nrows;
	}
      }
    }
  }

}

void gbp_fft(){
  data.sar_fft = malloc(variables.nrows*variables.ncols*sizeof(double complex));
  data.sar_img_shifted = malloc(variables.nrows*variables.ncols*sizeof(double complex));
  int i, j;
  for(i = 0; i < variables.ncols; i++){
    for(j = 0; j < variables.nrows; j++){
	data.sar_img_shifted[variables.nrows*i+j] = data.sar_image[variables.nrows*i+j]*pow(-1,i+j);
    }
  }

  fftw_plan fft = fftw_plan_dft_2d(variables.ncols, variables.nrows, data.sar_img_shifted, data.sar_fft, FFTW_FORWARD, FFTW_ESTIMATE);
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
  for(i = variables.ncols/2 - 5; i < variables.ncols/2+5; i++){
    for(j = variables.nrows/2 - 5; j < variables.nrows/2 + 5; j++){
	data.sar_fft[i*variables.nrows+j] = 0;
    }
  }

  fftw_plan ifft = fftw_plan_dft_2d(variables.ncols, variables.nrows, data.sar_fft, data.sar_image, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(ifft);
  fftw_destroy_plan(ifft);
}
