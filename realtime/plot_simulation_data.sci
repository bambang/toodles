function plot_simulation_data()
    
    dimensions = read("dimensions.dat", 6, 1);
    chirp_length = dimensions(1);
    scene_rows = dimensions(2);
    scene_cols = dimensions(3);
    radar_rows = dimensions(4);
    radar_cols = dimensions(5);
    signal_distance = dimensions(6);
    
    distance_step = signal_distance/chirp_length;
    end_row_distance = distance_step*radar_rows;
    end_col_distance = distance_step*radar_cols;
    
    printf("Chirp length: %i\n", chirp_length);
    printf("Scene Rows: %i Cols: %i\n", scene_rows, scene_cols);
    printf("Radar Rows: %i Cols: %i\n", radar_rows, radar_cols);
    
    scene_data = read("scene.dat", scene_cols, scene_rows);
    chirp_data = read("chirp.dat", chirp_length, 2);
    matched_chirp_data = read("matched.dat", chirp_length, 2);
    pulse_compressed_data = read("compressed.dat", chirp_length, 2);
    scene_with_waveform_data = read("scene_with_waveform.dat", radar_cols, 2*radar_rows);
    radar_image_data = read("radar_image.dat", radar_cols, 2*radar_rows);
    compressed_image_data = read("pulse_compressed_image.dat", radar_cols, 2*radar_rows);
    sar_image_data = read("sar_image.dat", radar_cols, 2*radar_rows);
    sar_fft_data = read("sar_fft.dat", radar_cols, 2*radar_rows);
    
    chirp_fft = read("chirpfft.dat", chirp_length, 2);
    matched_fft = read("matchedfft.dat", chirp_length, 2);
    
    complex_chirp_data = chirp_data(:,1)+%i*chirp_data(:,2);
    complex_matched_chirp = matched_chirp_data(:,1)+%i*matched_chirp_data(:,2);
    complex_pulse_compressed = pulse_compressed_data(:,1)+%i*pulse_compressed_data(:,2);
    complex_waveform_scene_data = zeros(radar_rows, radar_cols);
    complex_image_data = zeros(radar_rows, radar_cols);
    complex_compressed_data = zeros(radar_rows, radar_cols);
    complex_sar_data = zeros(radar_rows, radar_cols);
    complex_sar_fft = zeros(radar_rows, radar_cols);
    
    complex_chirp_fft = chirp_fft(:,1)+%i*chirp_fft(:,2);
    complex_matched_fft = matched_fft(:,1)+%i*matched_fft(:,2);
    
    for col = [1:radar_cols]
        for row = [1:radar_rows]
            complex_waveform_scene_data(row,col) = scene_with_waveform_data(col, 2*row-1) + %i*scene_with_waveform_data(col, 2*row);
            complex_image_data(row,col) =  radar_image_data(col,2*row-1) + %i*radar_image_data(col, 2*row);
            complex_compressed_data(row,col) = compressed_image_data(col, 2*row-1) + %i*compressed_image_data(col, 2*row);
            complex_sar_data(row,col) = sar_image_data(col,2*row-1) + %i*radar_image_data(col, 2*row);
            complex_sar_fft(row,col) = sar_fft_data(col, 2*row-1) + %i*sar_fft_data(col, 2*row);
        end
    end
    
    da = gda();
    da.x_label.text = "Time";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Chirp signal");
    plot(complex_chirp_data);
    
    da = gda();
    da.x_label.text = "Frequency";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Chirp FFT");
    plot(abs(complex_chirp_fft));
    
    da = gda();
    da.x_label.text = "Time";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Matched Chirp signal");
    plot(complex_matched_chirp);
    
    da = gda();
    da.x_label.text = "Frequency";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Matched Chirp FFT");
    plot(abs(complex_matched_fft));
    
    da = gda();
    da.x_label.text = "Time";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Pulse-compressed signal");
    plot(complex_pulse_compressed);
    
    da = gda();
    da.x_label.text = "Range (m)";
    da.y_label.text = "Crossrange (m)";
    figure("Figure_name", "Scene with waveform");
    plot3d([0:distance_step:end_row_distance-distance_step], [0:distance_step:end_col_distance-distance_step], abs(complex_waveform_scene_data));
    
    da = gda();
    da.x_label.text = "Range (m)";
    da.y_label.text = "Crossrange (m)";
    figure("Figure_name", "Radar image");
    plot3d([0:distance_step:end_row_distance-distance_step], [0:distance_step:end_col_distance-distance_step], abs(complex_image_data));
    
    da = gda();
    da.x_label.text = "Range (m)";
    da.y_label.text = "Crossrange (m)";
    figure("Figure_name", "Pulse-compressed Radar image");
    plot3d([0:distance_step:end_row_distance-distance_step], [0:distance_step:end_col_distance-distance_step], abs(complex_compressed_data));
    
    da = gda();
    da.x_label.text = "Range (m)";
    da.y_label.text = "Crossrange (m)";
    figure("Figure_name", "SAR image");
    plot3d([0:distance_step:end_row_distance-distance_step], [0:distance_step:end_col_distance-distance_step], abs(complex_sar_data));
//    
//    da = gda();
//    da.x_label.text = "Range frequency";
//    da.y_label.text = "Crossrange frequency";
//    figure("Figure_name", "SAR image FFT");
//    mesh(abs(fftshift(fft2(complex_sar_data))));
    
    da = gda();
    da.x_label.text = "Range frequency";
    da.y_label.text = "Crossrange frequency";
    figure("Figure_name", "SAR image FFT");
    mesh(abs(complex_sar_fft));
    
endfunction