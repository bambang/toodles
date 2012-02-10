function [sar_image, sar_fft, pulse_compressed, raw_data, waveform_in_scene] = sar_simulator(rangesize, crossrange_size, reflector_positions, start_frequency, bandwidth, platform_altitude, antenna_beamwidth)
    exec('pulse_compress.sci');
    exec('scene_generator.sci');
    exec('chirp.sci');
    exec('data_generator.sci');
    exec('chirp_matched.sci');
    exec('gbp.sci');
    exec('compress_raw_data.sci');
    
    timer();
    scene = scene_generator(rangesize, crossrange_size, reflector_positions);
    printf("Scene generation took %f s.\n", timer());
    [chirp_waveform, waveform_distance,waveform_time] = chirp(start_frequency, bandwidth);
    printf("Waveform generation took %f s.\n",timer());
    matched_chirp_waveform = chirp_matched(start_frequency, bandwidth);
    printf("Matched waveform generation took %f s.\n", timer());
    compressed_waveform = convol(chirp_waveform,matched_chirp_waveform);
    //compressed_waveform = pulse_compress(chirp_waveform, matched_chirp_waveform);
    printf("Pulse compression took %f s.\n",timer());
    [raw_data,waveform_in_scene] = data_generator(scene, antenna_beamwidth, platform_altitude, chirp_waveform);
    printf("Raw radar data generation took %f s.\n",timer());
    pulse_compressed = compress_raw_data(raw_data, matched_chirp_waveform);
    printf("Pulse compression of raw radar data took %f s.\n",timer());
    [sar_image,sar_fft] = gbp(pulse_compressed);
    printf("GBP and its FFT took %f s.\n",timer());
    
    // Calculations needed for plotting.
    [scene_rows, scene_cols] = size(scene);
    [image_range,image_crossrange] = size(sar_image);
    waveform_time = size(chirp_waveform, '*');
    distance_per_bin = waveform_distance/waveform_time;
    range_distance = distance_per_bin*image_range;
    crossrange_distance = distance_per_bin*image_crossrange;
    printf("Distance per bin: %f\n", distance_per_bin);
    printf("Range distance: %f\n", range_distance);
    printf("Crossrange distance: %f\n",crossrange_distance);
    
    range_vector = [0:distance_per_bin:range_distance-distance_per_bin];
    crossrange_vector = [0:distance_per_bin:crossrange_distance-distance_per_bin];
    
    printf("Range size: %i\n",size(range_vector,'*'));
    printf("Crossrange size: %i\n", size(crossrange_vector,'*'));
    [ix,iy] = size(raw_data);
    printf("Image range: %i Image crossrange: %i\n", ix, iy);
    
    // Plot the original scene.
    da = gda();
    da.x_label.text = "Range";
    da.y_label.text = "Crossrange";
    figure("Figure_name", "Scene");
    plot3d([1:scene_rows],[1:scene_cols], scene);
    
    // Plot the waveform.
    da = gda();
    da.x_label.text = "Time";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Waveform");
    [xw,yw] = size(chirp_waveform);
    plot(chirp_waveform(1,:));
    
    // Plot fft of waveform.
    da = gda();
    da.x_label.text = "Frequency";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Waveform FFT");
    plot(abs(fftshift(fft(chirp_waveform(1,:)))));
    
    // Plot pulse-compressed waveform.
    da = gda();
    da.x_label.text = "Time";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Pulse-compressed waveform");
    [xs,ys] = size(compressed_waveform);
    plot(compressed_waveform(1,:));
    
    // Plot fft of pulse-compressed waveform.
    da = gda();
    da.x_label.text = "Frequency";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Pulse-compressed waveform FFT");
    plot(abs(fftshift(fft(compressed_waveform(1,:)))));
    
    // Plot scene with waveform.
    da = gda();
    da.x_label.text = "Range (m)";
    da.y_label.text = "Crossrange (m)";
    figure("Figure_name", "Scene with waveform");
    //mesh([1:distance_per_bin:range_distance], [1:distance_per_bin:crossrange_distance],distance_per_bin*abs(waveform_in_scene));
    plot3d(range_vector, crossrange_vector, waveform_in_scene);
    
    // Plot radar image.
    [rrows,rcols] = size(raw_data);
    da = gda();
    da.x_label.text = "Range (m)";
    da.y_label.text = "Crossrange (m)";
    figure("Figure_name", "Radar image");
    //mesh([1:distance_per_bin:range_distance], [1:distance_per_bin:crossrange_distance],distance_per_bin*abs(raw_data));
    plot3d(range_vector, crossrange_vector, abs(raw_data));
    
    // Plot compressed radar image.
    da = gda();
    da.x_label.text = "Range (m)";
    da.y_label.text = "Crossrange (m)";
    figure("Figure_name", "Pulse-compressed radar image");
    //mesh([1:distance_per_bin:range_distance], [1:distance_per_bin:crossrange_distance],distance_per_bin*abs(pulse_compressed));
    plot3d(range_vector, crossrange_vector, abs(pulse_compressed));
    
    // Plot backprojected image.
    da = gda();
    da.x_label.text = "Range (m)";
    da.y_label.text = "Crossrange (m)";
    figure("Figure_name", "Backprojected image");
    plot3d(range_vector, crossrange_vector, abs(sar_image))
    //mesh([1:distance_per_bin:range_distance], [1:distance_per_bin:crossrange_distance],distance_per_bin*abs(sar_image));
    
    // Plot fft of backprojected image.
    da = gda();
    da.x_label.text = "Crossrange frequency";
    da.y_label.text = "Range frequency";
    figure("Figure_name", "Backprojected image 2D-FFT");
    mesh(abs(sar_fft));
    //[frow,fcol] = size(sar_fft);
    //plot3d([1:frow],[1:fcol],sar_fft);
    //da = gda();
    //da.x_label.text = "Range freq";
    //da.y_label.text = "Crossrange freq";
endfunction