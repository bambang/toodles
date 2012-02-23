function [chirp, matched, compressed_signal, swwf, radar_image, pc_image, sar_image, sar_fft, srows, scols] = sar_simulator()
    stacksize('max');
    exec('chirp_generator.sci');
    exec('matched_generator.sci');
    exec('compress_signals.sci');
    exec('insert_waveform_in_scene.sci');
    exec('radar_imager.sci');
    exec('pulse_compress_image.sci');
    exec('gbp.sci');
    
    bt = input("Bandwidth-time product: ");
    start_freq = input("Chirp start frequency: ");
    bandwidth = input("Chirp bandwidth: ");
    beamwidth = input("Antenna azimuth beamwidth: ");
    
    [chirp, waveform_distance] = chirp_generator(bt,start_freq, bandwidth);
    matched = matched_generator(bt, start_freq, bandwidth);
    compressed_signal = compress_signals(chirp, matched);
    
    printf("Signal waveform covers %f meters.\n", waveform_distance);
    srows = input("Enter area azimuth length (m): ")*size(chirp)/waveform_distance;
    scols = input("Enter area range (m): ")*size(chirp)/waveform_distance;
    altitude = input("Antenna altitude (m): ");
    
    scene = zeros(srows(1,2), scols(1,2));
    swwf = insert_waveform_in_scene(chirp, scene);

    radar_image = radar_imager(swwf, beamwidth, altitude);

    pc_image = pulse_compress_image(radar_image, matched);

    sar_image = gbp(pc_image);

    sar_fft = fftshift(fft2(sar_image));
    
    da = gda();
    da.x_label.text = "Time";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Chirp");
    plot(chirp)
    
    da = gda();
    da.x_label.text = "Frequency";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Chirp FFT");
    plot(abs(fft(chirp)));
    
    da = gda();
    da.x_label.text = "Time";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Matched chirp");
    plot(matched);
    
    da = gda();
    da.x_label.text = "Frequency";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Matched chirp FFT");
    plot(abs(fft(matched)));
    
    da = gda();
    da.x_label.text = "Time";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "Pulse-compression of chirp and its match");
    plot(compressed_signal);
    
    da = gda();
    da.x_label.text = "Frequency";
    da.y_label.text = "Amplitude";
    figure("Figure_name", "FFT of pulse-compression of chirp and its match");
    plot(abs(fft(compressed_signal)));
    
    da = gda();
    da.y_label.text = "Range";
    da.x_label.text = "Crossrange";
    figure("Figure_name", "Scene");
    mesh(abs(swwf));
    
    da = gda();
    da.y_label.text = "Range";
    da.x_label.text = "Crossrange";
    figure("Figure_name", "Radar image");
    mesh(abs(radar_image));
    
    da = gda();
    da.y_label.text = "Range";
    da.x_label.text = "Crossrange";
    figure("Figure_name", "Pulse-compressed radar image");
    mesh(abs(pc_image));
    
    da = gda();
    da.y_label.text = "Range";
    da.x_label.text = "Crossrange";
    figure("Figure_name", "SAR image");
    mesh(abs(sar_image));
    
    da = gda();
    da.y_label.text = "Range frequency";
    da.x_label.text = "Crossrange frequency";
    figure("Figure_name", "SAR image FFT");
    mesh(abs(sar_fft));
    
endfunction