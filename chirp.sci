function [chirp_signal,waveform_distance,chirp_time_vector] = chirp(start_frequency, bandwidth)
    // Generate chirp.
    // bandwidth = chirp_rate*time
    //btproduct = 1;
    //end_time_chirp = btproduct/bandwidth;
    //chirp_rate = bandwidth/end_time_chirp;
    //sample_frequency = 1.4*bandwidth;
    //sample_step = end_time_chirp/sample_frequency;
    //chirp_time_vector = [0:sample_step:end_time_chirp];
    //phase = 0;
    //end_frequency = start_frequency + chirp_rate*chirp_time_vector($)^2;
    //chirp_signal = amplitude*%e^(-%i*2*%pi*(phase+start_frequency*chirp_time_vector+chirp_rate*chirp_time_vector^2));
    
    // Generate pulse based on BT-product.
    //end_time_pulse = btproduct/bandwidth;
    //chirp_size = size(chirp_signal,'*');
    //pulse_signal = zeros(1,end_time_pulse*sample_frequency);
    //pulse_signal(1:chirp_size) = chirp_signal;
    //time_vector = [0:sample_step:end_time_pulse];

    //waveform_distance = end_time_chirp*3*10^8;
    
    //printf("Chirp rate: %i\n",chirp_rate);
    //printf("Number of samples in chirp: %i\n",end_time_chirp/sample_step);
    //printf("Number of samples in pulse: %i\n", end_time_pulse/sample_step);
    //printf("Start frequency: %i\n",start_frequency);
    //printf("End frequency: %f\n", end_frequency);
    //printf("Chirp end time: %f\n", end_time_chirp);
    //printf("Pulse end time: %f\n", end_time_pulse);
    //printf("Sampling rate: %f\n",sample_frequency);
    //printf("Sample step size: %f\n",sample_step);
    //printf("Distance of a waveform: %f\n", waveform_distance);
    
    chirp_rate = 100;
    end_time = bandwidth/chirp_rate;
    sample_frequency = 5*bandwidth;
    chirp_time_vector = [0:1/sample_frequency:end_time];
    chirp_signal = %e^(-%i*2*%pi*(start_frequency.*chirp_time_vector+chirp_rate.*chirp_time_vector.*chirp_time_vector));
    waveform_distance = (1/start_frequency)*3*10^8;
    
    printf("Chirp rate: %i\n",chirp_rate);
    printf("Start frequency: %f\n",start_frequency);
    printf("End frequency: %f\n", start_frequency+bandwidth);
    printf("Chirp end time: %f\n", chirp_time_vector($));
    printf("Sampling rate: %f\n",sample_frequency);
    printf("Distance of a waveform: %f\n", waveform_distance);
endfunction