function chirp_matched_vector = chirp_matched(start_frequency, bandwidth)
    chirp_rate = 100;
    end_time = bandwidth/chirp_rate;
    sample_frequency = 5*bandwidth;
    chirp_time_vector = [-end_time:1/sample_frequency:0];
    
    chirp_matched_vector = %e^(%i*2*%pi*(start_frequency.*chirp_time_vector+chirp_rate.*chirp_time_vector.*chirp_time_vector));
endfunction