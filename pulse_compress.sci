function compressed_signal = pulse_compress(signal, matched_signal)
    // This function calculates the pulse compression of a signal with its matched response.
    fftsignal = fft(signal);
    fftmatched_signal = fft(matched_signal);
    compressed_signal = ifft(fftsignal.*fftmatched_signal);
endfunction