function compressed = compress_signals(one, two)
    //compressed = ifft(fft(one).*fft(two));
    compressed = convol(one,two);
endfunction