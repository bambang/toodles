function compressed_data = compress_raw_data(raw_data, matched_waveform)
    [nrows, ncols] = size(raw_data);
    wavesize = size(matched_waveform, '*');
    
    compressed_data = zeros(nrows, ncols);

    for i = [1:ncols]
       con = convol(matched_waveform, raw_data(:,i));
       compressed_data(:,i) = con(1,1:size(compressed_data,1))';
    end
    
endfunction