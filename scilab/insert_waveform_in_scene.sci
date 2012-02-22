function inserted = insert_waveform_in_scene(waveform, scene)
    [nrows, ncols] = size(scene);
    wavesize = size(waveform, '*');
    inserted = zeros(nrows, ncols);
    
    for i = [1:wavesize]
        inserted(nrows/2-wavesize/2+i, ncols/2) = waveform(1, i);
    end
    
endfunction