function pcimage = pulse_compress_image(image, matched)
    [nrows, ncols] = size(image);
    pcimage = zeros(nrows, ncols);
    
    for i = [1:ncols]
       con = convol(matched, image(:,i));
       pcimage(:,i) = con(1,1:size(pcimage,1))';
    end
    
endfunction