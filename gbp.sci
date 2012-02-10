function [projected_image,projected_fft] = gbp(radar_image)
    // This is the Global Backprojection time-domain SAR imaging algorithm.
    // The input is compressed radar data in a 2D matrix.
    // A 2D projected matrix is returned by the function.
    [range_pixels, crossrange_pixels] = size(radar_image);
    projected_image = zeros(range_pixels, crossrange_pixels);
    
    for j = [1:crossrange_pixels]
        for k = [1:range_pixels]
            for l = [1:crossrange_pixels]
               range_index = sqrt((l-j)^2+k^2);
               if range_index <= range_pixels
                   projected_image(k,j) = projected_image(k,j) + radar_image(range_index,l);
               end
            end
        end
    end
    
    projected_fft = fftshift(fft2(projected_image));
    
endfunction