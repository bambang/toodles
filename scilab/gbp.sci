function projected_image = gbp(radar_image)
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
endfunction