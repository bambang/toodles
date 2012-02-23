function projected_image = gbp(radar_image)
    // This is the Global Backprojection time-domain SAR imaging algorithm.
    // The input is compressed radar data in a 2D matrix.
    // A 2D projected matrix is returned by the function.
    [range_pixels, crossrange_pixels] = size(radar_image);
    for j = [0:crossrange_pixels]
        for k = [0:range_pixels]
            for l = [0:crossrange_pixels]
               projected_image[k][j] = projected_image[k][j] + radar_image[l][sqrt(l-j)^2-k^2];
            end
        end
    end
endfunction