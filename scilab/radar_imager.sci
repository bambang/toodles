function radar_image = radar_imager(swwf, beamwidth, altitude)
    [prows,pcols] = size(swwf);
    beamcrossrange = round(altitude*tan(beamwidth));
    
    radar_image = zeros(prows,pcols);
    
    for col = 1:pcols
        for beamrange = [-beamcrossrange:beamcrossrange]
            for row = 1:prows
                if (col+beamrange) <= pcols
                    dist = sqrt(beamrange^2+row^2);
                    if dist <= prows
                        if col+beamrange >= 1
                            radar_image(dist,col) = radar_image(dist,col) + swwf(row,col+beamrange)/dist^2;
                        end
                    end
                end
            end
        end
    end
    
endfunction