function [data,waveform_in_scene] = data_generator(scene, beamwidth, altitude, waveform)
    // This function takes an existing generated scene, the antenna beamwidth,
    // the height of the SAR platform above the scene, and the waveform as input parameters.
    // It then generates the reflected radar data, based on the location of objects in the scene
    // and the antenna beamwidth.
    // The beamwidth causes the radar response from a single point to be smeared across a region,
    // and takes the form of a hyperbola.
    
    wtime = size(waveform, '*');
    [srows, scols] = size(scene);
    
    timer();
    
    // Insert waveform in scene.
    waveform_in_scene = zeros(wtime*srows, scols);
    for k = [1:scols]
        for i = [1:srows]
            if scene(i,k) ~=0
                for j = [1:wtime]
                    waveform_in_scene((i-1)*wtime + j,k) = waveform(1,j).*1000000;
                end
            end
        end
    end
    
    printf("Inserting the waveform in the scene took %f s.\n",timer());
    
    [prows,pcols] = size(waveform_in_scene);
    
    // 1) Calculate how many crossrange bins are illuminated by the beam at this crossrange position.
    // 2) For all points within the illuminated area, calculate the distance to the radar.
    // 3) Project these distances onto range distances
    // 4) Reduce the amplitude as needed
    // 5) Store results in the crossrange bin.
    // 6) Repeat until all range bins have been accounted for.
    
    // Calculate how many crossrange bins are illuminated by the beam.
    beamcrossrange = round(altitude*tan(beamwidth));
    printf("Antenna beam covers %i crossrange bins.\n", beamcrossrange);
    
    data = zeros(prows,pcols);

    timer();
    
    for col = 1:pcols
        for beamrange = [-beamcrossrange:beamcrossrange]
            for row = 1:prows
                if (col+beamrange) <= pcols
                    dist = sqrt(beamrange^2+row^2);
                    if dist <= prows
                        if col+beamrange >= 1
                            data(dist,col) = data(dist,col) + waveform_in_scene(row,col+beamrange)/dist^2;
                        end
                    end
                end
            end
        end
    end
    
    printf("Col compression took %f s.\n",timer());
    
endfunction