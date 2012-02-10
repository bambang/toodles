function scene = scene_generator(rows, cols, objects)
    // This function generates a 2D scene matrix with a number of objects.
    // The objects span a number of pixels each to emulate object radius.
    [orows, ocols] = size(objects);
    
    scene = zeros(rows, cols);
    
    for i = [1:orows]
        row = objects(i,1);
        col = objects(i,2);
        
        scene(row,col) = 1000;
    end
endfunction