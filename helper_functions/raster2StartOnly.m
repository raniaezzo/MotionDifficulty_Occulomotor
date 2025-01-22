function startMSMatrix = raster2StartOnly(raster_matrix)

% Assuming your original matrix is called 'raster_matrix'
% Initialize the new matrix as a copy of the original
startMSMatrix = raster_matrix;

% Loop through each row and check the condition for each column
for row = 1:size(raster_matrix, 1)
    for col = 2:size(raster_matrix, 2) % Start from column 2 since we need to check the previous column
        if raster_matrix(row, col - 1) == 1
            startMSMatrix(row, col) = 0; % Set current column to 0 if the previous column is 1
        end
    end
end

end