function [corr_h, corr_v, corr_d] = calculate_correlation(img)
%
% Calculates correlation coefficient for adjacent pixels
% (horizontal, vertical, diagonal).
% Ideal value (randomness) is close to 0.
%

    % Set number of pixel samples for test
    num_samples = 5000;
    
    [rows, cols] = size(img);
    img = double(img); % Required for corrcoef

    % --- Horizontal ---
    % Select random pixels, avoiding the last column
    x_coords = randi([1 rows], num_samples, 1);
    y_coords = randi([1 cols - 1], num_samples, 1);
    % Get values
    v1 = img(sub2ind([rows, cols], x_coords, y_coords));
    v2 = img(sub2ind([rows, cols], x_coords, y_coords + 1));
    % Calculate correlation
    C_h = corrcoef(v1, v2);
    corr_h = C_h(1, 2);

    % --- Vertical ---
    % Select random pixels, avoiding the last row
    x_coords = randi([1 rows - 1], num_samples, 1);
    y_coords = randi([1 cols], num_samples, 1);
    % Get values
    v1 = img(sub2ind([rows, cols], x_coords, y_coords));
    v2 = img(sub2ind([rows, cols], x_coords + 1, y_coords));
    % Calculate correlation
    C_v = corrcoef(v1, v2);
    corr_v = C_v(1, 2);

    % --- Diagonal ---
    % Select random pixels, avoiding the last row and column
    x_coords = randi([1 rows - 1], num_samples, 1);
    y_coords = randi([1 cols - 1], num_samples, 1);
    % Get values
    v1 = img(sub2ind([rows, cols], x_coords, y_coords));
    v2 = img(sub2ind([rows, cols], x_coords + 1, y_coords + 1));
    % Calculate correlation
    C_d = corrcoef(v1, v2);
    corr_d = C_d(1, 2);
    
end
