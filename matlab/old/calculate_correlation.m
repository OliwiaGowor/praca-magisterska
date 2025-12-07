function [corr_h, corr_v, corr_d] = calculate_correlation(img)
%
% Oblicza współczynnik korelacji dla sąsiednich pikseli
% (poziomo, pionowo, diagonalnie).
% Wartość idealna (losowość) jest bliska 0.
%

    % Ustaw liczbę próbek pikseli do testu
    num_samples = 5000;
    
    [rows, cols] = size(img);
    img = double(img); % Wymagane dla corrcoef

    % --- Poziomo (Horizontal) ---
    % Wybierz losowe piksele, unikając ostatniej kolumny
    x_coords = randi([1 rows], num_samples, 1);
    y_coords = randi([1 cols - 1], num_samples, 1);
    % Pobierz wartości
    v1 = img(sub2ind([rows, cols], x_coords, y_coords));
    v2 = img(sub2ind([rows, cols], x_coords, y_coords + 1));
    % Oblicz korelację
    C_h = corrcoef(v1, v2);
    corr_h = C_h(1, 2);

    % --- Pionowo (Vertical) ---
    % Wybierz losowe piksele, unikając ostatniego wiersza
    x_coords = randi([1 rows - 1], num_samples, 1);
    y_coords = randi([1 cols], num_samples, 1);
    % Pobierz wartości
    v1 = img(sub2ind([rows, cols], x_coords, y_coords));
    v2 = img(sub2ind([rows, cols], x_coords + 1, y_coords));
    % Oblicz korelację
    C_v = corrcoef(v1, v2);
    corr_v = C_v(1, 2);

    % --- Diagonalnie (Diagonal) ---
    % Wybierz losowe piksele, unikając ostatniego wiersza i kolumny
    x_coords = randi([1 rows - 1], num_samples, 1);
    y_coords = randi([1 cols - 1], num_samples, 1);
    % Pobierz wartości
    v1 = img(sub2ind([rows, cols], x_coords, y_coords));
    v2 = img(sub2ind([rows, cols], x_coords + 1, y_coords + 1));
    % Oblicz korelację
    C_d = corrcoef(v1, v2);
    corr_d = C_d(1, 2);
    
end