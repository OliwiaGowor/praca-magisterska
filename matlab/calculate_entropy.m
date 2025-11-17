function H = calculate_entropy(img)
%
% Oblicza entropię Shannona dla obrazu (skala szarości).
% WERSJA 2: Usunięto 'isuint8' dla lepszej kompatybilności.
%

    % Zawsze konwertuj na uint8, aby zagwarantować poprawny typ dla histcounts
    % Ta linia zastępuje sprawdzanie 'isuint8' i naprawia błąd.
    img = uint8(img); 

    % 1. Zlicz wystąpienia
    counts = histcounts(img, -0.5:1:255.5);
    
    % 2. Oblicz prawdopodobieństwo
    total_pixels = numel(img);
    probabilities = counts / total_pixels;
    
    % 3. Oblicz entropię
    probabilities = probabilities(probabilities > 0); % Usuń p=0
    H = -sum(probabilities .* log2(probabilities));

end