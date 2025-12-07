function C_attacked = apply_attack(C_orig, attack_type, intensity, img_tamper)
%
% Symuluje "Alteration Attack" na obrazie szyfrogramu.
% WYMAGA: Image Processing Toolbox (dla imnoise)
%
% C_orig: Oryginalny szyfrogram (uint8)
% attack_type: 'crop', 'tamper', 'salt_pepper', 'gaussian'
% intensity: Poziom ataku (np. 0.1 dla 10%)
% img_tamper: Obcy obraz do użycia przy ataku 'tamper'
%

    C_attacked = C_orig;
    [rows, cols] = size(C_orig);
    num_pixels = rows * cols;

    switch attack_type
        case 'crop'
            % Atak "Cropping": Zastąp 'intensity' % pikseli zerami
            num_attacked_pixels = round(num_pixels * intensity);
            % Wybierz losowe piksele do "usunięcia"
            idx_to_crop = randperm(num_pixels, num_attacked_pixels);
            C_attacked(idx_to_crop) = 0; % Zastąp zerami

        case 'tamper'
            % Atak "Tampering": Zastąp 'intensity' % pikseli innym obrazem
            if nargin < 4 || isempty(img_tamper)
                error('Atak "tamper" wymaga czwartego argumentu: img_tamper');
            end
            
            num_attacked_pixels = round(num_pixels * intensity);
            % Wybierz losowe piksele do zastąpienia w szyfrogramie
            idx_to_tamper = randperm(num_pixels, num_attacked_pixels);
            % Wybierz losowe piksele z obrazu-intruza
            idx_from_tamper = randperm(numel(img_tamper), num_attacked_pixels);
            
            C_attacked(idx_to_tamper) = img_tamper(idx_from_tamper);
            
        case 'salt_pepper'
            % Atak szumem "Sól i Pieprz"
            % 'intensity' to gęstość szumu (0.01 do 1.0)
            C_attacked = imnoise(C_orig, 'salt & pepper', intensity);

        case 'gaussian'
            % Atak szumem Gaussian
            % 'intensity' to wariancja szumu (np. 0.01 do 0.1)
            % Średnia (mean) jest 0, zgodnie z tekstem
            C_attacked = imnoise(C_orig, 'gaussian', 0, intensity);
            
        otherwise
            error('Nieznany typ ataku: %s', attack_type);
    end
end