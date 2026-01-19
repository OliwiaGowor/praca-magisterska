function H = calculate_entropy(img)
%
% Calculates Shannon entropy for an image (grayscale).
%

    % Ensure uint8 type for histcounts
    img = uint8(img); 

    % 1. Count occurrences
    counts = histcounts(img, -0.5:1:255.5);
    
    % 2. Calculate probabilities
    total_pixels = numel(img);
    probabilities = counts / total_pixels;
    
    % 3. Calculate entropy
    probabilities = probabilities(probabilities > 0); % Remove p=0
    H = -sum(probabilities .* log2(probabilities));

end
