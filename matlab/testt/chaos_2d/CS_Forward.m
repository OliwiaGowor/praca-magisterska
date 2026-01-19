function [Y, cskey] = CS_Forward(P, chaosX, chaosY, CR)
    % 2D compressed sensing (Entropy-style)
    % CR – compression ratio (0.25–0.5)
    
    [M,N] = size(P);
    
    % --- sparsifying basis (DCT)
    Psi_r = dctmtx(M);
    Psi_c = dctmtx(N);
    S = Psi_r * double(P) * Psi_c';
    
    % --- measurement matrix (chaos)
    Mr = round(CR*M);
    Mc = round(CR*N);
    
    Pr = mod(floor(chaosX(1:Mr*M)*1e14),M)+1;
    Pc = mod(floor(chaosY(1:Mc*N)*1e14),N)+1;
    
    Phi_r = eye(M);
    Phi_c = eye(N);
    
    Phi_r = Phi_r(Pr,:);
    Phi_c = Phi_c(Pc,:);
    
    % --- CS measurement (Eq. 8–10 in paper)
    Y = Phi_r * S * Phi_c';
    
    cskey.Phi_r = Phi_r;
    cskey.Phi_c = Phi_c;
    cskey.Psi_r = Psi_r;
    cskey.Psi_c = Psi_c;
    end
