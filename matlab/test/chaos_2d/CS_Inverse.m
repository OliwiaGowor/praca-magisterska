function P = CS_Inverse(Y, cskey)
    % Orthogonal inverse CS (practical reconstruction)
    
    Phi_r = cskey.Phi_r;
    Phi_c = cskey.Phi_c;
    Psi_r = cskey.Psi_r;
    Psi_c = cskey.Psi_c;
    
    % --- pseudo-inverse reconstruction
    S_hat = Phi_r' * Y * Phi_c;
    
    P = Psi_r' * S_hat * Psi_c';
    P = uint8(mod(round(P),256));
    end
