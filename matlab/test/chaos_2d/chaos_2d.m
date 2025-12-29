%% chaos_2d.m
clc; clear; close all;

%% 1. Load image
P = imread('cameraman.tif');
if size(P,3)==3, P = rgb2gray(P); end
P = uint8(P);
[M,N] = size(P);

%% 2. Parameters (Entropy 2024)
params.a  = 20;
params.b  = 30;
params.x0 = 0.123987;
params.y0 = 0.987321;
params.Pre = uint8(50);

%% 3. Encrypt
[C, keys] = EncryptImage_Entropy(P, params);

%% 4. Decrypt
P_rec = DecryptImage_Entropy(C, keys, params);

%% 5. Verify
disp(['Reconstruction error = ', ...
    num2str(sum(abs(double(P(:))-double(P_rec(:)))))]);

figure;
subplot(1,3,1); imshow(P);      title('Plain');
subplot(1,3,2); imshow(C);      title('Cipher');
subplot(1,3,3); imshow(P_rec);  title('Recovered');

%% =========================================================
%  Hyperchaotic map (Eq. 1)
% =========================================================
function [x,y] = HyperchaoticMap(x,y,a,b)
    x = sin(a*pi*x + b*y)^2;
    if y==0
        y = 0;
    else
        y = cos(b*pi/y + a*x)^2;
    end
end

%% =========================================================
%  Algorithm 3 — Encryption (Permutation + Double Diffusion)
% =========================================================
function [C, key] = EncryptImage_Entropy(P, params)

[M,N] = size(P);
L = M*N;

% --- generate hyperchaotic sequences
x = params.x0; y = params.y0;
X = zeros(1,L); Y = zeros(1,L);
for i = 1:L
    [x,y] = HyperchaoticMap(x,y,params.a,params.b);
    X(i) = x; 
    Y(i) = y;
end

% --- permutation indexes (position scrambling)
perm = mod(floor(X*1e14),L) + 1;
[~,perm] = unique(perm,'stable');
inv_perm(perm) = 1:L;

% --- keystream for value diffusion
K = uint8(mod(floor(Y*1e14),256));

% --- apply permutation
CR = 0.5; % compression ratio 
[Y, cskey] = CS_Forward(P, X, Y, CR);

Pcs = uint8(mod(round(Y),256));
Pvec = Pcs(:);
Pp = Pvec(perm);


% ---------- Forward diffusion (Eq. 14) ----------
C1  = zeros(L,1,'uint8');
Pre = params.Pre;
for i = 1:L
    C1(i) = bitxor(uint8(mod(double(Pp(i))+double(Pre),256)), K(i));
    Pre   = C1(i);
end

% ---------- Backward diffusion (Eq. 15) ----------
C2  = zeros(L,1,'uint8');
Pre = params.Pre;
for i = L:-1:1
    C2(i) = bitxor(uint8(mod(double(C1(i))+double(Pre),256)), K(i));
    Pre   = C2(i);
end

C = reshape(C2,M,N);

key.cskey = cskey;
key.CR = CR;
key.perm     = perm;
key.inv_perm = inv_perm;
key.K        = K;
end

%% =========================================================
%  Algorithm 4 — Decryption
% =========================================================
function P = DecryptImage_Entropy(C, key, params)

[M,N] = size(C);
L = M*N;

Cvec = C(:);
K    = key.K;
perm = key.perm;
invp = key.inv_perm;

% ---------- inverse backward diffusion ----------
C1  = zeros(L,1,'uint8');
Pre = params.Pre;
for i = L:-1:1
    C1(i) = uint8(mod(double(bitxor(Cvec(i),K(i))) - double(Pre),256));
    Pre   = Cvec(i);
end

% ---------- inverse forward diffusion ----------
Pp  = zeros(L,1,'uint8');
Pre = params.Pre;
for i = 1:L
    Pp(i) = uint8(mod(double(bitxor(C1(i),K(i))) - double(Pre),256));
    Pre   = C1(i);
end

% ---------- inverse permutation ----------
Yrec = reshape(Pp(invp), size(key.cskey.Phi_r,1), []);
P = CS_Inverse(double(Yrec), key.cskey);

end
