function [encrypted_img, keys_struct] = encrypt_circular_chaotic(plaintextim)
%
% == POCZĄTEK ORYGINALNEGO KODU SZYFROWANIA ==
% Logika jest identyczna jak w Twoim skrypcie.
% Zmieniono tylko wejście (teraz argument funkcji) i wyjście (zwracane zmienne).
%

% plaintextim = imread('cameraman.tif'); % Usunięte - teraz to argument
plaintext=double(plaintextim);
[rows,cols]=size(plaintext);
plainbits=zeros(rows,cols*8);
posrow=NaN(1,rows);
poscol=NaN(1,cols*8);
indexX=NaN(1,8);
indexY=NaN(1,8);
prbgmat=NaN(rows,cols);
prbgmatvec=NaN(1,rows*cols);
plainbitsshuffled=NaN(rows,cols);
% encrypted=NaN(rows,cols); % Niepotrzebne, alokowane przez bitxor

DD=10^10;
keys=NaN(1,9);
x1=NaN(1,8);
y1=NaN(1,8);
ten16=10^16;

% STEP 1 - Key initialization
hash=DataHash(plaintext,'SHA-256');
h1=hex2dec(hash(1:13));
h2=hex2dec(hash(14:26));
h3=hex2dec(hash(27:39));
h4=hex2dec(hash(40:52));
h5=hex2dec(hash(53:end));
x=(h1+h2+h5)/ten16;
y=(h3+h4+h5)/ten16;

a=5; b=5; N=1; A=0.84; B=0.75; C=1; D=1;

% Zapisz klucze do struktury, która zostanie zwrócona
keys_vec=[x,y,a,b,N,A,B,C,D];
keys_struct.keys = keys_vec;
keys_struct.rows = rows;
keys_struct.cols = cols;
keys_struct.DD = DD;

for i=1:49
    x=rem(x+ a + b*smht(x/N,A,B,C,D),N);
    y=rem(y+ a + b*smht(y/N,A,B,C,D),N);
end
xold=x;
yold=y;
x=rem(x+ a + b*smht(x/N,A,B,C,D),N);
y=rem(y+ a + b*smht(y/N,A,B,C,D),N);
x=rem(x+yold,1);
y=rem(xold+y,1);

% STEP 2 - Arrange (chaotically) the bit planes
for i=1:8
    x=rem(x+ a + b*smht(x/N,A,B,C,D),N);
    x1(i)=x;
end
[~,indexX]=sort(x1);
plainbits=logical([bitget(plaintext, indexX(1)), bitget(plaintext, indexX(2)),...
 bitget(plaintext, indexX(3)), bitget(plaintext, indexX(4)),...
 bitget(plaintext, indexX(5)),bitget(plaintext, indexX(6)),...
 bitget(plaintext, indexX(7)),bitget(plaintext, indexX(8))]);

% STEP 3 & 4 - Perform circshift on rows and cols
for j=1:rows
    x=rem(x+ a + b*smht(x/N,A,B,C,D),N);
    posrow(j)=floor(8*cols*x);
end
x=rem(x+ a + b*smht(x/N,A,B,C,D),N);
pos=floor(rows*x);
posrow=circshift(posrow,pos);

for j=1:cols*8
    x=rem(x+ a + b*smht(x/N,A,B,C,D),N);
    poscol(j)=floor(rows*x);
end
x=rem(x+ a + b*smht(x/N,A,B,C,D),N);
pos=floor(cols*8*x);
poscol=circshift(poscol,pos);

for j=1:rows
    plainbits(j,:)=circshift(plainbits(j,:),posrow(j));
end
for j=1:cols*8
    plainbits(:,j)=circshift(plainbits(:,j),poscol(j));
end

% STEP 5 - Append (chaotically) the bit planes
for i=1:8
    y=rem(y+ a + b*smht(y/N,A,B,C,D),N);
    y1(i)=y;
end
[~,indexY]=sort(y1);
plainbitsshuffled=...
    plainbits(:,(indexY(1)-1)*cols+1:indexY(1)*cols)+...
    plainbits(:,(indexY(2)-1)*cols+1:indexY(2)*cols)*2+...
    plainbits(:,(indexY(3)-1)*cols+1:indexY(3)*cols)*4+...
    plainbits(:,(indexY(4)-1)*cols+1:indexY(4)*cols)*8+...
    plainbits(:,(indexY(5)-1)*cols+1:indexY(5)*cols)*16+...
    plainbits(:,(indexY(6)-1)*cols+1:indexY(6)*cols)*32+...
    plainbits(:,(indexY(7)-1)*cols+1:indexY(7)*cols)*64+...
    plainbits(:,(indexY(8)-1)*cols+1:indexY(8)*cols)*128;

% STEP 6 - Compute the substitution matrix
for i=1:rows*cols
    y=rem(y+ a + b*smht(y/N,A,B,C,D),N);
    prbgmatvec(i)=rem(floor(DD*y),256);
end
y=rem(y+ a + b*smht(y/N,A,B,C,D),N);
shiftval=floor(rows*cols*y);
prbgmatvec=circshift(prbgmatvec,shiftval);
prbgmat=reshape(prbgmatvec,rows,cols);

% STEP 6 - Substitution through XOR
encrypted=bitxor(plainbitsshuffled,prbgmat);

%
% == KONIEC ORYGINALNEGO KODU SZYFROWANIA ==
%

% Zwróć wynik jako uint8 (Twój skrypt pomijał ten krok, ale imshow robi to niejawnie)
encrypted_img = uint8(encrypted);

end
