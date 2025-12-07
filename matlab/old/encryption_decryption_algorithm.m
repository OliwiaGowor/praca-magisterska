%{
Lazaros Moysis Youtube channel: https://www.youtube.com/@lazarosmoysis5095
RG: https://www.researchgate.net/profile/Lazaros-Moysis

The code implements the chaotic encryption algorithm proposed in the
following work:

Moysis, L., Lawnik, M., Alexan, W., Goudos, S. K., Baptista, M. S., &
Fragulis, G. F. (2025). Exploiting Circular Shifts for Efficient Chaotic
Image Encryption. IEEE Access.
https://ieeexplore.ieee.org/document/11009169 

Please cite this work if you use the code below.

The code is broken in sections, devided by the %% symbol. Run each section
separately using ctr+enter, or by clicking the 'run section' button.

The work is open access, so you can find information about the algorithm
and its steps in the pdf.

The algorithm uses the DataHash function by Jan, from:
https://uk.mathworks.com/matlabcentral/fileexchange/31272-datahash 

%}

%% Encryption Algorithm
clear
clc

% load image, and transform into grayscale if required

plaintextim = imread('cameraman.tif');
% plaintextim=rgb2gray(plaintextim);

% tansform into double
plaintext=double(plaintextim);

% extract row/column sizes
[rows,cols]=size(plaintext);

% initialize empty matrices, to save on execution time
plainbits=zeros(rows,cols*8);
posrow=NaN(1,rows);
poscol=NaN(1,cols*8);
indexX=NaN(1,8);
indexY=NaN(1,8);
prbgmat=NaN(rows,cols);
prbgmatvec=NaN(1,rows*cols);
plainbitsshuffled=NaN(rows,cols);
encrypted=NaN(rows,cols);

% Parameters for the PRNG to be used in the XOR stage
DD=10^10;
keys=NaN(1,9);
x1=NaN(1,8);
y1=NaN(1,8);
ten16=10^16;

% STEP 1 - Key initialization
% define initial conditions for the process
% these are the keys that are required for decryption 
hash=DataHash(plaintext,'SHA-256');
h1=hex2dec(hash(1:13));
h2=hex2dec(hash(14:26));
h3=hex2dec(hash(27:39));
h4=hex2dec(hash(40:52));
h5=hex2dec(hash(53:end));
x=(h1+h2+h5)/ten16;
y=(h3+h4+h5)/ten16;

% the rest of the chaotic map keys
% here, the 2 maps use the same parameters, you can change that for a
% higher key space
a=5;
b=5;
N=1;
A=0.84;
B=0.75;
C=1;
D=1;

% save all the keys in a vector
keys=[x,y,a,b,N,A,B,C,D];

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

% STEP 2 - Arrange (chaotically) the bit planes in an Mx8N matrix
for i=1:8
    x=rem(x+ a + b*smht(x/N,A,B,C,D),N);
    % save these in a separate vector for simplicity, cause for x we don't
    % need to save the previous values
    x1(i)=x; 
end

[~,indexX]=sort(x1);

% obtain the matrix of all the bit planes
plainbits=logical([bitget(plaintext, indexX(1)), bitget(plaintext, indexX(2)),...
 bitget(plaintext, indexX(3)), bitget(plaintext, indexX(4)),...
 bitget(plaintext, indexX(5)),bitget(plaintext, indexX(6)),... 
 bitget(plaintext, indexX(7)),bitget(plaintext, indexX(8))]); 

% STEP 3 & 4 - Perform circshift on the rows and cols of Mx8N matrix
% compute the row shift values
for j=1:rows
    x=rem(x+ a + b*smht(x/N,A,B,C,D),N);
    posrow(j)=floor(8*cols*x);
end

% circshift the row shift values
x=rem(x+ a + b*smht(x/N,A,B,C,D),N);
pos=floor(rows*x);
posrow=circshift(posrow,pos);

% compute the column shift values
for j=1:cols*8
    x=rem(x+ a + b*smht(x/N,A,B,C,D),N);
    poscol(j)=floor(rows*x);
end

% circshift the column shift values
x=rem(x+ a + b*smht(x/N,A,B,C,D),N);
pos=floor(cols*8*x);
poscol=circshift(poscol,pos);

% Do the circular shifting of the rows
for j=1:rows
    plainbits(j,:)=circshift(plainbits(j,:),posrow(j));
end

% Do the circular shifting of the columns
for j=1:cols*8
    plainbits(:,j)=circshift(plainbits(:,j),poscol(j));
end


% STEP 5 - Append (chaotically) the bit planes into an MxN grayscale image
for i=1:8
    y=rem(y+ a + b*smht(y/N,A,B,C,D),N);
    % save these in a separate vector for simplicity
    y1(i)=y; 
end

[~,indexY]=sort(y1);

% obtain the grayscale image
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

% first, create the required bytestream of length M*N
for i=1:rows*cols
    y=rem(y+ a + b*smht(y/N,A,B,C,D),N);
    prbgmatvec(i)=rem(floor(DD*y),256);
end

% circshift the bytestream of length M*N
y=rem(y+ a + b*smht(y/N,A,B,C,D),N);
shiftval=floor(rows*cols*y);

prbgmatvec=circshift(prbgmatvec,shiftval);

% reshape the bytestream into a matrix M*N
prbgmat=reshape(prbgmatvec,rows,cols);

% STEP 6 - Substitution through XOR

% The resulting ciphertext:
encrypted=bitxor(plainbitsshuffled,prbgmat);


% display the obtained images
figure
imshow(uint8(plaintext))
title('Plaintext')
figure
imshow(uint8(encrypted))
title('ciphertext')

%% Decryption Algorithm

% To decrypt the ciphertext, we reverse the encryption steps in a logical
% order.

% first, clear any other variables except the keys and the ciphertext.
% because the receiver only has access to these values
clearvars -except rows cols keys encrypted plaintext DD

x=keys(1);
y=keys(2);
a=keys(3);
b=keys(4);
N=keys(5);
A=keys(6);
B=keys(7);
C=keys(8);
D=keys(9);


for i=1:49
    x=mod(x+ a + b*smht(x/N,A,B,C,D),N);
    y=rem(y+ a + b*smht(y/N,A,B,C,D),N);
end
xold=x;
yold=y;
x=mod(x+ a + b*smht(x/N,A,B,C,D),N);
y=rem(y+ a + b*smht(y/N,A,B,C,D),N);
x=rem(x+yold,1);
y=rem(xold+y,1);



% first compute the indexing
for i=1:8
    y=mod(y+ a + b*smht(y/N,A,B,C,D),N);
    % save these in a separate vector
    y1(i)=y; 
end
[~,indexY]=sort(y1);

% STEP 6 & 7 - Substitution through XOR

% first, create the required bytestream of length M*N
for i=1:rows*cols
    y=mod(y+ a + b*smht(y/N,A,B,C,D),N);
    stream=mod(floor(DD*y),256);
    prbgmatvec(i)=stream;
end

% circshift the bytestream of length M*N
y=mod(y+ a + b*smht(y/N,A,B,C,D),N);
shiftval=floor(rows*cols*y);
prbgmatvec=circshift(prbgmatvec,shiftval);

% reshape the bytestream into a matrix M*N
prbgmat=reshape(prbgmatvec,rows,cols);

% finally, reverse the xor
plainbitsshuffled=bitxor(encrypted,prbgmat);


% completion of STEP 5
plainbits=[bitget(plainbitsshuffled, find(indexY==1)), bitget(plainbitsshuffled, find(indexY==2)),...
 bitget(plainbitsshuffled, find(indexY==3)), bitget(plainbitsshuffled, find(indexY==4)),...
 bitget(plainbitsshuffled, find(indexY==5)),bitget(plainbitsshuffled, find(indexY==6)),... 
 bitget(plainbitsshuffled, find(indexY==7)),bitget(plainbitsshuffled, find(indexY==8))]; 


% STEP 2 (partial) - Arrange (chaotically) the bit planes in an Mx8N matrix
% first, compute the plane arrangement indexes
% the rearrangement will be performed later, after we reverse the row and
% column shifting.
for i=1:8
    x=mod(x+ a + b*smht(x/N,A,B,C,D),N);
    % save these in a separate vector for simplicity, cause for x we don't
    % need to save the previous values
    x1(i)=x; 
end
[~,indexX]=sort(x1);

% STEP 3 & 4 - Perform circshift on the rows and cols of Mx8N matrix

% save the row shift values
for j=1:rows
    x=mod(x+ a + b*smht(x/N,A,B,C,D),N);
    posrow(j)=floor(8*cols*x);
end

% circshift the row shift values
x=mod(x+ a + b*smht(x/N,A,B,C,D),N);
pos=floor(rows*x);
posrow=circshift(posrow,pos);

% save the column shift values
for j=1:cols*8
    x=mod(x+ a + b*smht(x/N,A,B,C,D),N);
    poscol(j)=floor(rows*x);
end

% circshift the column shift values
x=mod(x+ a + b*smht(x/N,A,B,C,D),N);
pos=floor(cols*8*x);
poscol=circshift(poscol,pos);

% now do the circular shifting of the columns
for j=1:cols*8
    plainbits(:,j)=circshift(plainbits(:,j),-poscol(j));
end

% now do the circular shifting of the rows
for j=1:rows
    plainbits(j,:)=circshift(plainbits(j,:),-posrow(j));
end

decrypted=plainbits(:,(find(indexX==1)-1)*cols+1:find(indexX==1)*cols)+...
    plainbits(:,(find(indexX==2)-1)*cols+1:find(indexX==2)*cols)*2+...
    plainbits(:,(find(indexX==3)-1)*cols+1:find(indexX==3)*cols)*4+...
    plainbits(:,(find(indexX==4)-1)*cols+1:find(indexX==4)*cols)*8+...
    plainbits(:,(find(indexX==5)-1)*cols+1:find(indexX==5)*cols)*16+...
    plainbits(:,(find(indexX==6)-1)*cols+1:find(indexX==6)*cols)*32+...
    plainbits(:,(find(indexX==7)-1)*cols+1:find(indexX==7)*cols)*64+...
    plainbits(:,(find(indexX==8)-1)*cols+1:find(indexX==8)*cols)*128;

figure
imshow(uint8(decrypted))
title('Decrypted image')

