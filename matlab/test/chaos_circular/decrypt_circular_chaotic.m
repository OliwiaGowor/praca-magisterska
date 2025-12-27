function decrypted_img = decrypt_circular_chaotic(encrypted_img, keys_struct)
    % Unpack keys and parameters (as input arguments)
    encrypted = double(encrypted_img);
    keys = keys_struct.keys;
    rows = keys_struct.rows;
    cols = keys_struct.cols;
    DD = keys_struct.DD;
    
    x=keys(1); y=keys(2); a=keys(3); b=keys(4);
    N=keys(5); A=keys(6); B=keys(7); C=keys(8); D=keys(9);
    
    % Initialize empty matrices
    posrow=NaN(1,rows);
    poscol=NaN(1,cols*8);
    prbgmatvec=NaN(1,rows*cols);
    x1=NaN(1,8);
    y1=NaN(1,8);
    
    % Key iterations 
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
    
    % Calculate indexing y
    for i=1:8
        y=mod(y+ a + b*smht(y/N,A,B,C,D),N);
        y1(i)=y;
    end
    [~,indexY]=sort(y1);
    
    % STEP 6 & 7 - Reverse Substitution XOR
    for i=1:rows*cols
        y=mod(y+ a + b*smht(y/N,A,B,C,D),N);
        stream=mod(floor(DD*y),256);
        prbgmatvec(i)=stream;
    end
    y=mod(y+ a + b*smht(y/N,A,B,C,D),N);
    shiftval=floor(rows*cols*y);
    prbgmatvec=circshift(prbgmatvec,shiftval);
    prbgmat=reshape(prbgmatvec,rows,cols);
    
    plainbitsshuffled=bitxor(encrypted,prbgmat);
    
    % Reverse STEP 5
    plainbits=[bitget(plainbitsshuffled, find(indexY==1)), bitget(plainbitsshuffled, find(indexY==2)),...
     bitget(plainbitsshuffled, find(indexY==3)), bitget(plainbitsshuffled, find(indexY==4)),...
     bitget(plainbitsshuffled, find(indexY==5)),bitget(plainbitsshuffled, find(indexY==6)),...
     bitget(plainbitsshuffled, find(indexY==7)),bitget(plainbitsshuffled, find(indexY==8))];
    
    % Calculate indexing x (STEP 2)
    for i=1:8
        x=mod(x+ a + b*smht(x/N,A,B,C,D),N);
        x1(i)=x;
    end
    [~,indexX]=sort(x1);
    
    % Calculate shifts (STEP 3 & 4)
    for j=1:rows
        x=mod(x+ a + b*smht(x/N,A,B,C,D),N);
        posrow(j)=floor(8*cols*x);
    end
    x=mod(x+ a + b*smht(x/N,A,B,C,D),N);
    pos=floor(rows*x);
    posrow=circshift(posrow,pos);
    
    for j=1:cols*8
        x=mod(x+ a + b*smht(x/N,A,B,C,D),N);
        poscol(j)=floor(rows*x);
    end
    x=mod(x+ a + b*smht(x/N,A,B,C,D),N);
    pos=floor(cols*8*x);
    poscol=circshift(poscol,pos);
    
    % Reverse column shifts
    for j=1:cols*8
        plainbits(:,j)=circshift(plainbits(:,j),-poscol(j));
    end
    
    % Reverse row shifts
    for j=1:rows
        plainbits(j,:)=circshift(plainbits(j,:),-posrow(j));
    end
    
    % Reconstruct image (STEP 2 - final)
    decrypted=plainbits(:,(find(indexX==1)-1)*cols+1:find(indexX==1)*cols)+...
        plainbits(:,(find(indexX==2)-1)*cols+1:find(indexX==2)*cols)*2+...
        plainbits(:,(find(indexX==3)-1)*cols+1:find(indexX==3)*cols)*4+...
        plainbits(:,(find(indexX==4)-1)*cols+1:find(indexX==4)*cols)*8+...
        plainbits(:,(find(indexX==5)-1)*cols+1:find(indexX==5)*cols)*16+...
        plainbits(:,(find(indexX==6)-1)*cols+1:find(indexX==6)*cols)*32+...
        plainbits(:,(find(indexX==7)-1)*cols+1:find(indexX==7)*cols)*64+...
        plainbits(:,(find(indexX==8)-1)*cols+1:find(indexX==8)*cols)*128;
 
    decrypted_img = uint8(decrypted);
    
    end
