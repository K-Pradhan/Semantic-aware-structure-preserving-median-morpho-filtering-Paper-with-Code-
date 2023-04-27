function [FilteredImg] = AdaptiveMedianFilter(InputImage,GrayImg,BW)

[m , n, c] = size(InputImage); 
FilteredImg = InputImage; output = zeros(m*n,c); 
if(c>1)
Inputbnd1 = InputImage(:,:,1);Inputbnd2 = InputImage(:,:,2);Inputbnd3 = InputImage(:,:,3);
else
Inputbnd1 = InputImage;
end

w1=1;
W_imdt = -w1:w1;
[~,Y_imdt] = meshgrid(-w1:w1,-w1:w1);
w4=1;
W_enbr = -w4:w4;
[~,Y_enbr] = meshgrid(-w4:w4,-w4:w4);
w2 = 2; 
W_edge = -w2:w2;
[~,Y_edge] = meshgrid(-w2:w2,-w2:w2);
txt_end = 1;

TxtEnd = 25;
txtrwndw = 1;
w3=txtrwndw; W3 = -w3:w3;
[~,Y_texture] = meshgrid(-w3:w3,-w3:w3); 
TxtNbr = cell([1 TxtEnd]); TxtNbr{1} = Y_texture; TxtWnd = cell([1 TxtEnd]); TxtWnd{1} = W3; TxtWndBound = cell([1 TxtEnd]); TxtWndBound{1} = W3;
for txtrwndw=2:TxtEnd
           w3=txtrwndw;
           Y_texture = [Y_texture(:,1) Y_texture Y_texture(:,end)];
           Y_texture = [-w3*ones(1,2*w3+1); Y_texture; w3*ones(1,2*w3+1)];
          W3 = [-w3 W3 w3]; 
          TxtNbr{txtrwndw} = Y_texture; TxtWnd{txtrwndw} = W3; 
end

parfor pxl = 1:m*n

Nbr_imdt = pxl + W_imdt*m; 
Nbr_imdt = repmat(Nbr_imdt,2*w1+1,1); 
Nbr_imdt = Nbr_imdt + Y_imdt; 
Nbr_imdt = Nbr_imdt(Nbr_imdt>0 & Nbr_imdt<=m*n);

TxtrNbr =  Nbr_imdt; 

Nbr_enbr = pxl + W_enbr*m; 
Nbr_enbr = repmat(Nbr_enbr,2*w4+1,1); 
Nbr_enbr = Nbr_enbr + Y_enbr; 
Nbr_enbr = Nbr_enbr(Nbr_enbr>0 & Nbr_enbr<=m*n);

Nbr = pxl + W_edge*m; 
Nbr = repmat(Nbr,2*w2+1,1); 
Nbr = Nbr + Y_edge; 
Nbr = Nbr(Nbr>0 & Nbr<=m*n);

Filter_Wghts1 = []; Filter_Wghts2 = []; Filter_Wghts3 = [];
%###################################
if (  BW(pxl)==1)   
    midR = (max(GrayImg(Nbr)) + min(GrayImg(Nbr)))/2; 
     
     if ((median(GrayImg(Nbr_enbr)) <= midR)) 
         Nbr_New = Nbr(GrayImg(Nbr) <= midR);
         if(c>1)
         Filter_Wghts1 =  median(Inputbnd1(Nbr_New)); 
         Filter_Wghts2 =  median(Inputbnd2(Nbr_New));
         Filter_Wghts3 =  median(Inputbnd3(Nbr_New));
         else
         Filter_Wghts1 =  median(Inputbnd1(Nbr_New));
         end
     elseif((median(GrayImg(Nbr_enbr)) > midR))
         Nbr_New = Nbr(GrayImg(Nbr) > midR); 
          if(c>1)
         Filter_Wghts1 =  median(Inputbnd1(Nbr_New)); 
         Filter_Wghts2 =  median(Inputbnd2(Nbr_New));
         Filter_Wghts3 =  median(Inputbnd3(Nbr_New));
         else
         Filter_Wghts1 =  median(Inputbnd1(Nbr_New));
         end
     end
else
    txtrwndw = txt_end; 
    w3=txtrwndw; W3 = -w3:w3;
    [~,Y_texture] = meshgrid(-w3:w3,-w3:w3);
    while (nnz(BW(TxtrNbr))<= 2 && w3<TxtEnd) 
          txtrwndw = txtrwndw + 1;
          w3=txtrwndw;
          Y_texture = TxtNbr{txtrwndw};       
          W3 = TxtWnd{txtrwndw}; 
          TxtrNbr = pxl + W3*m; TxtrNbr = repmat(TxtrNbr,2*w3+1,1); TxtrNbr = TxtrNbr + Y_texture; % ImdNbr = [ImdNbr - w ImdNbr + w];
          TxtrNbr = TxtrNbr(TxtrNbr>0 & TxtrNbr<=m*n);
    end

    TxtrNbr_final = pxl + W3*m; TxtrNbr_final = repmat(TxtrNbr_final,2*w3+1,1); TxtrNbr_final = TxtrNbr_final + Y_texture; % ImdNbr = [ImdNbr - w ImdNbr + w];
    TxtrNbr_final = TxtrNbr_final(TxtrNbr_final>0 & TxtrNbr_final<=m*n);

    if(c>1)
       Filter_Wghts1 =  median(Inputbnd1(TxtrNbr_final));
       Filter_Wghts2 =  median(Inputbnd2(TxtrNbr_final));
       Filter_Wghts3 =  median(Inputbnd3(TxtrNbr_final));
    else
       Filter_Wghts1 =  median(Inputbnd1(TxtrNbr_final));
    end

end
    if(c>1)
        output(pxl,:) = [Filter_Wghts1; Filter_Wghts2; Filter_Wghts3];
    else
    FilteredImg(pxl) =  Filter_Wghts1;
    end
end

 if(c>1) 
     FilteredImg = reshape(output,[m n c]);
 end

 end
