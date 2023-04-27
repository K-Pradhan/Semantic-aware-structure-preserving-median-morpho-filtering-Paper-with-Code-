clear;clc;  close all;

itr = 3; w_sk = 5;

[InputImage] = imread('Input/fish.jpg');
figure;
imshow(InputImage);

% disp('Size of the Image');
[m,n,c] = size(InputImage);

nhood_gr = ones(3);
nhood = ones(3); 
nhood2 = ones(3); 

FilteredImage = InputImage; 


for fltr_rpt= 1:itr
  
if (fltr_rpt>1)
if(c>1)
GrayImg = rgb2gray(FilteredImage);
else
GrayImg = FilteredImage;  
end
GrayImg = im2double(GrayImg);
else
  if(c>1)
    GrayImg = rgb2gray(FilteredImage);
  else
    GrayImg = FilteredImage;  
  end
GrayImg = im2double(GrayImg);  
% GrayImg = medfilt2(GrayImg, [3 3]);
GrayImg1 = imopen(GrayImg,nhood2) ; 
GrayImg2 = imclose(GrayImg,nhood2); 
GrayImg = (GrayImg1+GrayImg2)/2;

end
eps = 0.001;


% for grd=1:1
w=1;
[X,Y] = meshgrid(-w:w,-w:w);
for pxl = 1:m*n
[I,J] = ind2sub([m,n],pxl);

X1 = X+I; Y1 =Y+J;
XY = (X1>0)&(X1<m) & (Y1>0)&(Y1<n);
X1 = X1(XY); Y1 = Y1(XY);

ImdNbr = sub2ind([m,n],X1,Y1);
Pxl_Range(pxl) = max(GrayImg(ImdNbr)) - min(GrayImg(ImdNbr)) + eps; 
Morph_GradImg(I,J) = Pxl_Range(pxl);
end

Morph_GradImg = mat2gray(Morph_GradImg);
 
% end

eps = 0.001;
lmu = mean(Pxl_Range+eps);
lsd = std(Pxl_Range+eps);
nmu = log(lmu/sqrt(1+lsd^2/lmu^2));
nsd = sqrt(log(1+lsd^2/lmu^2));
Morph_GradImg = mat2gray(Morph_GradImg); 
  
BW = imbinarize(Morph_GradImg+eps,exp(nmu+nsd)); 

w = w_sk; 
[X,Y] = meshgrid(-w:w,-w:w); W =-w:w;
parfor pxl = 1:m*n   

Nbr1 = pxl + W*m; Nbr1 = repmat(Nbr1,2*w+1,1); Nbr1 = Nbr1 + Y;
Nbr1 = Nbr1(Nbr1>0 & Nbr1<=m*n);

Range_Vals1 = Pxl_Range(Nbr1);

mn = mean(Range_Vals1); md = median(Range_Vals1);  sd = std(Range_Vals1);
skwns1 =  3 *(mn - md) /(sd); 

if ( (BW(pxl) ==  1) && skwns1 > 0.0) 
    BW(pxl) =  1; 
else
  BW(pxl) =  0;
end
end

[FilteredImage] = AdaptiveMedianFilter(FilteredImage,GrayImg,BW) ;
if(c>1)
FilteredImage1 = imopen(FilteredImage,nhood); 
FilteredImage2 = imclose(FilteredImage,nhood); 
FilteredImage = (FilteredImage1+FilteredImage2)/2;
else
FilteredImage1 = imopen(FilteredImage,nhood_gr); 
FilteredImage2 = imclose(FilteredImage,nhood_gr);
FilteredImage = (FilteredImage1+FilteredImage2)/2;
end

end
figure;
imshow(mat2gray(FilteredImage));

