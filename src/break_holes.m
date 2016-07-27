% Disclaimer:  IMPORTANT:  This software was developed at the National Institute of Standards and Technology by employees of the Federal Government in the course of their official duties. Pursuant to title 17 Section 105 of the United States Code this software is not subject to copyright protection and is in the public domain. This is an experimental system. NIST assumes no responsibility whatsoever for its use by other parties, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. We would appreciate acknowledgment if the software is used. This software can be redistributed and/or modified freely provided that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.



% I1 is the double precision image holding the pixel wights
% BW is the logical mask containing the holes to be broken
function BW = break_holes(I1, S, BW) 

I1 = double(I1);
filled_BW = imfill(BW, 'holes');
filled_BW(~S) = 0;
xor_img = xor(BW, filled_BW);
xor_img = bwareaopen(xor_img,3);

if max(xor_img(:)) == 0, return, end
    
% if there are holes to fill
% this is a loop to deal with nested holes
I1 = I1.^2;
I1(xor_img) = 0;
iter = 1;
while iter<100
    stats = regionprops(xor_img, 'Centroid');
    cent = zeros(numel(stats),2);
    for i = 1:numel(stats)
        cent(i,:) = round(stats(i).Centroid);
    end

    gd2 = graydist(I1, ~filled_BW);
    gd1 = graydist(I1, cent(:,1),cent(:,2));
    break_pixels = imregionalmin(gd1 + gd2)>0;
    break_pixels = bwmorph(break_pixels, 'diag');

    BW(break_pixels) = 0;
    
    filled_BW = imfill(BW, 'holes');
    filled_BW(~S) = 0;
    xor_img = xor(BW, filled_BW);
    xor_img = bwareaopen(xor_img,3);

    if max(xor_img(:)) == 0, break, end
    iter = iter + 1;
end




