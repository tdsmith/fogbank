% Disclaimer:  IMPORTANT:  This software was developed at the National Institute of Standards and Technology by employees of the Federal Government in the course of their official duties. Pursuant to title 17 Section 105 of the United States Code this software is not subject to copyright protection and is in the public domain. This is an experimental system. NIST assumes no responsibility whatsoever for its use by other parties, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. We would appreciate acknowledgment if the software is used. This software can be redistributed and/or modified freely provided that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.



function BW = filter_by_circularity(BW, circularity_threshold) 

stats = regionprops(BW, 'Area','Perimeter','PixelIdxList');
for i = 1:numel(stats)
    p = stats(i).Perimeter;
    a = stats(i).Area;
    circ = ((4*pi)*a)/(p^2);
    if circ < circularity_threshold
        % remove seed
        BW(stats(i).PixelIdxList) = 0;
    end
end

end