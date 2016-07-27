% Disclaimer:  IMPORTANT:  This software was developed at the National Institute of Standards and Technology by employees of the Federal Government in the course of their official duties. Pursuant to title 17 Section 105 of the United States Code this software is not subject to copyright protection and is in the public domain. This is an experimental system. NIST assumes no responsibility whatsoever for its use by other parties, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. We would appreciate acknowledgment if the software is used. This software can be redistributed and/or modified freely provided that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.



function T = percentile_computation(A,p)

assert(all(p >= 0 & p <= 1), 'Percentiles must be between zero and one inclusive.');

B = A(~isnan(A)); % ignore NaNs
if isempty(B)
    T = NaN(size(p));
    return;
end

B = sort(B(:));

indx = round(p*length(B)+1);

indx(indx<1) = 1;
indx(indx>length(B)) = length(B);

T = B(indx);

T = reshape(T, size(p));

end

