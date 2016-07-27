% Disclaimer:  IMPORTANT:  This software was developed at the National Institute of Standards and Technology by employees of the Federal Government in the course of their official duties. Pursuant to title 17 Section 105 of the United States Code this software is not subject to copyright protection and is in the public domain. This is an experimental system. NIST assumes no responsibility whatsoever for its use by other parties, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. We would appreciate acknowledgment if the software is used. This software can be redistributed and/or modified freely provided that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.



function BW = generate_border_mask(image, img_filter, prctile_threshold, threshold_direction, border_break_holes_flag, border_thin_mask_flag, foreground_mask)

if ~exist('foreground_mask','var')
	foreground_mask = ones(size(image));
end
if ~exist('img_filter','var')
    img_filter = 'none';
end
if ~exist('border_break_holes_flag','var')
    border_break_holes_flag = false;
end
if ~exist('border_thin_mask_flag','var')
    border_thin_mask_flag = false;
end

% input checking
assert(exist('image','var')>0, 'Image to generate seed from missing');
assert(exist('prctile_threshold','var')>0, 'Prctile threshold used to generate seed from missing');

assert(prctile_threshold >= 0 && prctile_threshold <= 100);
prctile_threshold = prctile_threshold./100;
assert(size(image,1) == size(foreground_mask,1) && size(image,2) == size(foreground_mask,2), 'border mask image must be the same size as the foreground mask');

image = double(image);

img_filter = lower(regexprep(img_filter, '\W', ''));
switch img_filter
	case 'gradient'
		image = imgradient(image, 'Sobel');
	case 'entropy'
		image = entropyfilt(image,true(5,5));
	case 'std'
		image = stdfilt(image,true(5,5));
    otherwise
		% do nothing
end

% threshold the image
P = percentile_computation(image(:), prctile_threshold);
switch threshold_direction
    case '>'
        BW = image > P(1);
    case '<'
        BW = image < P(1);
    case '>='
        BW = image >= P(1);
    case '<='
        BW = image <= P(1);
    otherwise
        error('invalid threshold operator');
end

% remove any pixels from the seed mask that are not part of the overall foreground
BW = BW & foreground_mask;

% morphological cleanup
BW = imdilate(BW,strel('disk',1));
BW = bwmorph(BW, 'bridge');
BW = bwmorph(BW, 'thin',1);
BW = bwmorph(BW, 'diag');

BW = fill_holes(BW, 10); % fill holes smaller han 10 pixels

if border_thin_mask_flag
    BW = bwmorph(BW, 'thin',inf);
end
BW = bwmorph(BW, 'diag');

if border_break_holes_flag
    % break the holes so they are no longer fillable
    BW = break_holes(image, foreground_mask, BW);
end

end