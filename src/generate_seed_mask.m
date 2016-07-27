% Disclaimer:  IMPORTANT:  This software was developed at the National Institute of Standards and Technology by employees of the Federal Government in the course of their official duties. Pursuant to title 17 Section 105 of the United States Code this software is not subject to copyright protection and is in the public domain. This is an experimental system. NIST assumes no responsibility whatsoever for its use by other parties, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. We would appreciate acknowledgment if the software is used. This software can be redistributed and/or modified freely provided that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.

% image the image to use in generating the seed mask
% img_filter filter operation to apply to the image
% 	(None, Gradient, Entropy, Std)
% prctile_threshold threshold defining foreground
% threshold_direction wether above or below the threshold is considered foreground
% 	(<,>=)
% min_obj_size minimum area in pixels of an object
% max_obj_size maximum area in pixels of an object
% circularity_threshold minimum circularity objects must have 
% cluster_distance
% foreground_mask (optional)
% border_mask (optional)


function BW = generate_seed_mask(image, img_filter, prctile_thresholdL, threshold_operatorL, prctile_thresholdR, threshold_operatorR, min_obj_size, max_obj_size, circularity_threshold, cluster_distance, foreground_mask, border_mask)

if ~exist('foreground_mask','var')
	foreground_mask = ones(size(image));
end
if ~exist('img_filter','var')
    img_filter = 'none';
end


% input checking
assert(exist('image','var')>0, 'Image to generate seed from missing');
assert(exist('prctile_thresholdL','var')>0, 'Prctile threshold used to generate seed from missing');
assert(exist('prctile_thresholdR','var')>0, 'Prctile threshold used to generate seed from missing');
assert(exist('min_obj_size','var')>0, 'Minimum object size missing');
assert(exist('max_obj_size','var')>0, 'Maximum object size missing');

assert(prctile_thresholdL >= 0 && prctile_thresholdL <= 100);
assert(prctile_thresholdR >= 0 && prctile_thresholdR <= 100);
prctile_thresholdL = prctile_thresholdL./100;
prctile_thresholdR = prctile_thresholdR./100;
assert(size(image,1) == size(foreground_mask,1) && size(image,2) == size(foreground_mask,2), 'Seed mask image must be the same size as the foreground mask');


image = double(image);
% modify the foreground mask to account for the border mask
if exist('border_mask', 'var') && ~isempty(border_mask)
	foreground_mask = foreground_mask & ~border_mask;
end

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
% P = percentile_computation(image(:), [prctile_thresholdL,prctile_thresholdR]);
P = percentile_computation(image(foreground_mask), [prctile_thresholdL,prctile_thresholdR]);
switch threshold_operatorL
    case '>'
        BW_L = image > P(1);
    case '<'
        BW_L = image < P(1);
    case '>='
        BW_L = image >= P(1);
    case '<='
        BW_L = image <= P(1);
    otherwise 
        error('Invalid threshold operator');
end
switch threshold_operatorR
    case '>'
        BW_R = image > P(2);
    case '<'
        BW_R = image < P(2);
    case '>='
        BW_R = image >= P(2);
    case '<='
        BW_R = image <= P(2);
    otherwise 
        error('Invalid threshold operator');
end

BW = BW_L & BW_R;

% remove any pixels from the seed mask that are not part of the overall foreground
BW = BW & foreground_mask;

% fill in the holes that are smaller than twice the minimum object size
BW = fill_holes(BW, min_obj_size*2);
% remove objects that are smaller than the minimum object size
BW = bwareaopen(BW, min_obj_size);

% remove the objects in BW above seed_max_object_size
CC = bwconncomp(BW);
for ii = 1:CC.NumObjects
	if numel(CC.PixelIdxList{ii}) > max_obj_size
		BW(CC.PixelIdxList{ii}) = 0; % remove the object by setting its pixels to zero
	end
end


if exist('circularity_threshold','var')
    % circularity check, only keeps cells below the circularity threshold
    BW = filter_by_circularity(BW, circularity_threshold);
end

if exist('cluster_distance','var')
    % cluster seed objects
    BW = cluster_objects(BW, cluster_distance, foreground_mask);
end


end