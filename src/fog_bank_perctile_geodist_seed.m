% Disclaimer:  IMPORTANT:  This software was developed at the National Institute of Standards and Technology by employees of the Federal Government in the course of their official duties. Pursuant to title 17 Section 105 of the United States Code this software is not subject to copyright protection and is in the public domain. This is an experimental system. NIST assumes no responsibility whatsoever for its use by other parties, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. We would appreciate acknowledgment if the software is used. This software can be redistributed and/or modified freely provided that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.


% [seed_image, Highest_cell_number] = 
%     fog_bank_perctile_geodist2(
%         grayscale_image, 
%         foreground_mask,
%         mask_matrix, 
%         min_peak_size, 
%         min_object_size, 
%         perc_binning,             <optional>
% 

% colors_vector
function [seed_image, nb_peaks] = fog_bank_perctile_geodist_seed(grayscale_image, foreground_mask, mask_matrix, seed_image, min_object_size, fogbank_direction, perc_binning)


% get the image to be segmented
grayscale_image = single(grayscale_image);
% img = grayscale_image;
[nb_rows, nb_cols] = size(grayscale_image);

if nargin == 5
    fogbank_direction = 1; % 1 is min to max 
    perc_binning = 5;
end
if nargin == 6
    perc_binning = 5;
end

assert(islogical(foreground_mask), 'fog_bank_perctile_geodist_seed:argCk','Invalid <foreground_mask>, Type');
assert(size(foreground_mask,1) == nb_rows && size(foreground_mask,2) == nb_cols, 'fog_bank_perctile_geodist_seed:argCk','Invalid <foreground_mask>, wrong size');
assert(islogical(mask_matrix),'fog_bank_perctile_geodist_seed:argCk','Invalid <mask_matrix> Type');
assert(size(mask_matrix,1) == nb_rows && size(mask_matrix,2) == nb_cols, 'fog_bank_perctile_geodist_seed:argCk','Invalid <mask_matrix> wrong size');
assert(size(seed_image,1) == nb_rows && size(seed_image,2) == nb_cols, 'fog_bank_perctile_geodist_seed:argCk','Invalid <seed_image> wrong size');
assert(min_object_size > 0, 'fog_bank_perctile_geodist:argCk','Invalid <min_object_size>');
assert(perc_binning > 0 && perc_binning < 100, 'fog_bank_perctile_geodist:argCk','Invalid <percentile_binning>');

if exist('colors_vector', 'var')
    img = grayscale_image;
end

% transform background to nan and get minimum value on cell area
grayscale_image(~foreground_mask) = NaN;


% P_vec = 0:perc_binning:100;
P_vec = (0:perc_binning:100)/100;
Y = percentile_computation(grayscale_image(:),P_vec);
% Y = prctile(grayscale_image(:),P_vec);
if fogbank_direction
    Y = sort(Y,'ascend');
    % ensure that the first iteration comprises of only the seed pixels
    min_val = min(grayscale_image(seed_image>0)) - 1;
    grayscale_image(seed_image>0) = min_val;
    Y = [min_val, Y];
else
    Y = sort(Y,'descend');
    % ensure that the first iteration comprises of only the seed pixels
    max_val = max(grayscale_image(seed_image>0)) + 1;
    grayscale_image(seed_image>0) = max_val;
    Y = [max_val, Y];
end


% based on the distance transform matrix, gradually drop a fog from the sky down to the ground passing through
% all the mountains in between and keeping them separated from one another.

% Start dropping the fog 
for n = 1:numel(Y)
    % get the binary image containing the pixels that are to be assigned a label at this fog level
    if fogbank_direction
        image_b = grayscale_image <= Y(n) & mask_matrix;
    else
        image_b = grayscale_image >= Y(n) & mask_matrix;
    end
    
    % assign non zero pixels in image_b the label of the closest connected peak in seed_image
    seed_image = assign_nearest_connected_label(seed_image, image_b);
    if exist('colors_vector', 'var')
        disp_I = superimpose_big_image(img, seed_image,0,colors_vector, 0);
        imwrite(disp_I, ['FB_prog_' sprintf('%03d',n) '.tif']);
    end
    
end

% assign any un assigned pixels to the nearest body
seed_image = assign_nearest_connected_label(seed_image, foreground_mask);
nb_peaks = max(seed_image(:));

% Scout all the pixels in the image looking for the non background pixels
indx = find(seed_image);
objects_size = zeros(nb_peaks, 1);
for i = 1:numel(indx)
    % Increment the size of the pixel
    objects_size(seed_image(indx(i))) = objects_size(seed_image(indx(i))) + 1;
end

% Delete cells with size less than threshold
[seed_image, nb_peaks] = check_cell_size_renumber(seed_image, nb_peaks, objects_size, min_object_size);
% [seed_image, nb_peaks] = check_cell_size(seed_image, nb_peaks, objects_size, min_object_size);

% Check for objects connectivity in the image after performing the separation by the sticking fog.
seed_image = check_body_connectivity(seed_image, nb_peaks);

if exist('colors_vector', 'var')
    disp_I = superimpose_big_image(img, seed_image,0,colors_vector, 0);
    imwrite(disp_I, ['FB_prog_' sprintf('%03d',numel(Y)+1) '.tif']);
end


end



function dataType = get_min_required_datatype(maxVal)
if maxVal <= intmax('uint8')
        dataType = 'uint8';
    elseif maxVal <= intmax('uint16')
        dataType = 'uint16';
    elseif maxVal <= intmax('uint32')
        dataType = 'uint32';
    else
        dataType = 'double';
end
end


function [img, highest_cell_number] = check_cell_size(img, Highest_cell_number, cell_size, cell_size_threshold)

% Create a renumber_cells vector that contains the renumbering of the cells with size > min_size
renumber_cells = zeros(Highest_cell_number+1, 1);
highest_cell_number = 0;
for i = 1:Highest_cell_number
    renumber_cells(i+1) = 0;
    % if cell i is a cell with size > min_size, give it a new number
    if cell_size(i) > cell_size_threshold
        renumber_cells(i+1) = i;
    end
end

% Delete small cells 
BW = img > 0;
img = renumber_cells(img+1);
% assign deleted pixels the label of the nearest connected body
img = assign_nearest_connected_label(img, BW);

end



function [img, highest_cell_number] = check_cell_size_renumber(img, Highest_cell_number, cell_size, cell_size_threshold)

% Create a renumber_cells vector that contains the renumbering of the cells with size > min_size
renumber_cells = zeros(Highest_cell_number+1, 1);
highest_cell_number = 0;
for i = 1:Highest_cell_number
    % if cell i is a cell with size > min_size, give it a new number
    if cell_size(i) > cell_size_threshold
        highest_cell_number = highest_cell_number + 1;
        renumber_cells(i+1) = highest_cell_number;
    end
end

% Delete small cells 
BW = img > 0;
img = renumber_cells(img+1);
% assign deleted pixels the label of the nearest connected body
img = assign_nearest_connected_label(img, BW);

end





