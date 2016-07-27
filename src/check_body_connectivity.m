% Disclaimer:  IMPORTANT:  This software was developed at the National Institute of Standards and Technology by employees of the Federal Government in the course of their official duties. Pursuant to title 17 Section 105 of the United States Code this software is not subject to copyright protection and is in the public domain. This is an experimental system. NIST assumes no responsibility whatsoever for its use by other parties, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. We would appreciate acknowledgment if the software is used. This software can be redistributed and/or modified freely provided that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.

% Check for objects connectivity in the image
function image_out = check_body_connectivity(image_in, Highest_cell_number)

% check for the mex file binaries that are required to run the mex file
mex_file_name = ['check_body_connectivity_mex.' mexext];
if ~exist(mex_file_name, 'file') && exist('check_body_connectivity_mex.c','file')
    % check the compiler configuration, do not compile unless there is a selected C Compiler
    if ~isempty(mex.getCompilerConfigurations('C','Selected'))
        eval(['mex ' which('check_body_connectivity_mex.c')]);
    end
end
if exist(mex_file_name, 'file')
    image_out = check_body_connectivity_mex(image_in);
    return;
end

% Initialize the output
image_out = image_in;

% Get the image dimensions
[nb_rows, nb_cols] = size(image_out);

% print_update(1,1,Highest_cell_number);
% Scout all the objects in image_out
for k = 1:Highest_cell_number
%     print_update(2,k,Highest_cell_number);
    % Create the binary image where all objects are set to 0 except for object i
    image_b = image_out == k;

    % label the bodies in the binay image
    [labeled_image, nb_objects] = bwlabel(image_b);

    % if only one body found: continue
    if nb_objects == 1, continue, end

    % if more than 1 body is found initialize the size vector for all the bodies found
    objects_size = zeros(nb_objects, 1);

    % compute the size of each body
    indx = find(labeled_image);
    for i = 1:numel(indx), objects_size(labeled_image(indx(i))) = objects_size(labeled_image(indx(i))) + 1; end

    % Keep the body with the biggest size
    [max_size, winner_body] = max(objects_size);         %#ok<ASGLU>
    
    % Create the matrix body_neighbors that holds all the numbers of the neighbors of each body
    body_neighbors = zeros(Highest_cell_number, nb_objects);
    
    % Scout all labeled_image looking for the neighbors (in image_out) of each body except the winner body
    for j = 2:nb_cols-1
        for i = 2:nb_rows-1
            
            % if pixel(i,j) is a background pixel in labeled_image or belongs to the winner body: continue
            if labeled_image(i,j) == 0 || labeled_image(i,j) == winner_body, continue, end
            
            % Check if the left neighbor pixel is not the background and is not object k in image_out
            pixel = image_out(i,j-1);
            if pixel > 0 && pixel ~= k
                body_neighbors(pixel, labeled_image(i,j)) = body_neighbors(pixel, labeled_image(i,j)) + 1;
            end
            
            % Check if the top neighbor pixel is not the background and is not object k in image_out
            pixel = image_out(i-1,j);
            if pixel > 0 && pixel ~= k
                body_neighbors(pixel, labeled_image(i,j)) = body_neighbors(pixel, labeled_image(i,j)) + 1;
            end

            % Check if the right neighbor pixel is not the background and is not object k in image_out
            pixel = image_out(i,j+1);
            if pixel > 0 && pixel ~= k
                body_neighbors(pixel, labeled_image(i,j)) = body_neighbors(pixel, labeled_image(i,j)) + 1;
            end

            % Check if the top neighbor pixel is not the background and is not object k in image_out
            pixel = image_out(i+1,j);
            if pixel > 0 && pixel ~= k
                body_neighbors(pixel, labeled_image(i,j)) = body_neighbors(pixel, labeled_image(i,j)) + 1;
            end
        end
    end
    
    % Find the dominant neighbor of each body
    [max_neighbor, neighbor_winner] = max(body_neighbors);
    
    % Re-Scout all image_out and change the number of the bodies of object k to the neighbor_winner except for the
    % winner body which must remain number k
    for i = 1:numel(labeled_image)
        % if pixel(i,j) is a background pixel in labeled_image or belongs to the winner body: continue
        if labeled_image(i) == 0 || labeled_image(i) == winner_body, continue, end

        % if no neighbor is found for body with number labeled_image(i,j), delete the pixel and continue
        if max_neighbor(labeled_image(i)) == 0, image_out(i) = 0; continue, end

        % Otherwise renumber the body to the winner_neighbor
        image_out(i) = neighbor_winner(labeled_image(i));
    end
end
% print_update(3,Highest_cell_number,Highest_cell_number);


