% Disclaimer:  IMPORTANT:  This software was developed at the National Institute of Standards and Technology by employees of the Federal Government in the course of their official duties. Pursuant to title 17 Section 105 of the United States Code this software is not subject to copyright protection and is in the public domain. This is an experimental system. NIST assumes no responsibility whatsoever for its use by other parties, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. We would appreciate acknowledgment if the software is used. This software can be redistributed and/or modified freely provided that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.



function BW = cluster_objects(BW, cluster_distance, valid_traversal_mask)

if nargin == 2
    valid_traversal_mask = true(size(BW));
end

[BW, nb_obj] = bwlabel(BW);
cluster_labels = zeros(nb_obj,1);
% generate the centroids image to use as the seeds for geodesic distance
stats = regionprops(BW, 'Centroid', 'BoundingBox');

centroids = zeros(size(stats,1),2);
for i = 1:numel(stats)
    centroids(i,:) = stats(i).Centroid;
end
centroids = round(centroids);


BWd = bwdistgeodesic(valid_traversal_mask, centroids(:,1), centroids(:,2));
clusters = bwlabel(BWd <= cluster_distance/2);
% assign clusters labels
for i = 1:size(centroids,1)
    cent = [centroids(i,2), centroids(i,1)];
    cluster_label = clusters(cent(1), cent(2));
    
    cluster_labels(i) = cluster_label;
end

cluster_labels = [0;cluster_labels];
% relabel the image
BW = cluster_labels(BW+1);


end