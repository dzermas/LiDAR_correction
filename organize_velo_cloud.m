function [ slices_ind ] = organize_velo_cloud( velo )
%ORGANIZE_VELO_CLOUD Takes the unorganized velodyne measurements as input and 
%   outputs an index array that corresponds each point to a slice of the
%   velodyne sensor.

% The method is based on the fact that the velo(i,2)(y coordinate) goes from positive
% values to negative, and when it reaches positive again then a full cycle
% has been reached. Same solution can be reached by computing angles (ex.
% atan2(y,x)), but binary computations are cheaper.

flag = 0;
slice = 1;
slices_ind = zeros(length(velo),1);
for i=1:size(velo,1)-1
    slices_ind(i) = slice;
     
    % Zero values blow the logic, so I am changing them by a tiny bit
    if velo(i+1,2) == 0
        velo(i+1,2) = -0.0001;
    end
    
    if velo(i,2)*velo(i+1,2) < 0
        flag = flag + 1;
        if mod(flag,2) == 0
            slice = slice + 1;
        end
    end
end
slices_ind(end) = slices_ind(end-1);

end

