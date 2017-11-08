close all; dbstop error; clc; clearvars

%% Load OXTS data for each frame
base_dir = 'PathToData\2011_10_03\2011_10_03_drive_0027_sync';
% frames
start_frame = 1; 
end_frame = 100;
% Points too far away can be ignored, the borders are set here
Left = -120; Right = 120; Front = 120; Back = -120; Low_Height = -5;
% 10Hz for the KITTI dataset
framerate = 10;

frames = start_frame:end_frame;
nFrames = length(frames);
% Using the devkit toolbox
oxts = loadOxtsliteData(base_dir,frames);

%% The oxts include in order:
% 1.  lat:   latitude of the oxts-unit (deg)
% 2.  lon:   longitude of the oxts-unit (deg)
% 3.  alt:   altitude of the oxts-unit (m)
% 4.  roll:  roll angle (rad),    0 = level, positive = left side up,      range: -pi   .. +pi
% 5.  pitch: pitch angle (rad),   0 = level, positive = front down,        range: -pi/2 .. +pi/2
% 6.  yaw:   heading (rad),       0 = east,  positive = counter clockwise, range: -pi   .. +pi
% 7.  vn:    velocity towards north (m/s)
% 8.  ve:    velocity towards east (m/s)
% 9.  vf:    forward velocity, i.e. parallel to earth-surface (m/s)
% 10. vl:    leftward velocity, i.e. parallel to earth-surface (m/s)
% 11. vu:    upward velocity, i.e. perpendicular to earth-surface (m/s)
% 12. ax:    acceleration in x, i.e. in direction of vehicle front (m/s^2)
% 13. ay:    acceleration in y, i.e. in direction of vehicle left (m/s^2)
% 14. ay:    acceleration in z, i.e. in direction of vehicle top (m/s^2)
% 15. af:    forward acceleration (m/s^2)
% 16. al:    leftward acceleration (m/s^2)
% 17. au:    upward acceleration (m/s^2)
% 18. wx:    angular rate around x (rad/s)
% 19. wy:    angular rate around y (rad/s)
% 20. wz:    angular rate around z (rad/s)
% 21. wf:    angular rate around forward axis (rad/s)
% 22. wl:    angular rate around leftward axis (rad/s)
% 23. wu:    angular rate around upward axis (rad/s)
% 24. pos_accuracy:  position accuracy (north/east in m)
% 25. vel_accuracy:  velocity accuracy (north/east in m/s)
% 26. navstat:       navigation status (see navstat_to_string)
% 27. numsats:       number of satellites tracked by primary GPS receiver
% 28. posmode:       position mode of primary GPS receiver (see gps_mode_to_string)
% 29. velmode:       velocity mode of primary GPS receiver (see gps_mode_to_string)
% 30. orimode:       orientation mode of primary GPS receiver (see gps_mode_to_string)

% Translations of variables to the "LiDAR point clouds correction acquired 
% from a moving car based on CAN-bus data" paper
% Dtheta = wu
% Dx = vf

%% For each frame, correct the points here...
for frame=start_frame:end_frame
    %% Read raw points -----------
    fid = fopen(sprintf('%s/velodyne_points/data/%010d.bin',base_dir,frame),'rb');
    velo = fread(fid,[4 inf],'single')';
    fclose(fid);
    velo = velo(1:5:end,1:3); % remove every 5th point for display speed
    idx = ((sqrt((velo(:,1)).^2 + (velo(:,2)).^2)) < 3); % filter false points in the blind spot
    velo(idx,:) = [];
    % remove all points too far left, right and behind the camera
    idx = (velo(:,2) < Left) | (velo(:,2) > Right) | (velo(:,1) < Back) | (velo(:,1) > Front) | (velo(:,3) < Low_Height);
    velo(idx,:) = [];
    nPoints = size(velo,1);
    
    % Assign scan line for the inclination angle omega
    [ slices_ind ] = organize_velo_cloud( velo );
        
    %% Car state at each frame
    Dtheta = oxts{frame+nFrames-end_frame}(23) / framerate;
    Dx = oxts{frame+nFrames-end_frame}(9) / framerate;
    
    % Initialize Visualization
    figure, hold on, axis equal
    dot_size = 2;
    car_radius = 5;
    % plot axis (x-red, y-green, z-blue)
    plot3([0 car_radius], [0 0], [0 0], 'LineWidth',2, 'Color', [1 0 0])
    plot3([0 0], [0 car_radius], [0 0], 'LineWidth',2, 'Color', [0 1 0])
    plot3([0 0], [0 0], [0 car_radius], 'LineWidth',2, 'Color', [0 0 1])
    xlabel('X'), ylabel('Y'), zlabel('Z')
    
    % Initialize omega angle based on the first scan line
    flag_scan = 1;
    % 26.8 / 64 = 0.4 degrees for omega
    % We compute one omega for each scan line using atan2(z,x)
    [omega,~] = cart2pol(velo(1,1),velo(1,3));
    omega = wrapTo2Pi(omega);
    % Each point has a different orientation, so we need to compute the
    % proper LiDAR position for each point and then correct the point with
    % respect to the correct LiDAR state.
    velo_corr = zeros(size(velo));
    %% For each point in the frame...
    for pt=1:nPoints
        % alpha is the azimuth angle, we compute it based on the atan2(y,x)
        [alpha,~] = cart2pol(velo(pt,1),velo(pt,2));
        alpha = wrapTo2Pi(alpha);
        % Distance from point to LiDAR
        d = norm(velo(pt,:));
        
        % See if we changed scan line. If yes, then the omega needs to
        % change as well
        if flag_scan ~= slices_ind(pt)
            [omega,~] = cart2pol(velo(pt,1),velo(pt,3));
            omega = wrapTo2Pi(omega);
            flag_scan = slices_ind(pt);
            disp('Change scan line')
        end
        % LiDAR position correction when the lidar finishes a full rotation
        % using the Runge-Kutta equations
        DxAlpha = -Dx * ( alpha  / (2*pi) );
        DthetaAlpha = -Dtheta * ( alpha  / (2*pi) );
              
        thetaAlpha = -DthetaAlpha;
        Xa = -DxAlpha * cos( thetaAlpha - DthetaAlpha / 2 ); 
        Ya = -DxAlpha * sin( thetaAlpha - DthetaAlpha / 2 );
        % Plot the motion of the LiDAR within a frame
        % plot(Xa, Ya, 'm.')
        
        % Point correction using the new position of the LiDAR
        Ralpha = [cos(thetaAlpha) -sin(thetaAlpha) 0;
                  sin(thetaAlpha)  cos(thetaAlpha) 0;
                  0                0               1]; %yaw rotation
        velo_corr(pt,:) = ([Xa; Ya; 0] + ( Ralpha * ...
            [cos(omega)*cos(alpha); cos(omega)*sin(alpha); sin(omega)] * d) )';
    end
    
    scatter3(velo_corr(:,1),velo_corr(:,2),velo_corr(:,3), 5, 'b')
    scatter3(velo(:,1),velo(:,2),velo(:,3), 5, 'r')
    legend('Xaxis','Yaxis','Zaxis','corrected','original')
    pause
    close all
end


