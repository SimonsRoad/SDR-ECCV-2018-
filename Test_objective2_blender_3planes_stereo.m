% This script is to evaluate the effectiveness of the proposed objective
% function. The synthetic dataset "blender_3planes_stereo" is used.
%%
clc;clear;close all;
% add /wrapper
addpath([pwd,'/wrapper']);

%% choose a sequence and TUNE paremeters.
% blender
load(fullfile('data', 'observation_3planes_stereo_small_baseline.mat'));
disp('Data loaded ...');

% arbitrarily set a reference view ( the numbers are observations' id )
start = 22;
ending = 40;
RV_id = start;

% parameters (Feel free to tune.)
w = 35;
Lnorm = 'l2';
DisplayOptimizationProcess = 0 == 0;

%% info about the observations
width = cam_left.width;
height = cam_left.height;
numObservation = ending - start + 1;
% Time
T_start  = TimeStamps_left{start};
T_end    = TimeStamps_left{ending};
T_RV     = TimeStamps_left{RV_id};
% Absolute pose of the RV
T_w_rv   = Poses_left{RV_id};

disp(['The reference view (RV) is captured at ', num2str(T_RV), ' s']);
disp(['Observations start at ', num2str(T_start), ' s, ends at ', num2str(T_end), ' s.']);
disp(['All observations last for ', num2str(T_end - T_start), 's']);

% load groundtruth info at the reference view. Note that the image is not
% required by the proposed algorithm. It is just for visualization.
% THE FOLLOWING TWO IMAGES MUST BE CONSISTENT WITH THE RV_id, OTHERWISE THE
% VISUALIZATION DOES NOT MAKE SENSE!!!
RV_Image    = imread('data/3planes_stereo_small_baseline/left/image_raw/0.209000.png');
RV_DepthMap = imread('data/3planes_stereo_small_baseline/left/depth_raw/0.209000.png');

% % get events in the RV
RV_EMaps = EMaps_left{RV_id};
x = EventLists_left{RV_id}';
numEvents = size(x,2);

% compute transformations Ts_left_rv (i.e the relative pose of each 
% observation w.r.t the reference view).
Ts_left_rv = cell(numObservation,1);
for i = 1:numObservation
    R_left_w  = Poses_left{start+i-1}(:,1:3)';
    t_left_w  = -R_left_w * Poses_left{start+i-1}(:,4);
    R_left_rv = R_left_w * T_w_rv(:,1:3);
    t_left_rv = t_left_w + R_left_w * T_w_rv(:,4);
    Ts_left_rv{i} = [R_left_rv, t_left_rv];
end

%% Illustrate the objective function of each event in the reference frame.
for i = 1:100:numEvents % I set the step as 100 here to skip neighbour events. Feel free to tune.
    % boundary check
    if( x(1,i) < 40 || x(1,i) > width - 40 || x(2,i) < 40 || x(2,i) > height - 40 ) 
        continue;
    end
    
    % groundtruth depth at each given event location.
    depth = double(RV_DepthMap(floor(x(2,i)),floor(x(1,i)))) / 5000.0;
    d_gt = 1 / depth;

    % compute cost at each hypothesis inverse depth
    d_iter = 0.2 : 0.01: 1;% mapping scope: 1 m ~ 5 m. Feel free to tune.
    
    % The matlat mex interface. Please check out the wrapper/*.cpp files for details.
    % By the way, the SAE here refers to the Time-Surface Map in the paper.
    Cr = EPTAM_mapping_mex('ComputeObjective', x(:,i), d_iter, numObservation, width, height, ...
                            SAEs_smooth_left(start:ending), SAEs_smooth_right(start:ending), ...
                            Ts_left_rv,...
                            cam_left.P, cam_right.P,...
                            w, Lnorm);

    % display
    figure;
    subplot(2,3,1);
    plot(d_iter, Cr, 'r-','LineWidth', 5);
    y_lim = 1e+7;%max(Cr);
    ylim([0, y_lim]);
    xlim([0.2, 1]);
    title('Objective function');
    hold on;
    m = y_lim;
    yy = 0:m/20:m;
    xx = d_gt * ones(size(yy,2),1);
    plot(xx,yy,'k--','Linewidth',3);%draw groundtruth inverse depth.
    xlabel('Inverse depth ({m^{-1}})');
    ylabel('Energy');

    subplot(2,3,2);
    imshow(RV_Image);
    hold on;
    plot(x(1,i), x(2,i), 'ro', 'Linewidth', 2)
    title('The selected event in RV (Image)');

    subplot(2,3,3);
    imshow(RV_EMaps);
    hold on;
    plot(x(1,i), x(2,i), 'ro', 'Linewidth', 2)
    title('The selected event in RV (Event Map)');

    % compute the reprojections of the groundtruth 3D point
    p_left = Ts_left_rv{end}(:,1:3) * cam_left.cam2World_rect(x(:,i), d_gt) + Ts_left_rv{end}(:,4);
    x1_s = cam_left.P * [p_left;1];
    x1_s = x1_s / x1_s(3);
    x2_s = cam_right.P * [p_left;1];
    x2_s = x2_s / x2_s(3);

    subplot(2,3,4);
    image(SAEs_left{ending});
    hold on;
    plot(x1_s(1), x1_s(2),'sr','MarkerSize',w/2);
    title('Time-Surface Map (left)');

    subplot(2,3,5);
    image(SAEs_right{ending});
    hold on;
    plot(x2_s(1), x2_s(2),'sr', 'MarkerSize', w/2);
    title('Time-Surface Map (right)');
    
    %
    disp(['This is the ', num2str(i), ' th event.']);
end