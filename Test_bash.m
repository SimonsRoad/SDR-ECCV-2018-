% % This script shows the whole story of our ECCV 2018 paper.
% It includes two parts:

%% PART I: Iverse Depth Estimation.
clc;clear;close all;
addpath([pwd,'/core']);
numRV = [];

% parameters
Observation_dir = [pwd,'/data/'];
sequenceName = 'observation_3planes_stereo_small_baseline';
RV_ids = 2:10:80;

obScope = 10;% regarding hwo many neighbour observations are used for the reconstruction.

% reconstruct events at the selected referen view (RV_ids).
[EMap_RVs, InvDepthMap_RVs, cam_left] = StereoDVSMapping([Observation_dir, sequenceName, '.mat'], [], RV_ids, obScope);

result_name = ['result_', sequenceName,'.mat'];
save(result_name);

%% PART II: Fusion.
% clc;clear;close all;
% sequenceName = 'observation_3planes_stereo_small_baseline';
% addpath([pwd,'/core']);
addpath([pwd,'/visualization']);
useMedianFilter = false;
N = 4;
load(['result_', sequenceName, '.mat']);

% visualization parameter
varTheshold   = 80;
costThreshold = 90;
% The following paramters are manualy chosen to find a good viewpoint for visualization.
numEstimations = length(InvDepthMap_RVs);
viewID = floor(numEstimations/2);
vip_rv = 2:2:8;
VisRange = [2, 3];
ViewAngle = [0 0 -1];

% 3D pointcloud
Estimations = StereoDVS_VisualizeTuning(InvDepthMap_RVs, cam_left, varTheshold, costThreshold);
[RV_fusion, pc_view] = StereoDVSFusion2(Estimations, viewID, vip_rv, N, cam_left, useMedianFilter);
StereoDVS_VisualizePointcloud(pc_view(1:3,:), VisRange, ViewAngle, [], []);

% Depth map
StereoDVS_VisualizeDepthMap(Estimations{viewID}.invDepthMap, VisRange, [], []);