function [EMap_RVs, InvDepthMap_RVs, cam] = StereoDVSMapping(Observations, numRV, RV_ids, obScope)
% Introduction:
% numRV: when set to zero, the RV_indices will be directly obtained from
% RV_ids.

    addpath('/home/zhouyi-anu/workspace/project/EPTAM_matlab/wrapper');
    % load observations
    load(Observations);
    cam = cam_left;% we use this camera outsi
    
    % parameters (fixed)
    w = 25;
    Lnorm = 'l2';
    outlier_percentile = 100;
    if isempty(obScope)
        observationScope = 15;
    else
        observationScope = obScope;
    end

    % choose RVs
    if ~isempty(numRV)
        RV_step = floor(numObservations / (numRV + 1));
        RV_indices = RV_step:RV_step:numObservations-1;
    else
        RV_indices = RV_ids;
        numRV = numel(RV_ids);
    end
    
    disp('****************************');
    disp('loaded data');
    disp('data information:');
    disp(['- data name: ', Observations]);
    disp(['- number of observations:', num2str(numObservations)]);
    disp(['- width: ', num2str(width)]);
    disp(['- height: ', num2str(height)]);
    disp(['- w: ', num2str(w)]);
    disp(['- Lnorm: ', Lnorm]);
    disp(['- RV_indices: ', num2str(RV_indices)]);
    disp('****************************');
    
    % allocate results
    EMap_RVs = cell(numRV, 1);
    InvDepthMap_RVs = cell(numRV, 1);
    
    % estimate inverse depth map for each RV
    for i = 1:numRV
        start = max(RV_indices(i) - observationScope + 1, 1);
        ending = min(RV_indices(i) + observationScope - 1, numObservations);
        numObservation = ending - start + 1;
        T_start  = TimeStamps_left{start};
        T_end    = TimeStamps_left{ending};
        T_RV     = TimeStamps_left{RV_indices(i)};
        T_w_rv   = Poses_left{RV_indices(i)};
        disp(['The RV is captured at ', num2str(T_RV)]);
        disp(['Observations start at ', num2str(T_start), ' s, ends at ', num2str(T_end), ' s.']);
        disp(['All involved observations take ',num2str(T_end - T_start), ' s.']);
        
        %
        EMap_RVs{i} = EMaps_left{i};
        
        % get events in the RV
        x = EventLists_left{RV_indices(i)}';
        numEvents = size(x,2);
        
        % compute poses Poses_left_rv in advance
        Ts_left_rv = cell(numObservation,1);
        for j = 1:numObservation
            R_left_w  = Poses_left{start+j-1}(:,1:3)';
            t_left_w  = -R_left_w * Poses_left{start+j-1}(:,4);
            R_left_rv = R_left_w * T_w_rv(:,1:3);
            t_left_rv = t_left_w + R_left_w * T_w_rv(:,4);
            Ts_left_rv{j} = [R_left_rv, t_left_rv];
        end
    
        % estimation
        tic;
        InvDMapList = EPTAM_mapping_mex('EstimateInvDepthMap_parallel', numEvents, x, numObservation, width, height, ...
                                         SAEs_smooth_left(start : ending), SAEs_smooth_right(start : ending), ...
                                         dSAEs_du_left(start : ending), dSAEs_dv_left(start : ending), ...
                                         dSAEs_du_right(start : ending), dSAEs_dv_right(start : ending), ...
                                         Ts_left_rv, ...
                                         cam_left.P, cam_right.P, ...
                                         w, Lnorm );

        disp(['Finish ', num2str(i), 'th RV''s estimation']);
        toc
        
        % invDepthMap is a height * width * 3 matrix, the first level
        % (:,:,1) is invDepth; the second level (:,:,2) is variance; 
        %the third level (:,:,3) is remaining cost.
        InvDMapList = vec2mat(InvDMapList, 5)';
        % clean
        InvDMapList(:, InvDMapList(3,:) < 0.10) = [];% < 10 m
        InvDMapList(:, InvDMapList(3,:) > 3) = [];% > 0.3 m
        % occlusion rejection (by checking the remaining cost)
        occlusion_threshold = prctile( InvDMapList(5,:), outlier_percentile );
        InvDMapList(:, InvDMapList(5,:) > occlusion_threshold) = [];
        % variance rejection
        var_threshold = prctile( InvDMapList(4,:), outlier_percentile );
        InvDMapList(:, InvDMapList(4,:) > var_threshold) = [];
        
        % recover 3D
        numPoint = size(InvDMapList,2);
        X_3D = zeros(5,numPoint);
        for k = 1:numPoint
            pixel = InvDMapList(1:2,k);
            rho = InvDMapList(3,k);
            X_3D(1:3,k) = cam_left.cam2World_rect(pixel, rho);
            X_3D(4,k) = InvDMapList(4,k);%var
            X_3D(5,k) = InvDMapList(5,k);%remaining cost
        end
        
        % result
        InvDepthMap_RVs{i}.pose = T_w_rv;
        InvDepthMap_RVs{i}.p3D = X_3D;
        InvDepthMap_RVs{i}.invDepthMap = PointCloud2InvDepthMap( X_3D, cam_left );% width * height * 3 map
        InvDepthMap_RVs{i}.EMap = EMaps_left{RV_indices(i)};
        InvDepthMap_RVs{i}.timestamp = TimeStamps_left{RV_indices(i)};
        InvDepthMap_RVs{i}.obsID = RV_indices(i); 
    end
end