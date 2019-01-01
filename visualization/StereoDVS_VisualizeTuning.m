function InvDepthMap_RVs = StereoDVS_VisualizeTuning(InvDepthMap_RVs, cam_left, varThreshold, costThreshold)
    %
    numRV = length(InvDepthMap_RVs);
    for i = 1:numRV
        var_trunc = prctile( InvDepthMap_RVs{i}.p3D(4,:), varThreshold );
        InvDepthMap_RVs{i}.p3D(:, InvDepthMap_RVs{i}.p3D(4,:) > var_trunc ) = [];
        cost_trunc = prctile( InvDepthMap_RVs{i}.p3D(5,:), costThreshold );
        InvDepthMap_RVs{i}.p3D(:, InvDepthMap_RVs{i}.p3D(5,:) > cost_trunc ) = [];
        
        % update InvDepthMap
        InvDepthMap_RVs{i}.invDepthMap = PointCloud2InvDepthMap(InvDepthMap_RVs{i}.p3D, cam_left);
    end
end