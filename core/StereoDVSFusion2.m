function [RV_fusion, pc_view] = StereoDVSFusion2(Estimations, ViewID, VIP_RVs, N, cam_left, useMedianFilter)
    % number of estimations
    numEstimations = length( Estimations );
    numVIP_RV = length(VIP_RVs);
    RV_fusion = cell(numVIP_RV,1);
    pc_view = [];
    
    % pose of the ViewID_th view.
    T_w_v = Estimations{ViewID}.pose;
    if size(T_w_v, 1) == 3
        T_w_v = [T_w_v;0 0 0 1];
    end
    invT_w_v = inv(T_w_v);
    
    % Perform fusion on each vip rv with its 2 * N neighbouring estimations.
    for it = 1:numVIP_RV
        i = VIP_RVs(it);
        nbRVs = max(1, i - N):min(numEstimations,i + N);
        nbRVs(nbRVs == i) = [];
        RV_fusion{it} = Estimations{i};
        for j = nbRVs
            RV_fusion{it} = InvDepthMapFusion(RV_fusion{it}, Estimations{j}, cam_left);
        end
        % median filter?
        if useMedianFilter
            A = RV_fusion{it}.invDepthMap(:,:,1);
            fun = @(x) nanmedian(x(:));
            B = nlfilter(A,[1 1]*15,fun).*~isnan(A);
            B(B==0) = NaN;
            RV_fusion{it}.invDepthMap(:,:,1) = B;
        end
        % update pc
        RV_fusion{it}.p3D = InvDepthMap2Pointcloud(RV_fusion{it}.invDepthMap, cam_left);
        
        % accumulate all pointcloud to ViewID_th view.
        T_w_i = Estimations{i}.pose;
        pc_k = Estimations{i}.p3D(1:3,:);
        if size(T_w_i, 1) == 3
            T_w_i = [T_w_i;0 0 0 1];
        end
        T_v_k = invT_w_v * T_w_i;
        pc = T_v_k(1:3,1:3) * pc_k + T_v_k(1:3,4);
        pc = [pc; Estimations{i}.p3D(4:5,:)];
%         StereoDVS_VisualizePointcloud(pc, [], []);
        pc_view = [pc_view, pc];
%         StereoDVS_VisualizePointcloud(pc_view, [], []);
        disp(['Fusion: the ',num2str(it),' RV']);
    end
end

function RV_fusion = InvDepthMapFusion(estimation1, estimation2, cam_left)
    RV_fusion = estimation1;
    %
    T_w_1 = estimation1.pose;
    if size(T_w_1, 1) == 3
        T_w_1 = [T_w_1;0 0 0 1];
    end
    invT_w_1 = inv(T_w_1);
    %
    T_w_2 = estimation2.pose;
    if size(T_w_2, 1) == 3
        T_w_2 = [T_w_2;0 0 0 1];
    end
    T_1_2 = invT_w_1 * T_w_2;
    %
    p3D_1 = T_1_2(1:3,1:3) * estimation2.p3D(1:3,:) + T_1_2(1:3,4);
    p3D_1 = [p3D_1;estimation2.p3D(4:5,:)];
    numPoint = size(p3D_1,2);
    for j = 1:numPoint
        x_1 = cam_left.world2Cam_rect(p3D_1(1:3,j));
        rho  = 1.0 / p3D_1(3,j);
        var  = p3D_1(4,j);
        cost = p3D_1(5,j);
        % fusion
        RV_fusion.invDepthMap = Fusion(RV_fusion.invDepthMap, x_1, rho, var, cost);
    end
end

% % Kalman fusion
function Map = Fusion(Map, x, rho, var, cost)
u = floor(x(1));
v = floor(x(2));
[height, width, channel] = size(Map);
if u < 1 || u > width - 1 || v < 1 || v > height - 1
    return;
end

candidate = [[u;v;Map(v,u,1);Map(v,u,2);Map(v,u,3)],         [u+1;v;Map(v,u+1,1);Map(v,u+1,2);Map(v,u+1,3)], ...
             [u;v+1;Map(v+1,u,1);Map(v+1,u,2);Map(v+1,u,3)], [u+1;v+1;Map(v+1,u+1,1);Map(v+1,u+1,2);Map(v+1,u+1,3)]];

group1 = find( isnan( candidate(3,:) ) ); %non-occupied pixels.
group2 = find( ~isnan( candidate(3,:) ) );%occupied pixels.
group3 = [];%compatible
group4 = [];%non-compatible
for i = group1
    Map(candidate(2, i), candidate(1, i), 1) = rho;
    Map(candidate(2, i), candidate(1, i), 2) = var;
    Map(candidate(2, i), candidate(1, i), 3) = cost;
end

% compatibles and non-compatibles
for i = group2
    if ChiSquareTest(rho, candidate(3,i), var, candidate(4,i))
        group3 = [group3, i];
    else
        group4 = [group4, i];
    end
end

% fusion
if numel(group3) > 0
    Mus  = [rho, candidate(3, group3)];
    Vars = [var, candidate(4, group3)];
    [rho_new, var_new] = FuseGaussian(Mus, Vars);
    Map(candidate(2, group3), candidate(1, group3), 1) = rho_new;
    Map(candidate(2, group3), candidate(1, group3), 2) = var_new;
    Map(candidate(2, group3), candidate(1, group3), 3) = cost;
end

% replace the non-compatible ones whose cost and var are higher than the
% new observation.
for i = group4
    if candidate(4, i) > var %&& candidate(5, i) > cost
        Map(candidate(2, i), candidate(1, i), 1) = rho;
        Map(candidate(2, i), candidate(1, i), 2) = var;
        Map(candidate(2, i), candidate(1, i), 3) = cost;
    end
end
end

% % Fuse Gaussion Distribution
function [mu, var] = FuseGaussian(Mus, Vars)
    numDistribution = numel(Mus);
    nominator = 0;
    denominator = 0;
    for i = 1:numDistribution
        nominator = nominator + Mus(i) / Vars(i);
        denominator = denominator + 1 / Vars(i);
    end
    mu = nominator / denominator;
    var = 1 / denominator;
end

% % Chi square test
function bCompatible = ChiSquareTest(rho1, rho2, var1, var2)
    delta_rho = rho1 - rho2;
    compatibility = delta_rho^2/var1 + delta_rho^2/var2;
    if compatibility < 5.99
        bCompatible = true;
    else
        bCompatible = false;
    end
end