function pc = InvDepthMap2Pointcloud(invDepthMap, cam)
    % invDepthMap h * w * 3
    [height, width, channel] = size(invDepthMap);
    validPoint = ~isnan(invDepthMap(:,:,1));
    pc = zeros(3, sum(validPoint(:)));
    p_index = 0;
    for y = 1:height
        for x = 1:width
            if ~isnan(invDepthMap(y,x,1))
                pixel = [x; y];
                rho = invDepthMap(y,x,1);
                p_index = p_index + 1;
                if rho == 0
                    continue;% this happen because of the median filter.
                end
                pc(:,p_index) = cam.cam2World_rect(pixel, rho);
            end
        end
    end
    pc(:, pc(3,:) < 0.5 ) = [];
    pc(:, pc(3,:) > 6 ) = [];
end