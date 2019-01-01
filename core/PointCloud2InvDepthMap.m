function InvDepthMap = PointCloud2InvDepthMap( pc, cam )
    width = cam.width;
    height = cam.height;
    InvDepthMap = NaN(height, width, 3);
    numPoint = size(pc, 2);
    for i = 1:numPoint
        point = pc(:,i);
        pixel = cam.world2Cam_rect(point(1:3));
        u = floor(pixel(1));
        v = floor(pixel(2));
        rho = 1.0 / point(3);
        var = point(4);
        cost = point(5);
        if rho > 3 || rho < 0.16
            continue;
        end
        
        % assign rho, var, cost to neighbouring four interger pixel
        % coordinate
        % u, v
        if InsideImage(u, v, width, height)
            InvDepthMap(v, u, 1) = rho;
            InvDepthMap(v, u, 2) = var;
            InvDepthMap(v, u, 3) = cost;
        end
        % u + 1, v
        if InsideImage(u + 1, v, width, height)
            InvDepthMap(v, u + 1, 1) = rho;
            InvDepthMap(v, u + 1, 2) = var;
            InvDepthMap(v, u + 1, 3) = cost;
        end
        % u, v + 1
        if InsideImage(u, v + 1, width, height)
            InvDepthMap(v + 1, u, 1) = rho;
            InvDepthMap(v + 1, u, 2) = var;
            InvDepthMap(v + 1, u, 3) = cost;
        end
        % u + 1, v + 1
        if InsideImage(u + 1, v + 1, width, height)
            InvDepthMap(v + 1, u + 1, 1) = rho;
            InvDepthMap(v + 1, u + 1, 2) = var;
            InvDepthMap(v + 1, u + 1, 3) = cost;
        end
    end
end

function bFlag = InsideImage(u, v, width, height)
    if u < 1 || u > width || v < 1 || v > height
        bFlag = false;
    else
        bFlag = true;
    end
end