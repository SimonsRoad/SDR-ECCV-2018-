% Camera model used in the project. This is specific to the equidistant distortion model.
classdef Camera
    properties
        width
        height
        f
        fx
        fy
        cx
        cy
        distortion_model
        D
        K
        invK
        R
        P
        undistortionMap_X
        undistortionMap_Y
    end
    
    methods
        function cam = Camera( calibFileName, UndistortionMapX_FileName, UndistortionMapY_FileName )
            calibInfo = YAML.read( calibFileName );
            cam.width = calibInfo.width;
            cam.height = calibInfo.height;
            cam.fx = calibInfo.K(1,1);
            cam.fy = calibInfo.K(2,2);
            cam.cx = calibInfo.K(1,3);
            cam.cy = calibInfo.K(2,3);
            cam.f = ( calibInfo.K(1,1) + calibInfo.K(2,2) ) / 2;
            cam.distortion_model = calibInfo.distortion_model;
            cam.D = calibInfo.D;
            cam.K = calibInfo.K;
            cam.invK = inv(cam.K);
            cam.R = calibInfo.R;
            cam.P = calibInfo.P;
            if ~isempty(UndistortionMapX_FileName) || ~isempty(UndistortionMapY_FileName)
                cam.undistortionMap_X = importdata(UndistortionMapX_FileName);
                cam.undistortionMap_Y = importdata(UndistortionMapY_FileName);
            else
                cam.undistortionMap_X = [];
                cam.undistortionMap_Y = [];
            end
        end
        
        % back-project a rectified 2D pixel to 3D, given the inverse depth
        % inv_d = 1/z;
        function p = cam2World_rect( obj, pixel, inv_d )
            z = 1 / inv_d;
            x_ss = [pixel; 1; 1];
            P_tilde = [obj.P;0,0,0,z];
            p_s = z * P_tilde \ x_ss;
            p = p_s(1:3) / p_s(4);
        end
        
        function x = world2Cam_rect( obj, p )
            x_hom = obj.P * [p;1];
            x = x_hom / x_hom(3);
        end        
    end
end