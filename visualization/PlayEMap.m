% % This script is used for visualizing the EMap or SAEs.
% 
% %
target_sequence = 8;

% parameters
Observation_dir = '/home/zhouyi-anu/workspace/project/EPTAM_matlab/observations/';

switch target_sequence
    case 1  
        sequenceName = 'observation_3planes_stereo_small_baseline';
    case 2
        sequenceName = 'observation_primitives_stereo';
    case 3
        sequenceName = 'observation_rpg_DAVIS_stereo_desk2';
    case 4
        sequenceName = 'observation_rpg_DAVIS_stereo_boxes2';
    case 5
        sequenceName = 'observation_rpg_DAVIS_stereo_monitor2';
    case 6
        sequenceName = 'observation_rpg_DAVIS_stereo_recycler';
    case 7
        sequenceName = 'observation_indoor_flying3_transfer';
    case 8
        sequenceName = 'observation_indoor_flying1_transfer';
    otherwise
        disp('The selected sequence does not exist');
end

% reconstruction objects set by numRV
load([Observation_dir, sequenceName, '.mat']);
disp('loaded data');

% %% To generate the configuration figure in the paper.
% addpath('/home/zhouyi-anu/workspace/project/std_matlab_plot/altmany-export_fig-83ee7fd');
% saveDir = '/home/zhouyi-anu/workspace/paper/pic/configuration';
% for i = 100:500
% figure(1), imshow(EMaps_left{i});
% export_fig('Rectified_EMap_left','-png','-transparent');
% movefile('Rectified_EMap_left.png', saveDir);
% 
% figure(2), imshow(EMaps_right{i});
% export_fig('Rectified_EMap_right','-png','-transparent');
% movefile('Rectified_EMap_right.png', saveDir);
% end
%%
% figure;
plotRow = 4;
plotCol = 2;
bShowEpipolarLine = false;
xx = 1:width;
yy1 = 40 * ones(1,length(xx));
yy2 = 80 * ones(1,length(xx));
yy3 = 120 * ones(1,length(xx));
yy4 = 160 * ones(1,length(xx));

for i = 1:100
    %
    subplot(plotRow,plotCol,1);
    imshow(EMaps_left{i});
    title('events in the left DVS');
    %
    subplot(plotRow,plotCol,2);
    imshow(EMaps_right{i});
    title('events in the right DVS');
    %
    subplot(plotRow,plotCol,3);
    imagesc(SAEs_smooth_left{i},'CDataMapping','scaled');
    hold on;
    if bShowEpipolarLine
    plot(xx,yy1,'r-');
    plot(xx,yy2,'r-');
    plot(xx,yy3,'r-');
    plot(xx,yy4,'r-');
    end
    title('SAE in the left DVS');
    
    %
    subplot(plotRow,plotCol,4);
    imagesc(SAEs_smooth_right{i},'CDataMapping','scaled');
    hold on;
    if bShowEpipolarLine
    plot(xx,yy1,'r-');
    plot(xx,yy2,'r-');
    plot(xx,yy3,'r-');
    plot(xx,yy4,'r-');
    end
    title('SAE in the right DVS');
    
%     subplot(plotRow,plotCol,5);
%     rv_SAE = SAEs_smooth_left{i};
%     age_threshold = prctile(rv_SAE(rv_SAE(:) > 0),50);
%     rv_SAE(rv_SAE(:) < age_threshold) = 0;
%     rv_SAE = cv.erode(uint8(rv_SAE));
%     imshow(rv_SAE);
%     title('Eroded SAE in the left DVS');
%     
%     subplot(plotRow,plotCol,6);
%     rv_SAE = imsharpen(rv_SAE,'Radius',2,'Amount',3);
%     bw = imbinarize(rv_SAE,0.7);
%     imshow(bw);
%     title('Sharpened and binarized SAE in the left DVS');
%     
%     %
    subplot(4,2,5);
    imagesc(SAEs_positive_smooth_left{i},'CDataMapping','scaled');
    hold on;
    if bShowEpipolarLine
    plot(xx,yy1,'r-');
    plot(xx,yy2,'r-');
    plot(xx,yy3,'r-');
    plot(xx,yy4,'r-');
    end
    title('Positive SAE in the left DVS');
%     
%     %
    subplot(4,2,6);
    imagesc(SAEs_positive_smooth_right{i},'CDataMapping','scaled');
    hold on;
    if bShowEpipolarLine
    plot(xx,yy1,'r-');
    plot(xx,yy2,'r-');
    plot(xx,yy3,'r-');
    plot(xx,yy4,'r-');
    end
    title('Positive SAE in the right DVS');
%     
%     %
    subplot(4,2,7);
    imagesc(SAEs_negative_smooth_left{i},'CDataMapping','scaled');
    hold on;
    if bShowEpipolarLine
    plot(xx,yy1,'r-');
    plot(xx,yy2,'r-');
    plot(xx,yy3,'r-');
    plot(xx,yy4,'r-');
    end
    title('Negative SAE in the left DVS');
%     
%     %
    subplot(4,2,8);
    imagesc(SAEs_negative_smooth_right{i},'CDataMapping','scaled');
    hold on;
    if bShowEpipolarLine
    plot(xx,yy1,'r-');
    plot(xx,yy2,'r-');
    plot(xx,yy3,'r-');
    plot(xx,yy4,'r-');
    end
    title('Negative SAE in the right DVS');
    
    disp(i);
%     pause(0.1);
end