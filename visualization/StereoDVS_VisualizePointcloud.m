function StereoDVS_VisualizePointcloud(pointcloud, VisRange, ViewAngle, filename, saveDir)
% PLEASE ADD PATH TO WHERE EXPORT_FIG PACKAGE IS INSTALLED IF THIS SCRIPT
% CANNOT BE EXECUTED !!!
% https://ch.mathworks.com/matlabcentral/fileexchange/23629-export-fig
% addpath('PATH TO altmany-export_fig-83ee7fd');

% visualize pointcloud.
pointcloud(:,pointcloud(3,:) > VisRange(2)) = [];
pointcloud(:,pointcloud(3,:) < VisRange(1)) = [];
min_d = max(VisRange(1), min(pointcloud(3,:)));
max_d = min(max(pointcloud(3,:)), VisRange(2));

color_index = floor( (pointcloud(3,:)-min_d) / (max_d - min_d) * 63) + 1;
colorMap = flipud(colormap('jet'));
color = colorMap(color_index,:);
figure, pcshow(pointcloud', color, 'MarkerSize', 10), view(ViewAngle);
% figure, pcshow(pointcloud', 'MarkerSize', 10), view([0 0 -1]);
% colorbar;
% caxis([0.0 6]);
hold on;

% draw camera
% ax = 0:0.02:1;
% camX_axis = [ax;zeros(size(ax));zeros(size(ax))];
% camY_axis = [zeros(size(ax));ax;zeros(size(ax))];
% camZ_axis = [zeros(size(ax));zeros(size(ax));ax];
% plot3(camX_axis(1,:),camX_axis(2,:),camX_axis(3,:), 'r-');
% plot3(camY_axis(1,:),camY_axis(2,:),camY_axis(3,:), 'g-');
% plot3(camZ_axis(1,:),camZ_axis(2,:),camZ_axis(3,:), 'b-');


% xlim([-5,5]);
% ylim([-3,3]);
% zlim([0.3,6]);
axis off;

if ~isempty(filename) && ~isempty(saveDir)
export_fig('pointcloud','-png','-transparent'); 
movefile('pointcloud.png', saveDir);
% save
savefig(filename);
movefile(filename, saveDir);
end
end

