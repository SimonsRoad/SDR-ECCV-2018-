function StereoDVS_VisualizeDepthMap(invDepthMap, DepthRange, filename, saveDir)
 % PLEASE ADD PATH TO WHERE EXPORT_FIG PACKAGE IS INSTALLED IF THIS SCRIPT
% CANNOT BE EXECUTED !!!
% https://ch.mathworks.com/matlabcentral/fileexchange/23629-export-fig
% addpath('PATH TO altmany-export_fig-83ee7fd');

% visualize DepthMap, uncertainty map, remaining cost map.
[height, width, channel] = size(invDepthMap);
rhoMap  = invDepthMap(:,:,1);
dMap = 1./rhoMap;
dMap(dMap(:) > DepthRange(2)) = DepthRange(1);
dMap(dMap(:) < DepthRange(1)) = DepthRange(1);
dMap(isnan(dMap(:))) = DepthRange(1);
varMap = invDepthMap(:,:,2);

% draw rhoMap
min_d = DepthRange(1);
max_d = DepthRange(2);
color_index = floor( (dMap(:) - min_d) / (max_d - min_d) * 63) + 1;
colorMap = flipud(colormap(jet));
colorMap(1,:) = [0 0 0];
colorDMap = vec2mat(color_index, height)';
figure, imagesc(colorDMap);
colormap(colorMap);
axis off;
% colorbar;
hold on;
title('InvDepth Map');
% save depth map
file1 = [filename,'_depthMap'];
if ~isempty(filename) && ~isempty(saveDir)
export_fig(file1,'-png','-transparent', '-native'); movefile([file1,'.png'], saveDir);
end

end