function handle = plot_CM(CM,titleString,Ncolors,zerocolor)
% Plots a correlation matrix CM with color options.
%
% Function written by Jason E. Hill, Ph.D., modified from code provided by 
% Elena Allen in association with Vince Calhoun's research group for use in
% their "toy model" FCS simulator
%
% Date: 24 JAN 2017

if nargin < 4
    zerocolor = [1; 1; 1];
else
    if ischar(zerocolor)
        zerocolor = squeeze(double(label2rgb(0, jet, zerocolor)))/255;
    end
end
if nargin < 3
    Ncolors = 256;
end
if nargin < 2
    titleString = [];
end

nR = size(CM,1);
if size(CM,2) ~= nR
    error('Correlation Matrix should be square!');
end
    
CMdiag0 = CM;
for r = 1:nR
    CMdiag0(r,r) = 0.0;
end

CMmax = max(CMdiag0(:));
CMmin = min(CMdiag0(:));

CMmap = make_CM_colormap(Ncolors,zerocolor);

% plot correlation matrix using imagesc()
figure
handle = imagesc(1:nR,1:nR,CM);
colormap(CMmap)
climit = max(abs(CMmin),abs(CMmax));
caxis([-climit,climit])
colorbar
set(gca, 'ticklength', [0 0])
axis xy
hold on;
xlimits = get(gca, 'xlim');
ylimits = get(gca, 'ylim');
xgridlines = xlimits(1):xlimits(2); % where to draw grid lines x 
ygridlines = ylimits(1):ylimits(2); % where to draw grid lines x 
plot(repmat(xgridlines,2,1), repmat(ylimits',1,length(xgridlines)), 'k'); % vertical lines
plot(repmat(xlimits',1,length(ygridlines)), repmat(ygridlines,2,1), 'k'); % horizontal lines
set(gca,'XTick',1:nR);
set(gca,'XTickLabel',1:nR);
set(gca,'YTick',1:nR);
set(gca,'YTickLabel',1:nR);
axis square; 
axis ij
set(gca, 'XTick', [], 'YTick', [])%, 'XColor', [1 1 1], 'YColor', [1 1 1])
child = get(gca, 'Children');
set(child(strcmp(get(child, 'Type'),'line')), 'Color', 'w');
if ~isempty(titleString)
title(titleString)
end
    
