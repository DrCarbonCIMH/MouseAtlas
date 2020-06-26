function Plot_cormat(T,names,lim1, lim2)
% plot a correlation matrix
% T = correlation matrix calculated by Calc_Cormat.m
% names = names of ROIs
% lim1/lim2 = lower/upper limit of color range (optional)

cmap=colormap('jet');
if nargin <3 % no color limits set
    Tmax = max(max(abs(T)));
    lim1 = -Tmax;
    lim2 = Tmax;
end
% reshape the names to also include numbering
for n=1:length(names)
    names{n}=[num2str(n) ': ' strrep(names{n}, '_', '\_')];
    xnames{n}=num2str(n);
end

% setup figure for plotting
figure; imagesc(T,[lim1 lim2]); colormap(cmap); colorbar; axis image;
set(gca,'XTick',[0.5:1:size(T,1)],'XTickLabel',[xnames],'YTick',[0.5:1:size(T,1)],'YTickLabel',names);
grid on

for x=1:size(T,1)
    for y=x:size(T,2)
        if (x == y)
            xv=[x- 0.5 x-0.5 x+.5 x+.5];yv=[y-.5 y+.5 y+.5 y-.5];
            patch(xv,yv,[1 1 1])
        end
    end
end
