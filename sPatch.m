function h = sPatch(s, days, col, CI)
% h = sPatch(measurements, days)
% s -> is a >2 vector with slopes, s
% days -> is a 2 element vector with the initial and final days to be ploted
% in x
% color -> is the color of the line

if nargin<4
    CI=0;
end

smean = nanmean(s);
sstd = nanstd(s);

hold on
X =[days(1) days(1) days(2), days(2)];

if CI
    n = sum(not(isnan(s)));
    sSEM = sstd./sqrt(n);
    CI =sSEM*1.96;
    Y =[days(1)*(smean-CI) days(1)*(smean+CI) days(2)*(smean-CI) days(2)*(smean+CI)];
else
    Y =[days(1)*(smean-sstd) days(1)*(smean+sstd) days(2)*(smean-sstd) days(2)*(smean+sstd)];
end
h(1) = patch( X, Y, [.73 .73 1],'edgeColor', [.78,.78, 1], 'LineWidth',1 );
h(2) = plot( days , [days(1)*smean; days(2)*smean], '-', 'color', col, 'linewidth', 1.5 );

end
