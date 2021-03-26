function h = JitterPlot (thisX,y,  jitter, shape, col, markersize, media )
% h = JitterPlot (thisX,y,  jitter, shape, col, markersize, media )
%   thisX el centro donde quieres que esten dispersos los puntos
%   los valores en el eje y
%   jitter es cuanto quieres que se dispersen los valores en x
%   shape puede ser 'o', '.-', ':', etc.
%   col es el color, puede ser 'g','y','b',[.1 .2 1], [R G B], etc.
%  media = 1 pone la nanmean
%  media = 2 pone la nanmedian

if nargin<4
    shape='.';
    col='k';
end
if nargin<6
    markersize=5;
end
        
if nargin<7
    media=0;
end

if media == 1
     plot([thisX-jitter thisX+jitter], [nanmean(y(:)) nanmean(y(:))], '-', 'color', col*.8, 'linewidth', 2 )
     hold on
elseif media == 2
	plot([thisX-jitter thisX+jitter], [nanmedian(y(:)) nanmedian(y(:))], '-', 'color', col*.8, 'linewidth', 2 )
     hold on
end
        J=(rand(size(y))-0.5)*jitter;
        hold on
        h = plot( thisX+J, y, shape, 'color', col, 'markersize', markersize, 'markerfacecolor', col );
        hold on


end
