function h = JitterPlot (thisX, y, jitter, shape, col, MarkSize)
% h = JitterPlot (thisX, y, jitter, shape, col)
%   thisX el centro donde quieres que esten dispersos los puntos
%   los valores en el eje y
%   jitter es cuanto quieres que se dispersen los valores en x
%   shape puede ser 'o', '.-', ':', etc.
%   col es el color, puede ser 'g','y','b',[.1 .2 1], [R G B], etc.

if nargin<4
    shape='.';
    col='k';
elseif nargin<6
    MarkSize = 5;
end

        J=(rand(size(y))-0.5)*jitter;
        h = plot( thisX+J, y, shape, 'color', col, 'markersize', MarkSize );


end
