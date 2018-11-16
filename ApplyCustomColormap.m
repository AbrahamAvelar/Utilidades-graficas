function ApplyCustomColormap(data, fixedToWhite, bins, symetric, nans )
% ApplyCustomColormap(data, fixedToWhite)
% Adaptado de internet, de uno que hizo Brandon Eidson en Matlab Answers
% bins es por cuanto multiplicar el numero de colores que hay entre el min
% y el maximo, pero ya toma en cuenta cuantos datos hay. para todo fin
% pr'actico conviene ponerle 1 o 10 dependiendo de cuantos datos tengas

if nargin<4
    nans=0;
    symetric=0;
end

L = length(data);             %number of datapoints
%data = 3.6*rand(L); % create example data set with values ranging from 0 to 3.6

indexValue = fixedToWhite;     % value for which to set a particular color

topColor = [1 0 0];         % color for maximum data value (red = [1 0 0])
indexColor = [0 0 0];       % color for indexed data value (white = [1 1 1])
bottomcolor = [0 1 0];      % color for minimum data value (blue = [0 0 1])

% Calculate where proportionally indexValue lies between minimum and
% maximum values

largest = max(max(data));
smallest = min(min(data));
index = L*abs(indexValue-smallest)/(largest-smallest);

% Create color map ranging from bottom color to index color
% Multipling number of points by 100 adds more resolution
customCMap1 = [linspace(bottomcolor(1),indexColor(1),bins*index)',...
            linspace(bottomcolor(2),indexColor(2),bins*index)',...
            linspace(bottomcolor(3),indexColor(3),bins*index)'];

% Create color map ranging from index color to top color
% Multipling number of points by 100 adds more resolution
customCMap2 = [linspace(indexColor(1),topColor(1),bins*(L-index))',...
            linspace(indexColor(2),topColor(2),bins*(L-index))',...
            linspace(indexColor(3),topColor(3),bins*(L-index))'];

% if symetric
%     menordif=min(abs(largest), abs(smallest));
%     smallest = indexValue - menordif;
%     largest  = indexValue + menordif;
% end

        
        
        
customCMap = [customCMap1;customCMap2];  % Combine colormaps


if nans
    customCMap(1,:) = [.5 .5 .5];
end


colormap(customCMap)
%psudo = pcolor(data);
colorbar

end