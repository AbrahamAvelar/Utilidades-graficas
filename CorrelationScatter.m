function [plotHandle r p n m] = CorrelationScatter(x,y,labelx, labely, namesx, namesy, textos)
% [plotHandle r p] = CorrelationScatter(x,y,labelx, labely, namesx, namesy,textos)
error(nargchk(0,7,nargin))
if length(x) ~= length(y)
    warning('X and Y have different sizes');
end

if nargin>4 %if names are provided, they must be paired in both vectors
    i=FindGeneList(namesx,namesy);
    y=y(i);
end

xvalid = find(~isnan(x));
yvalid = find(~isnan(y));

validas = intersect(xvalid,yvalid);
[r p]=corrcoef(x(validas), y(validas));
r=r(2);
p=p(2);

plotHandle = plot(x(validas),y(validas),'ok');
%legend(strcat(' r= ',num2str(r),' p= ',num2str(p)),'location','best')
titulo=strcat(' r= ',num2str(r),' p= ',num2str(p), '  n=', num2str(length(validas)));
FastLabels(labelx, labely, titulo);

n=length(validas);
m=robustfit(x(validas),y(validas));
if nargin>6
    text(x(validas)+.001,y(validas),textos(validas))
end
grid on
axis square
end

