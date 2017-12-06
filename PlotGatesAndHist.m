function [ h ] = PlotGatesAndHist(BgDataAll, GateArray1, GateArray2, x, y, z)
% [ h ] = PlotGatesAndHist(BgDataAll,GateArray1, GateArray2,x, y, z)
% Hace una figura de scatter plot con los Gates en x&y y luego unos
% histogramas de z separando a x de y
%
% Funciona con los fcs ya cargados de LoadCytFiles.
% Dependiendo de tus datos, tal vez te convenga cambiar el xlim en este archivo.
%
% Funciones que necesita:
% GatingMouse HistLine->smooth

for i=1:length(BgDataAll)%6:8%1:length(BgDataAll)
    figure(i)
    a=size(BgDataAll(i).dat);
    if a(2)==10

        [NewData gated] = GatingMouse(BgDataAll(i).dat, x,y,BgDataAll(1).info(x).Name,BgDataAll(1).info(y).Name, 'log','log', 0,GateArray1);
        [NewData ungated] = GatingMouse(BgDataAll(i).dat, x,y,BgDataAll(1).info(x).Name,BgDataAll(1).info(y).Name, 'log','log',0, GateArray2);
        
        tam=min(length(gated), length(ungated));%es para que sean del mismo tamaño los dos histogramas
        
        subplot(1,2,1) %SCATTER PLOT DE X/Y CON SUS GATES, ESCALA LOGARITMICA
        plot(BgDataAll(i).dat(:,x),BgDataAll(i).dat(:,y),'ok
        ','MarkerSize',1)
        hold on
        plot(BgDataAll(i).dat(gated,x),BgDataAll(i).dat(gated,y),'og','MarkerSize',1)
        plot(BgDataAll(i).dat(ungated,x),BgDataAll(i).dat(ungated,y),'ob','MarkerSize',1)
        set(gca, 'yscale','log','xscale','log')
        xlabel(BgDataAll(i).info(x).Name)
        ylabel(BgDataAll(i).info(y).Name)
        title(strcat(BgDataAll(i).filename,'-',mat2str(i)))
        plot(GateArray1(:,1), GateArray1(:,2), 'g-', 'linewidth', 1)
        plot(GateArray2(:,1), GateArray2(:,2), 'b-', 'linewidth', 1)

        meangated=mean(log2(BgDataAll(i).dat(gated,z)));
        meangated2=mean(log2(BgDataAll(i).dat(ungated,z)));
        size(BgDataAll(i).dat);
        size(gated);
        meangated-meangated2;
        difscores(i) = meangated-meangated2;
        medias(i) = meangated;

        [heigth1 center1]=hist(log2(BgDataAll(i).dat(gated(1:tam),z)),100);
        [heigth2 center2]= hist(log2(BgDataAll(i).dat(ungated(1:tam),z)),100);
        
        subplot(1,2,2) % Histogramas 
        hold on
        plot([meangated meangated],[-1 2000],'g')
        plot([meangated2 meangated2],[-1 2000],'b')
        xlim([0 8])
        ylim([0 max([heigth1 heigth2])+30 ])
         HistLine(center1, heigth1,'g')
        title('Gated')
        
        hold on
        xlabel(BgDataAll(i).info(z).Name)
        plot([meangated meangated],[-1 2000],'g')
        plot([meangated2 meangated2],[-1 2000],'b')
        title('UnGated')
        xlabel(BgDataAll(i).info(z).Name)
        HistLine(center2, heigth2, 'b')
        xlim([0 8])
        ylim([0 max([heigth1 heigth2])+30])        
    end
end

end