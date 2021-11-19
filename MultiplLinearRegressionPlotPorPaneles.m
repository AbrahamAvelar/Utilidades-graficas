function MultiplLinearRegressionPlotPorPaneles(EXPbgdataScriptEGG,cepas,cepa,refs,colores,tfijo,subplotNo,subplotind)%
% handle =MultiplLinearRegressionPlot(EXPbgdataScriptEGG,cepas,cepa,refs,colores,tfijo,subplotNo)%
%
% MultiplLinearRegressionPlot( EXPbgdataScriptEGG,cepas,cepa,refs)
% EXPbgdataScriptEGG es un plato en específico
% cepa es el hash con los nombres de las mutantes
% strains = el índice de la réplica que quier

% subplot(1,2,1)

intervalos = [find(EXPbgdataScriptEGG.tOut==0) length(EXPbgdataScriptEGG.tOut)+1];
nombres=fieldnames(cepas);
for i = cepa %[cepas.refComp cepas.hap3]%1:96
    w = cepas.(str2mat(nombres(i)));
    LOGRoC = log(EXPbgdataScriptEGG.mut(:,w)./EXPbgdataScriptEGG.ref(:,w)) ;
    diaCs=0;
    toplotYfixt=[]; toplotXfixt=[];
    for dia = 1:length(EXPbgdataScriptEGG.Tdays)-1
        if subplotind==1
            subplot( subplotNo,length(EXPbgdataScriptEGG.Tdays)-1, dia )
        else
        	subplot( subplotNo,length(EXPbgdataScriptEGG.Tdays)-1, (subplotind-1)*(length(EXPbgdataScriptEGG.Tdays)-1)+dia )
        end
        start=intervalos(dia)+1;
        finish=intervalos(dia+1)-1;
        t = EXPbgdataScriptEGG.tOut(start:finish);
        RoC = LOGRoC(start:finish);
        T = EXPbgdataScriptEGG.Tdays(dia);
        S = EXPbgdataScriptEGG.S(w)  ;
        A = EXPbgdataScriptEGG.A(w)  ;
        G = EXPbgdataScriptEGG.G(w)  ;
        icStop=EXPbgdataScriptEGG.icS(1,w);%intervalos de confianza
        icSbot=EXPbgdataScriptEGG.icS(2,w);
        if length(diaCs)<=dia;
            diaCs = [diaCs diaCs(dia)+length(t)];
            C=EXPbgdataScriptEGG.C(diaCs(dia)+1:diaCs(dia+1));
        end
        unos=ones(length(t),1);
        if ismember(w,refs);
            yfit=unos*A;
            yfitfixedt=A;
            ynoG=A;
            ysoloS=0;
            yfitTOPCI=A;
            yfitBOTCI=A;
        else
            yfit =unos*A + unos*(S*T) + (G*t*24)';%CHECAR EN QUï¿½ UNIDADES ESTï¿½ LA G DE NOAM DIAS? U HORAS (*24)?
            yfitfixedt=A +(S*T)+G/tfijo*24; %+ G/3;% 1/3 en t= 8 horas de tiempo fijo
            yfitTOPCI = unos*A + unos*(icStop*T) + (G*t*24)'; %
            yfitBOTCI = unos*A + unos*(icSbot*T) + (G*t*24)'; %
            ynoG=A+S*T;
            ysoloS=S*T;
            yfitnotGTOPCI = A + (icStop*T) ; %
            yfitnotGBOTCI = A + (icSbot*T) ; %
            plot( [0 0], [yfitnotGBOTCI yfitnotGTOPCI],'color',[.5 .5 .5] ) ;
            hold on
        end
    	hora=t*24;
         for k=1:length(t)
             %plot([t(k)+T;t(k)+T]*24, [yfit(k)';yfit(k)'+C(k)'], '-','color',colorin(floor(hora(k)*salto),:), 'color',[.5 .5 .5])%colorin(floor(hora(k)*salto),:), 'MarkerFaceColor', [1/16*T 1/16*T 1/16*T] );
             plot([t(k);t(k)], [yfit(k)';yfit(k)'+C(k)'], '-', 'color',[.52 .52 .52])
             hold on
         end
        plot(t, RoC,'.','color', colores(i,:) ,'linewidth',.5,'MarkerFaceColor', colores(i,:) );% ) [1/16*T 1/16*T 1/16*T]
        hold on
        plot([0  t],[ynoG yfit'],'-','color',colores(i,:)); %PlotConError(T+t,yfit,yfitTOPCI, yfitBOTCI,colores(i,:),'s-',3,.0051);%Misma Cosa, con o sin error
        plot(t,yfit,'x','color',colores(i,:),'MarkerSize',5,'markerfacecolor',colores(i,:)); %PlotConError(T+t,yfit,yfitTOPCI, yfitBOTCI,colores(i,:),'s-',3,.0051);%Misma Cosa, con o sin error
        toplotXfixt = [toplotXfixt T];
        toplotYfixt = [toplotYfixt ynoG];
        ylim([-2 1.25])
        xlim([-.1 .7])
        title( strcat('T_',num2str(dia),'=',num2str(floor(T*100)/100) ))
        set(gca, 'xtick', 0:1/3:1,'xticklabel', (0:1/3:1)*24 )
        if dia ~= 1
            set(gca,'yticklabel',[])
        end
        box off
    end
%     plot(toplotXfixt, toplotYfixt,'o-','markersize',5,'color',colores(i,:))
%     As(i)=A;
end
%     plot( [-0.015,-0.015], [0 As(i)],'-','color',colores(i,:))
%     plot( [-0.03,0.0], [As(i) As(i)],'-','color',colores(i,:))
%     plot( [-0.03,0.0], [0 0],'-','color',colores(i,:))
%     plot( [-0.03,0.0], [0 0],'-','color',colores(i,:))
%xlim([-1 16])
%ylim([-3 1.5])
%title('ln(mutFP/refFP) = A_x_w_i+S_w_i*T_(_1_._._n_)+G_w_i*t_(_1_._._m_)+C_T_(_1_._._n_)_,_t_(_1_._._m_)' );
xlabel('t, hours')
%ylabel('ln(RFP/CFP)')
% set(gca, 'ytick',-3:1:1,'xtick', 0:3:15, 'linewidth', 1)

