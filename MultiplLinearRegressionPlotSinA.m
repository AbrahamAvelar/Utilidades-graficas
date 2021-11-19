function MultiplLinearRegressionPlotSinA( EXPbgdataScriptEGG,cepas,cepa,refs,colores,tfijo)%
% handle = MultiplLinearRegressionPlot(EXPbgdataScriptEGG,cepas,cepa,refs,colores)
%
% MultiplLinearRegressionPlot( EXPbgdataScriptEGG,cepas,cepa,refs)
% EXPbgdataScriptEGG es un plato en específico
% cepa es el hash con los nombres de las mutantes
% strains = el índice de la réplica que quier


%subplot(1,2,1)
intervalos = [find(EXPbgdataScriptEGG.tOut==0) length(EXPbgdataScriptEGG.tOut)+1];
nombres=fieldnames(cepas);
for i = cepa %[cepas.refComp cepas.hap3]% 1:96
    w = cepas.(str2mat(nombres(i)));
    LOGRoC = log(EXPbgdataScriptEGG.mut(:,w)./EXPbgdataScriptEGG.ref(:,w)) ;
    diaCs=0;
    toplotYfixt=[]; toplotXfixt=[];
    for dia = 1:length(EXPbgdataScriptEGG.Tdays)-1
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
            plot( [T T], [yfitnotGBOTCI yfitnotGTOPCI]-A,'color',[.5 .5 .5] ) ;
            hold on
        end
        hora=t*24;
        plot(t+T, RoC-A,'.','color', colores(i,:) ,'linewidth',.5,'MarkerFaceColor', colores(i,:) );% ) [1/16*T 1/16*T 1/16*T]
        hold on
        plot([T  T+t],[ynoG yfit']-A,'-','color',colores(i,:)); %PlotConError(T+t,yfit,yfitTOPCI, yfitBOTCI,colores(i,:),'s-',3,.0051);%Misma Cosa, con o sin error
        plot(T+t,yfit-A,'x','color',colores(i,:),'MarkerSize',5,'markerfacecolor',colores(i,:)); %PlotConError(T+t,yfit,yfitTOPCI, yfitBOTCI,colores(i,:),'s-',3,.0051);%Misma Cosa, con o sin error
        toplotXfixt = [toplotXfixt T];
        toplotYfixt = [toplotYfixt ynoG];
    end
    plot(toplotXfixt, toplotYfixt-A,'o-','markersize',5,'color',colores(i,:))
    As(i)=A;
end
%     plot( [-0.015,-0.015], [0 As(1)],'-','color',colores(i,:))
%     plot( [-0.03,0.0], [As(1) As(1)],'-','color',colores(i,:))
%     plot( [-0.03,0.0], [0 0],'-','color',colores(i,:))
%     plot( [-0.03,0.0], [0 0],'-','color',colores(i,:))
xlim([-1 16])
%ylim([-3 1.2])
title('ln(mutFP/refFP) = A_x_w_i+S_w_i*T_(_1_._._n_)+G_w_i*t_(_1_._._m_)+C_T_(_1_._._n_)_,_t_(_1_._._m_)' );
xlabel('T, days')
ylabel('ln(RFP/CFP)')
% set(gca, 'ytick',-3:1:1,'xtick', 0:3:15, 'linewidth', 1)

