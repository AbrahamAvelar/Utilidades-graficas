function h=MultiplLinearRegressionPlot( OutputModelASGC,pozo,isref,color,PrintValues,tfijo)%
% 
% h=MultiplLinearRegressionPlot( OutputModelASGC,pozo,isref,color,tfijo)%
% OutputModelASGC es un plato en específico del PRIMER Output de ModelASGC,
% contiene A, S, G, C, icS, icG, etc.
% pozo es la posición del plato que se quiere graficar
% isref = 1 si es una competencia entre referencias. 0 si es con una mutante
% color es un vector de 3 números entre 0 y 1 RGB, por ejemplo [.5 .2 1]
% PrintValues = 1 imprime cuánto vale A, S, G y C
% PrintValues = 2 imprime cuánto vale S
% PrintValues = 0 no imprime valores
% tfijo es el factor de las unidades en las que se calculó la G respecto a
% los días, por ejemplo si se calculó en horas tfijo=24, si se hizo en días
% tfijo=1.

%subplot(1,2,1)
intervalos = [find(OutputModelASGC.tOut==0) length(OutputModelASGC.tOut)+1];
%nombres=fieldnames(cepas);
for i = pozo %[cepas.refComp cepas.hap3]% 1:96
    w = i; % cepas.(str2mat(nombres(i)));
    LOGRoC = log(OutputModelASGC.mut(:,w)./OutputModelASGC.ref(:,w)) ;
    diaCs=0;
    toplotYfixt=[]; toplotXfixt=[];
    for dia = 1:length(OutputModelASGC.Tdays)-1
        start=intervalos(dia)+1;
        finish=intervalos(dia+1)-1;
        t = OutputModelASGC.tOut(start:finish);
        RoC = LOGRoC(start:finish);
        T = OutputModelASGC.Tdays(dia);
        S = OutputModelASGC.S(w)  ;
        A = OutputModelASGC.A(w)  ;
        G = OutputModelASGC.G(w)  ;
        icStop=OutputModelASGC.icS(1,w);%intervalos de confianza
        icSbot=OutputModelASGC.icS(2,w);
        if length(diaCs)<=dia;
            diaCs = [diaCs diaCs(dia)+length(t)];
            C=OutputModelASGC.C(diaCs(dia)+1:diaCs(dia+1));
        end
        unos=ones(length(t),1);
        if isref %ismember(w,refs);
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
            plot( [T T], [yfitnotGBOTCI yfitnotGTOPCI],'color',[.5 .5 .5] ) ;
            hold on
        end
        hora=t*24;
        plot(t+T, RoC,'.','color', color ,'linewidth',.5,'MarkerFaceColor', color );% ) [1/16*T 1/16*T 1/16*T]
        hold on
        plot([T  T+t],[ynoG yfit'],'-','color',color); %PlotConError(T+t,yfit,yfitTOPCI, yfitBOTCI,color,'s-',3,.0051);%Misma Cosa, con o sin error
        plot(T+t,yfit,'x','color',color,'MarkerSize',5,'markerfacecolor',color); %PlotConError(T+t,yfit,yfitTOPCI, yfitBOTCI,color,'s-',3,.0051);%Misma Cosa, con o sin error
        toplotXfixt = [toplotXfixt T];
        toplotYfixt = [toplotYfixt ynoG];
    end
    plot(toplotXfixt, toplotYfixt,'o-','markersize',5,'color',color)
    As(i)=A;
end
    plot( [-0.015,-0.015], [0 As(1)],'-','color',color)
    plot( [-0.03,0.0], [As(1) As(1)],'-','color',color)
    plot( [-0.03,0.0], [0 0],'-','color',color)
    h=plot( [-0.03,0.0], [0 0],'-','color',color);
xlim([-1 19])
%ylim([-3 1.2])
title('ln(mutFP/refFP) = A_x_w_i+S_w_i*T_(_1_._._n_)+G_w_i*t_(_1_._._m_)+C_T_(_1_._._n_)_,_t_(_1_._._m_)' );
%text(0,0,strcat('ln(mutFP/refFP) = ',num2str(A),'+',num2str(S), '*T_(_1_._._n_)+', num2str(G),'*t_(_1_._._m_)+C_T_(_1_._._n_)_,_t_(_1_._._m_)') );
if PrintValues == 1
    text(0,-2,strcat('A=',num2str(A),' S=',num2str(S),' G=',num2str(G)) );
elseif PrintValues ==2
    text(0,0,strcat(' S=',num2str(S)));
end
xlabel('T, days')
ylabel('ln(RFP/CFP)')
% set(gca, 'ytick',-3:1:1,'xtick', 0:3:15, 'linewidth', 1)
