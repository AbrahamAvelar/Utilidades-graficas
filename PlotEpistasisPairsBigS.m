function PlotEpistasisPairsBigS( Muts, namesx, namesy, coeffsBigS, coeffsicSdw, coeffsicSup, mediaX, errXup, errXdw, mediaY, errYup, errYdw )
% PlotEpistasisPairsBigS( MutQuery, MutTarge, namesx, namesy, coeffsBigS,..
% coeffsicSdw, coeffsicSup, mediaX, errXup, errXdw, mediaY, errYup, errYdw)
%
% Basado en PlotEpistasisPairs y PlotSEpistasis
    MQ=FindGeneList(Muts, namesx);
    MT=FindGeneList(Muts, namesy);
    tamx=length(MQ);%tamano de x
    tamy=length(MT);
    %con=0;
    %figure()
    for PL=1:tamx
        X = mediaX(MQ(PL));
        Xup = errXup(MQ(PL));
        Xdw = errXdw(MQ(PL));
      for W=PL+1:tamy
      %    con=con+1;
      %    subplot(tamx,tamy,con)
      %      figure()
            Y   = mediaY(MQ(W));
            Yup = errYup(MQ(W));
            Ydw = errYdw(MQ(W));
            XY  = coeffsBigS ( MQ(PL),MT(W) );
            XYup= coeffsicSup( MQ(PL),MT(W) );
            XYdw= coeffsicSdw( MQ(PL),MT(W) );
            YX  = coeffsBigS ( MQ(W),MT(PL) );
            YXup= coeffsicSup( MQ(W),MT(PL) );
            YXdw= coeffsicSdw( MQ(W),MT(PL) );
        
            ejex = [1 2];% 2.8 3.2];
            ejey = 1+[X Y];% XY YX];
            yerrTop= 1+[Xup Yup];% XYup YXup]+1;
            yerrBott= 1+[Xdw Ydw];% XYdw YXdw]+1;
            EXPmul = (1+X)*(1+Y);
            PlotConError( ejex,ejey,yerrBott,yerrTop, [.3 .9 .3],'o',3,.25 );
            plot(3, mean(1+[XY YX]),'o','color',[.3 .9 .3],'markersize',3,'Markerfacecolor',[.3 .9 .3])
            if YX+1<EXPmul
                PlotConError( ejex,ejey,yerrBott,yerrTop, [.9 .3 .3],'o',3,.25 );
                plot(3, mean(1+[XY YX]),'o','color',[.9 .3 .3],'markersize',3,'Markerfacecolor',[.9 .3 .3] )
            end
            PlotConError([3.05 2.95], 1+[XY YX],[XYdw YXdw]+1,[XYup YXup]+1,  [.6 .6 .6], '.',3, .005);
            
            hold on
            plot([0 4],[1 1],':','color',[.5 .5 .5])
            plot([0 4],[EXPmul EXPmul], '--','color',[.5 .5 1], 'linewidth', 2 );
            ylim([0.9 1.06])
            
            text( (1:2)-.38, ones(2,1)*.89, [namesx(MQ(PL)),namesy(MT(W))])%,'rotation',270 );
            text( 3-.69, .89, strcat( namesx(MQ(PL)),'-',namesy(MT(W))) )%,'rotation',270 );
            
            %text( (1:3)-.03, ones(3,1)*.899, [namesx(MQ(PL)),namesy(MT(W)),strcat( namesx(MQ(PL)),'-',namesy(MT(W)))] ,'rotation',270 );
            
            
%            set(gca,'xtick', 1:3,'xticklabel', [namesx(MQ(PL)),namesy(MT(W))  strcat( namesx(MQ(PL)),'/',namesy(MT(W))) ] )
            xlim([.5 3.5])
            box off
       end
    end
end