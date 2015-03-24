function parsave_storm(savefname,consumerstats,vacstats,BUILDTIME,VACLAND,...
    RENT,RETURN,LOTTYPE,BASELAYER,Rpop,Rvacrate,Rvaclots,numlt,Rleftoverpop,...
    avgrentdynms,rentdynstats,farmsoldinfo,avgrent,avgfarmsize,stdfarmsize,...
    DELTA,survivalrate,LOCWGHT,REGWGHT,PCTSEARCH,zeta,HIBETA,MIDBETA,LOWBETA,...
    POPGROW,ccost,newhouseset,bidtot,meandisp,maxdevdist,setupmap,zonedensity,...
    vacrate,Farminfo,oldincome,Realreturn,Realavgret,Exptrentdiff,...
    Avgexptdiff,htincome,numtotbids,totfarmrecord,htperyear,Newbidlevel,...
    totbrokerrecord,Farmdist2dev,Bprojmap,Bmodelmap,WTAlandmap,WTPlandmap,...
    Landprojmap,Landmodmap,bidshare,buildshare,landdemand,EXPTHOUSE,lotchoice,...
    Exptprofit,Exptret,Realexptret,Realavgexptret,idealset,profset,...
    avgbrokervar,carrycost,Lottype,CONINFO,PREFMAP,TSI,IMPACT,DAMAGE,LANDINFO,...
    lotlocate)

cd X:\model_results\CHALMS_coast_032315

save(savefname,...
    'consumerstats','vacstats','BUILDTIME','VACLAND','RENT','RETURN',...
            'LOTTYPE','BASELAYER','Rpop','Rvacrate','Rvaclots',...
            'numlt','Rleftoverpop','avgrentdynms','rentdynstats',...
            'farmsoldinfo','avgrent','avgfarmsize','stdfarmsize','DELTA',...
            'survivalrate','LOCWGHT','REGWGHT','PCTSEARCH','zeta','HIBETA',...
            'MIDBETA','LOWBETA','POPGROW','ccost','newhouseset','bidtot',...
            'meandisp','maxdevdist','setupmap','zonedensity',...
            'vacrate','Farminfo','oldincome','Realreturn',...
            'Realavgret','Exptrentdiff','Avgexptdiff','htincome',...
            'numtotbids','totfarmrecord','htperyear','Newbidlevel',...
            'totbrokerrecord','Farmdist2dev','Bprojmap','Bmodelmap',...
            'WTAlandmap','WTPlandmap','Landprojmap','Landmodmap',...
            'bidshare','buildshare','landdemand','EXPTHOUSE','lotchoice',...
            'Exptprofit','Exptret','Realexptret','Realavgexptret','idealset',...
            'profset','avgbrokervar','carrycost','Lottype','CONINFO','PREFMAP',...
            'TSI','IMPACT','DAMAGE','LANDINFO','lotlocate')
        
cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\base-chalms-code

end
        
        