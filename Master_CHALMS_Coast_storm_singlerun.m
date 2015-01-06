%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%    MASTER CODE    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
tic

EXPTRUNS=1;
MRUNS=1;
% parpool(max(EXPTRUNS,12))
stream=RandStream.create('mrg32k3a','Seed',13);
RandStream.setGlobalStream(stream);

% stream.Substream=86;
% repeatstate1=stream.State;
% stream.Substream=35;
% repeatstate2=stream.State;
% stream.Substream=4;
% repeatstate3=stream.State;

repeatstate=cell(3,EXPTRUNS);
rsind=[86 35 4];
for ic=1:EXPTRUNS
    for ir=1:3
        stream.Substream=rsind(ir);
        repeatstate{ir,ic}=stream.State;
    end
end

parmfit=zeros(MRUNS,7);

% poolobj=parpool(12);
% % addAttachedFiles(poolobj,{'RefLandscape_Coast_batch.m','CHALMS_Coast_batch.m',...
% %     'load_expmntlparms.m','loadempdata.m','HouseMarketInitial_coast_batch.m',...
% %     'HouseMarketDynamic_coast_batch.m','LandMarket_coast_base.m','parsave.m',...
% %     'FarmerModule_Coast_base.m','BrokersModule_Coast_0v3.m','distmat.m',...
% %     'dist2coast.mat','load_farmmap.m','load_DIST2CBD_east.m','load_distmat.m'});
% addAttachedFiles(poolobj,{'load_expmntlparms_storm.m','loadempdata.m',...
%     'parsave_storm.m','distmat.m','load_farmmap.m','load_DIST2CBD_east.m',...
%     'load_distmat.m'});


%%
for erun=1:EXPTRUNS
    
    %     rndstr.SubStream=erun;
   for mrun=1:MRUNS
        %%
        rndstr=RandStream.getGlobalStream;
        cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\base-chalms-code
        rndstr.Substream=mrun;
        
%         disp([erun mrun])
        
        % load experimental parameters file
        [am0,am_slope,ampref_max,ampref_min,maxPflood,highrisk,stormfreq,maxdam,...
            Cmit,miteff,AVGFARMRETURN,STDFARMRETURN,coastvalue,midvalue,...
            inlandvalue,milecost,milestraveled]=load_expmntlparms_storm(EXPTRUNS);
        %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        %     % M-files needed to run
        %     1. EmpDataInput_coast_base.m
        %     2. LandMarket_coast_base.m
        %     3. HouseMarketInitial_coast_base.m
        %     4. HousemarketDynamic_coast_base.m
        %     5. CHALMS_Coast_base.m
        %     6. FarmerModule_Coast_base.m
        %     7. BrokerModule_coast_base.m
        %     8. Reflandscape_Coast_base.m
        %     9. DIST2CBD.mat
        %     10. master_dist.mat
        %     11. MasterRun_CHALMS_Coast_base.m
        %     12. FARMMAP.mat
        %     13. distmat.m
        %     14. distcalc.m
        %     15. GetResults_CHALMS_Coast_base.m
        %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        %%% Housing Layer %%%
        HTYPE=2;
        LTYPE=4;
        HT=HTYPE*LTYPE;
        
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %@@@@@@@@@@@@@@@@@@@@    INITIAL CONDITIONS    @@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        discount=0.05;
        %     ccost_base=([208273.7676
        %     270773.7676
        %     333273.7676
        %     219872.3592
        %     282372.3592
        %     344872.3592
        %     268120.5986
        %     360620.5986
        %     423120.5986])*discount;
        ccost_base=([201624.1197
            326624.1197
            208273.7676
            333273.7676
            217748.2394
            342748.2394
            232872.3592
            357872.3592])*discount;
        ccost=ccost_base;
        
        resnum=[3963; 10230; 7892; 3947; 2096; 1255];
        rtstart=[1 2; 3 7; 8 17; 18 27; 28 37; 38 47; 48 50];
        inflate=207.342/113.6;
        presqftcost=100*ones(1,HT);
        streetcost=[7000; 7000; 7000; 15000; 15000; 15000; 25000; 25000; 25000]; %1987 dollars
        sewercost=[8000; 8000; 8000; 18000; 18000; 18000; 5000; 5000; 5000];
        incomenum=round(116783.*[0.50 .35 0.15]);
        inspan=[40000 69999; 80000 119999; 120000 200000];
        returnmeans=[201.17 150.17 51.82 570.94];
        sizemans=[52 65 70 31];
        
        [restimedata,avgrestime,stdrestime,strtcst,swrcst,infracost,...
            incomedata,parmhat,avgincome,stdincome]=loadempdata(HT,rtstart,resnum,...
            presqftcost,streetcost,sewercost,inflate,incomenum,inspan);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%    PARAMETERS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        DELTA=0.0527;
        survivalrate=0.9499;
        LOCWGHT=0.2955;
        REGWGHT=1-LOCWGHT;
        PCTSEARCH=0.5887;
        zeta=0.5;   %dampening on epsilon-- must be >0, controls the rate of price ...
        %increase/preception of market power in land market
        learnDELTA=0.5;
        TMAX=30;
        TSTART=10;
        
        %%% Landscape Layer %%%
        NWIDTH=80;
        NLENGTH=80;
        NCELLS=NLENGTH*NWIDTH;
        mile2acre=640;                          %1 sq mile = 640 acres
        % 1600 acres= 2.5mi^2
        % acre2cell=4096;                         %megacell is 64x64 acres
        cellside=1;   %linear distance
        cell2mile=0.0395;   %cell side equals 0.0395 mi
        cell2ft=cell2mile*5280;
        celldiag=(cellside*sqrt(2));          %miles diag to neighboring cell's center
        acre2sqft=43560;
        avgdist2nei=mean([cellside celldiag]);
        margtc=milecost(erun)*milestraveled(erun)*cell2mile;        %Assumed: 250 travel days a year, roundtrip
        %%% Zones Layer %%%
        NZONES=25;
        
        %%% Broker Layer %%%
        Nbrokers=(NLENGTH/5)^2;
        
        %%% Agricultural Layer %%%
        Nfarmers=64;
        
        %%% Price Projection Models %%%
        NUMMODEL=20;
        FARMNUMCLASS=6;
        POPNUMCLASS=5;
        BROKERNUMCLASS=6;
        MAXMEANMODEL=10;
        MAXCYCLEMODEL=10;
        MAXPROJECT=10;
        % DELTA=1/50;
        
        %%% Farmers %%%
        FARMPROD=450;
        FARMCOST=225;
        PRODSTD=20;
        COSTSTD=50;
        % distance models for farmers
        NUMMODELDIST=100;
        maxcoeff=200;
        mincoeff=-200;
        % survivalrate=0.1;
        selectivity=NUMMODELDIST*survivalrate;
        nextgen=10;
        recombo=40;
        
        %%% Developer %%%
        Ndevelopers=1;
        inmaxrent=10000;
        CLOSECELL=30;
        
        %%% Consumers %%%
        Nconstart=270;
        
        % RESSPAN=7;
        maxwage=max(incomedata);
        minwage=min(incomedata);
        himaxwage=maxwage;
        himinwage=100000;
        midmaxwage=99999;
        midminwage=60000;
        lowmaxwage=59999;
        lowminwage=minwage;
        WAGECLASS=3;
        HIBETA=[0.20 0.24];     %approximate max and mins for housepref come from calvert county census data
        MIDBETA=[0.25 0.29];
        LOWBETA=[0.3 0.5];
        incomeg=0.05;
        incomesigma=1;
        incomegrow=0.005;
        
        %%% Population %%%
        popg=0.05;
        popsigma=1;
        pop2dem=2.91;   %persons per household, Calvert County quickfacts
        POPGROW=0.10;   %ditto (19% actually)
        popurb=0.7;
        popag=0.3;
        
        tax=0.05;
        thetavac=0.2;
        thetadev=0.2;
        
        %Flags
        epflag=0;
        oldbuildflag=0;
        newbuildflag=0;
        termrunflag=0;
        
        % Parameters for developer's risk aversion models
        alpha_gain=3;   % skewedness factor (Ligmann-Zielinska, 2009)
        alpha_loss=2.5;
        w_gain=0.5;
        w_loss=0.5;
        proftarget=5000*discount;
        
        %     stormthresh=15;
        %     Cdam=0.5;
        %     Cmit=3000;
        %     miteff=1;
        
        testtime=[10 15 20 25 30];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%    VARIABLES    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%% Landscape Layer %%%
        BASELAYER=zeros(NCELLS,1);  %developed or undeveloped
        travelcost=cell(NCELLS,1);
        IMPACT=cell(NCELLS,TMAX);
        DAMAGE=num2cell(zeros(NCELLS,1));
        MITIGATE=cell(NCELLS,1);
        Cdam=cell([],TMAX);
        
        COAST=zeros(NLENGTH,NWIDTH);
        SCAPE=zeros(NLENGTH,NWIDTH);
        
        %%% Zones Layer - for use with zoning settings %%%
        ZONES=cell(NZONES,3);   %[cellid,min zone,max zone]
        Zones=zeros(NLENGTH,NWIDTH);
        zoning=zeros(NZONES,2);     %[(min lotsize) (max lotsize)]
        ZONEMAP=zeros(NLENGTH,NWIDTH);
        
        nzoneslong=sqrt(NZONES);
        nzoneswide=sqrt(NZONES);
        extralengthz=rem(NLENGTH,nzoneslong);
        extrawidthz=rem(NWIDTH,nzoneswide);
        izonelength=(NLENGTH-extralengthz)/nzoneslong;
        izonewidth=(NWIDTH-extrawidthz)/nzoneswide;
        
        %%% Housing Layer %%%
        % %Lottype=[id,location_index,lotsize,housesize,ltype,ccost,amlevel,travelcost,buildtime,brokerid]
        % %lotchoice=[id,location_index,ltype,occ/vac,consumer_id,residence_time,sell_price,mitchoice]
        Lottype=cell([],10);
        lotchoice=cell([],8);
        HOUSESIZE=zeros(NLENGTH,NWIDTH);
        AMLEVEL=zeros(NLENGTH,NWIDTH);
        LOTS=zeros(NLENGTH,NWIDTH);
        
        %%% Broker Layer %%%
        BROKER=cell(Nbrokers,4);
        HBROKER=zeros(NLENGTH,NWIDTH);
        UUcell=zeros(NLENGTH,NWIDTH);
        
        %%% Agricultural Layer %%%
        LANDINFO=cell(3,TMAX);    %[farmerid,land value,PLAND]
        
        %%% Farmers %%%
        Farminfo=cell(Nfarmers,2);   %[acres prod costs value_acre]
        farmretinfo=zeros(Nfarmers,1);  %baseline agriculturual return
        wtaland=zeros(Nfarmers,TMAX);   %farmer willingness to accept price
        Paskland=zeros(Nfarmers,TMAX);  %farmer asking price for land
        transprice=zeros(Nfarmers,TMAX);    %transaction price for land
        Plandproj=zeros(Nfarmers,TMAX);     %farmer projected price for land
        maxrent=zeros(Nfarmers,TMAX);
        maxreturn=zeros(Nfarmers,TMAX);
        landsold=zeros(NZONES,TMAX);
        
        sellrecord=zeros(Nfarmers,1);   %time of land sale
        pdrrecord=zeros(Nfarmers,1);
        buyrecord=zeros(Nfarmers,1);    %sale price for land
        
        fitness=zeros(Nfarmers,NUMMODELDIST,TMAX);
        
        % Price Projection models
        landproj = zeros(Nfarmers,NUMMODEL);
%         landerror = zeros(Nfarmers,NUMMODEL);
        landbestSAVE=zeros(Nfarmers,TMAX);
        ilandbestSAVE=zeros(Nfarmers,TMAX);
        landprojSAVE=zeros(Nfarmers,TMAX);
        landprojbestSAVE=zeros(Nfarmers,NUMMODEL);
        landmodelbestSAVE=zeros(Nfarmers,TMAX);
        landmodelSAVE=zeros(Nfarmers,TMAX);
        learnlandproj = zeros(Nfarmers,NUMMODEL);
        learnlanderror = zeros(Nfarmers,NUMMODEL);
        learnlandbestSAVE=zeros(Nfarmers,TMAX);
        ilearnlandbestSAVE=zeros(Nfarmers,TMAX);
        learnlandprojSAVE=zeros(Nfarmers,NUMMODEL);
        learnlandmodelSAVE=zeros(Nfarmers,TMAX);
        learnlandmodelbestSAVE=zeros(Nfarmers,TMAX);
        learnlandprojbestSAVE=zeros(Nfarmers,TMAX);
        learnwtaland=zeros(Nfarmers,TMAX);
        
        
        landprojdiff=zeros(Nfarmers,TMAX);
        landpctdiff=zeros(Nfarmers,TMAX);
        landpricevar=zeros(Nfarmers,TMAX);
        projdiff=zeros(Nfarmers,TMAX);
        pctdiff=zeros(Nfarmers,TMAX);
        pricevar=zeros(Nfarmers,TMAX);
        
        %%% Developer %%%
        Profit=zeros([],1);
        Rmin=zeros([],1);      %minimum rent
        Paskhouse=zeros([],1);  %house asking price
        Ccost=zeros([],1);
        Pland=zeros([],1);
        Plot=zeros([],1);
        PROFIT=zeros(NCELLS,TMAX);
        RETURN=cell(NCELLS,TMAX);
        RENT=zeros(NCELLS,TMAX);
        RENTPROJ=cell(NCELLS,TMAX);
        subRENTPROJ=zeros(HT,1);
        subRETURN=zeros(HT,1);
        MAXRETURN=zeros(NCELLS,1);  %maximum return on development per cell
        RETIND=zeros(NCELLS,1); %index of house type that yields max return
        EUIND=zeros(HT,NCELLS);
        MAXEUIND=zeros(NCELLS,1);   %index of house type ordered by expected utility
        Uret=zeros(NCELLS,HT);
        PCHANGE=zeros(NCELLS,HT);
        RETURNPROJ=zeros(NCELLS,1);
        LOTPRICE=zeros(NCELLS,1);
        LOTVALUE=zeros(NCELLS,1);
        INITIALPASK=zeros(NCELLS,1);
        house2cell=zeros(NCELLS,1);
        RENTGRAD=zeros(NCELLS,1);
        RENTPROJLAND=zeros(NCELLS,1);
        maxRENTPROJ=zeros(NCELLS,1);
        CCOST=zeros(NCELLS,1);
        ZZ=zeros(NCELLS,1);
        htset=(1:HT)';
        pctbuildnew=zeros(HT,TMAX);
        pctbuildold=zeros(HT,TMAX);
        newhouses=zeros(HT,TMAX);
        newbuild=zeros(HT,TMAX);
        newacres=zeros(HT,TMAX);
        oldhouses=zeros(HT,TMAX);
        bidtot=zeros(HT,TMAX);
        bidlevel=zeros(HT,TMAX);
        newhouseset=zeros(HT,TMAX);
        numnewhouses=zeros(1,TMAX);
        numoldhouses=zeros(1,TMAX);
        numnewacres=zeros(1,TMAX);
        Pdevbid=zeros(Nfarmers,TMAX);
        wtpland=zeros(Nfarmers,TMAX);
        wtahouse=zeros([],1);
        dist2vac=zeros(NLENGTH,NWIDTH);
        testdist2dev=zeros(NLENGTH,NWIDTH);
        subprofit=zeros(NLENGTH,NWIDTH);
        subrentproj=zeros(NLENGTH,NWIDTH,HT);
        epsilon=zeros(1,TMAX);
        landdemand=zeros(Ndevelopers,TMAX);
        regionaldist=zeros(HT,TMAX);
        regionalrent=zeros(HT,TMAX);
        simlotrange=zeros(HT,2);
        simlots_income=zeros(HT,1);
        simlots_util=zeros(HT,1);
        simlots_alpha=zeros(HT,1);
        simlots_beta=zeros(HT,1);
        simlots_gamma=zeros(HT,1);
        simlots_ampref=zeros(HT,1);
        vac_ccost=zeros(HT,TMAX);
        vac_rent=zeros(HT,TMAX);
        vac_land=zeros(TMAX,1);
        profits=zeros(HT,TMAX);
        budget_lt=zeros(HT,TMAX);
        carrycost=zeros(TMAX,1);
        BUDGET=zeros(TMAX,1);
        LANDBUDGET=zeros(1,TMAX);
        WTPMAP=zeros(NCELLS,TMAX);
        highRETURN=zeros(HT,1); %projected max return for usein developer's risk aversion models
        lowRETURN=zeros(HT,1);  %projected min return for usein developer's risk aversion models
        potgain=zeros(HT,1);    %potential gain based on highRETURN
        potloss=zeros(HT,1);
        EU_dev=cell(NCELLS,TMAX);
        subEU_dev=zeros(HT,1);
        MAXRET=zeros(NCELLS,TMAX);
        maxcount=zeros(NCELLS,1);
        
        bestPOPSAVE=zeros(Ndevelopers,TMAX);
        ibestPOPSAVE=zeros(Ndevelopers,TMAX);
        nprojSAVE=zeros(Ndevelopers,TMAX);
        
        %%% Consumers %
        avghousemp=zeros(1,TMAX);
        % %CONINFO=[income,searchtime,consumer_good,housesize,lotsize,proximity,subrisk,occ/vac,utility]
        CONINFO=cell(Nconstart,9);
        damcoef=cell(Nconstart,TMAX);   %coefficient that sets expected damages from storms for every lot
        searchtimemin=2;
        searchtimemax=6;
        
        %%% Population %%%
        POP=zeros(1,TMAX);
        
        %%% Brokers %%%
        avgconinfo=cell(Nbrokers,6);  %[income,utility,alpha,beta,gamma,ampref]
        brokerlotinfo=cell(Nbrokers,7);  %[avg_price lotsize housesize #of_bidders %above_Pask #lots(lt) amlevel]
%         AVGUTIL=cell([],1);
        brkravgstats=zeros(Nbrokers,5); %[avg_income avg_alpha avg_beta avg_gamma avg_ampref]
        houseinfo=zeros(HT,7,Nbrokers,TMAX);      %[avg_price lotsize housesize #of_bidders %above_Pask #lots(lt) amlevel]
        EXPTHOUSE=zeros(NCELLS,TMAX);
        MINBIDLEVEL=cell([],1);
        bproj=zeros(HT,NUMMODEL);
        brokerproj = zeros(HT,NUMMODEL,Nbrokers);
        brokererror = zeros(HT,NUMMODEL,Nbrokers,TMAX);
        brokerabserror = zeros(HT,NUMMODEL,Nbrokers,TMAX);
        brokerbestSAVE=zeros(Nbrokers,HT,TMAX);
        ibrokerbestSAVE=zeros(Nbrokers,HT,TMAX);
        brokerprojSAVE=zeros(Nbrokers,HT,TMAX);
        brokermodelSAVE=zeros(Nbrokers,HT,TMAX);
        numlt=zeros(HT,TMAX);
        newopenlots=zeros(HT,TMAX);
        avgbrokervar=zeros(HT,TMAX);
        probloss=zeros(HT,TMAX);
        probeven=zeros(HT,TMAX);
        probover=zeros(HT,TMAX);
        probunder=zeros(HT,TMAX);
        overvalue=zeros(HT,TMAX);
        undervalue=zeros(HT,TMAX);
        maxvalue=zeros(HT,TMAX);
        minvalue=zeros(HT,TMAX);
        avgover=zeros(HT,TMAX);
        avgunder=zeros(HT,TMAX);
        diffbrokererror=zeros(HT,NUMMODEL,Nbrokers,TMAX);
        brokerbestdiffSAVE=zeros(Nbrokers,HT,TMAX);
        brokerbestabsSAVE=zeros(Nbrokers,HT,TMAX);
        phat=zeros(HT,2);
        EUrankret=zeros(NCELLS,HT);
        EUlandret=zeros(NCELLS,HT);
        bcheck=zeros(HT,Nbrokers);
        
        % Price Projection models
        learnbproj=zeros(HT,NUMMODEL);
        learnbrokerproj = zeros(HT,NUMMODEL,Nbrokers);
        learnbrokererror = zeros(HT,NUMMODEL,Nbrokers,TMAX);
        learnabserror = zeros(HT,NUMMODEL,Nbrokers,TMAX);
        learnbrokerbestSAVE=zeros(Nbrokers,HT,TMAX);
        ilearnbrokerbestSAVE=zeros(Nbrokers,HT,TMAX);
        learnbrokerprojSAVE=zeros(Nbrokers,HT,TMAX);
        learnbrokermodelSAVE=zeros(Nbrokers,HT,TMAX);
        difflearnerror=zeros(HT,NUMMODEL,Nbrokers,TMAX);
        learnbestdiffSAVE=zeros(Nbrokers,HT,TMAX);
        learnbestabsSAVE=zeros(Nbrokers,HT,TMAX);
        brkrprojdiff=zeros(HT,Nbrokers,TMAX);
        brkrpctdiff=zeros(HT,Nbrokers,TMAX);
        brkrpricevar=zeros(HT,Nbrokers,TMAX);
        
        % Input sample rents
        testrents=10000*[0.6721    0.6999    0.8438    0.9093    0.9526
            0.8860    0.9349    1.1483    1.2432    1.2900
            1.2645    1.3129    1.5328    1.6199    1.6786
            1.9035    1.9304    1.9204    1.9743    2.0447
            2.0947    2.1348    2.1311    2.1762    2.2347
            2.4736    2.5188    2.5519    2.6037    2.6382
            02.35     2.37      2.42      2.46      2.5
            2.55      2.6       2.63      2.67      2.7
            2.83      2.89      2.94      2.96      3.0];
        pcoeffs=zeros(HT,2);
        
        %%%% RESULTS %%%%
        devcells=zeros(2,TMAX);         %[#dev_cells; pct_dev_land]
        vacstats=zeros(HT,TMAX);
        avgrent=zeros(HT,TMAX);
        vacrate=zeros(1,TMAX);
        nohouserate=zeros(1,TMAX);
        oldincome=zeros(1,TMAX);
        consumerstats=zeros(4,TMAX);     %[avg_income avg_util housemp income_out]
        vacantlots=zeros(1,TMAX);
        leftoverpop=zeros(1,TMAX);
        agrland=zeros(1,TMAX);
        vacland=cell(1,TMAX);
        BUILDTIME=zeros(NCELLS,1);
        BIDLEVELMAP=zeros(NCELLS,TMAX);
        VACLAND=zeros(NCELLS,TMAX);
        AVGRENT=zeros(NCELLS,TMAX);
        RETURNMAP=zeros(NCELLS,TMAX);
        LOTTYPE=zeros(NCELLS,TMAX);
        INCOME=zeros(NCELLS,TMAX);
        BASEMAP=zeros(NCELLS,1);
        PREFMAP=zeros(NCELLS,TMAX);
        SUBRISKMAP=zeros(NCELLS,TMAX);
        OutIncome=zeros([],1);
        
        Rpop=zeros(1,TMAX);
        Rvacrate=zeros(1,TMAX);
        Rvaclots=zeros(1,TMAX);
        Rleftoverpop=zeros(1,TMAX);
        idealset=zeros(HT,TMAX);
        profset=zeros(HT,TMAX);
        Exptrentdiff=zeros(NCELLS,TMAX);
        Avgexptdiff=zeros(HT,TMAX);
        Realreturn=zeros(NCELLS,TMAX);
        Realexptret=zeros(NCELLS,TMAX);
        Realavgret=zeros(HT,TMAX);
        Realavgexptret=zeros(HT,TMAX);
        Newbidlevel=zeros(NCELLS,TMAX);
        Avgnewbid=zeros(HT,TMAX);
        
        numtotbids=zeros(HT,TMAX);
        htincome=zeros(HT,TMAX);
        htperyear=zeros(HT,TMAX);
        Exptprofit=zeros(HT,TMAX);
        Exptret=zeros(HT,TMAX);
        
        totbrokerrecord=zeros(3,HT,TMAX,Nbrokers);
        
%         totfarmrecord=zeros(4,TMAX,Nfarmers);   %[wtaland; wtpland; landprojSAVE; landmodelSAVE]
%         Allfarmdist=zeros(Nfarmers,1);
        Farmdist2dev=zeros(Nfarmers,TMAX);
        Dynltmap=zeros(NCELLS,TMAX);
        Dynrentmap=zeros(NCELLS,TMAX);
        Dynmaxretmap=zeros(NCELLS,TMAX);
        Dynretltmap=zeros(NCELLS,TMAX);
        Dynmaxprofmap=zeros(NCELLS,TMAX);
        Dynprofltmap=zeros(NCELLS,TMAX);
        Dynplandmap=zeros(NCELLS,TMAX);
        Dyneultmap=zeros(NCELLS,TMAX);
        ZONED=zeros(NCELLS,1);
        
        subdynlt=zeros(NCELLS,1);
        subdynrent=zeros(NCELLS,1);
        Bmodelmap=zeros(NCELLS,TMAX);
        Bprojmap=zeros(NCELLS,TMAX);
        
        savedState=rndstr.State;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%    Load Reference Landscape    %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% RefLandscape_Coast_batch
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %@                             LANDSCAPE LAYERS                           @
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        
        % rand('state',86);
        % s=rand('state');
        % randn('state',86)
        % sn=randn('state');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%    Physical Landscape    %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Map layer input
        %%% Coast is the eastern edge of the region
        % COAST(:,1:5)=1;
        icoast=find(COAST==1);
        SCAPE(COAST~=1)=1;
        iscape=(SCAPE==1);
        iscapelist=find(iscape==1);
        
        coastdist=NWIDTH+1-cumsum(SCAPE,2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Distance matrices   %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        distfname='C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\data_files\DIST2CBD_east.mat';
        Sdist2cbd=load_DIST2CBD_east(distfname);
        %         load DIST2CBD_east
        
        % icenterrow=1;
        % % icentercol=round(NWIDTH/2);
        % icentercol=NWIDTH;
        %
        % for col=1:NWIDTH
        %     dist2hznnei(1:NLENGTH,col)=abs(col-icentercol).*ones(NLENGTH,1);
        % end
        %
        % for row=1:NWIDTH
        %     dist2vrtnei(row,1:NWIDTH)=abs(row-icenterrow).*ones(1,NLENGTH);
        % end
        %
        % for col=1:NWIDTH
        %     for row=1:NLENGTH
        %         DISTANCE(row,col)=sqrt(dist2hznnei(row,col)^2+dist2vrtnei(row,col)^2);
        %     end
        % end
        % am_int=500000;
        travelcost(iscapelist)=num2cell(margtc*Sdist2cbd.dist2cbd(iscapelist));
        travelcost(icoast)=num2cell(10000*ones(length(icoast),1));
        % coastprox=num2cell(reshape(10*(max(max(coastdist))+1-coastdist),NCELLS,1));
        coastprox=num2cell(reshape(am0(erun)*1./exp(am_slope(erun)*coastdist),NCELLS,1));
        % coastprox=num2cell(reshape(0.1*(max(max(coastdist))+1-coastdist),NCELLS,1));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%     Impact Surface    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Probability of a hurricane of a given category making landfall in a
        %%% given year, from Costanza et al. (2008), Ambio, 37(4)
       
        % Mid-Atlantic = Delaware, Maryland, New Jersey, North Carolina, Pennsylvania, Virginia
        ma_storm=[0.013 0 0 0 0; 0.0065 0.0065 0 0 0; 0.013 0 0 0 0; 0.1364 0.0844 ...
            0.0714 0.0065 0; 0.0065 0 0 0 0; 0.0584 0.013 0.0065 0 0];
        nc_storm=[0.1364 0.0844 0.0714 0.0065 0];
        fl_storm=[0.2792 0.2078 0.1753 0.0390 0.0130];
        tx_storm=[0.1494 0.1104 0.0779 0.0455 0];
        
        stormprob=cell(4,1);    %[mid-Atlantic average, NC, FL, TX]
        stormprob(1)=mat2cell(mean(ma_storm,1),1,5);
        stormprob(2)=mat2cell(nc_storm,1,5);
        stormprob(3)=mat2cell(fl_storm,1,5);
        stormprob(4)=mat2cell(tx_storm,1,5);
        Psevere=stormprob{4};
%         Psevere=stormprob{stormfreq(erun)};
%         Psevere=stormfreq(erun)*[0.013 0 0 0 0; 0.0065 0.0065 0 0 0; 0.013 0 0 0 0; 0.1364 0.0844 ...
%             0.0714 0.0065 0; 0.0065 0 0 0 0; 0.0584 0.013 0.0065 0 0];
        
        % maxPflood=0.7;
        % highrisk=30;
        risksurf=1./(1+exp((coastdist-highrisk(erun))/10));
        Pflood=num2cell(reshape(maxPflood(erun)*(risksurf+(1-max(max(risksurf)))),NCELLS,1));
        % Pflood=mat2cell([(1:NCELLS)' reshape(maxPflood*(risksurf+(1-...
        %     max(max(risksurf)))),NCELLS,1)],ones(NCELLS,1),2);
        TSI=num2cell(ones(NCELLS,1));
        % siterisk=0.1*randn(size(SCAPE));
        % siterisk=zeros(size(SCAPE));
        
%         Pimpact=cat(1,Pflood{:})*mean(Psevere,1);

        % Percent damage with storm categories 1-5
        housedam=maxdam(erun)*10.23749-0.23462*(coastdist*cell2ft/1000)+...
            0.001649*(coastdist*cell2ft/1000).^2;
%         housedam=maxdam(erun)*[0.05 0.1 0.25 0.5 1];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ZONES Layer    %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % nzlong=sqrt(NCELLS/NZONES);
        % zoneblock=zeros(nzlong);
        % for iz=1:NZONES
        %     startpt=floor((iz-1)/sqrt(NZONES))*(NLENGTH*nzlong)+nzlong*(iz-1)+NLENGTH*(0:1:nzlong-1)+1;
        %     endpt=floor((iz-1)/sqrt(NZONES))*(NLENGTH*nzlong)+nzlong*(iz-1)+NLENGTH*(0:1:nzlong-1)+nzlong;
        %     for ii=1:nzlong
        %         zoneblock(:,ii)=startpt(ii):endpt(ii);
        %     end
        %     ZONES{iz}=reshape(zoneblock,nzlong^2,1);
        % end
        
        zonemarklong=1;
        zonemarkwide=1;
        for ii=1:nzoneswide
            for jj=1:nzoneslong
                Zones(zonemarklong:zonemarklong+izonelength-1,...
                    zonemarkwide:zonemarkwide+izonewidth-1)=...
                    ii*nzoneslong-(nzoneslong-jj);
                
                if jj==nzoneslong && extralengthz > 0
                    Zones(izonelength*jj+1:izonelength*jj+extralengthz,...
                        zonemarkwide:zonemarkwide+izonewidth-1)=...
                        ii*nzoneslong-(nzoneslong-jj);
                end
                if ii==nzoneswide && extrawidthz > 0
                    Zones(zonemarklong:zonemarklong+izonelength-1,...
                        izonewidth*ii+1:izonewidth*ii+extrawidthz)=...
                        ii*nzoneslong-(nzoneslong-jj);
                end
                zonemarklong=mod(izonelength*jj+1,izonelength*nzoneslong);
            end
            if jj==nzoneslong && ii==nzoneswide && extralengthz > 0
                Zones(izonelength*jj+1:izonelength*jj+extralengthz,...
                    izonewidth*ii+1:izonewidth*ii+extrawidthz)=...
                    ii*nzoneslong-(nzoneslong-jj);
            end
            zonemarkwide=mod(izonewidth*ii+1,izonewidth*nzoneswide);
        end
        for iz=1:NZONES
            ZONES{iz,1}=find(Zones==iz);
        end
        % for n=1:NZONES
        %     str= sprintf('izone%d = find(ZONES == %d);',n,n);
        %     eval(str);
        % end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ZONING   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        house2acre=[4 2 1 0.5]';
        zs=min(max(1./house2acre,0),2);
        
        for zt=1:NZONES
            ZONES(zt,2:3)=num2cell([min(zs) max(zs)]); %no zoning
            %     zoning(zt,:)=[1 1];
            %     if isempty(find(zt==[1:2 6:8 11:13 16:18 21:22],1))==0
            %         zoning(zt,:)=[min(zs) max(zs)]; %no zoning
            %         ZONEMAP(ZONES==zt)=zoning(zt,1);
            %         ZONED(ZONES==zt)=0;
            %     elseif isempty(find(zt==[3:5 9:10 14:15 19:20 23:25],1))==0
            % %         zoning(zt,:)=[5 max(zs)];
            %         zoning(zt,:)=[min(zs) max(zs)];
            %         ZONEMAP(ZONES==zt)=zoning(zt,1);
            %         ZONED(ZONES==zt)=1;
            %     end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%    Broker Layer    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        nbrokerlong=(NLENGTH/5);
        nbrokerwide=(NWIDTH/5);
        extralengthb=rem(NLENGTH,nbrokerlong);
        extrawidthb=rem(NWIDTH,nbrokerwide);
        ibrokerlength=(NLENGTH-extralengthb)/nbrokerlong;
        ibrokerwidth=(NWIDTH-extrawidthb)/nbrokerwide;
        brokermarklong=1;
        brokermarkwide=1;
        for ii=1:nbrokerwide
            for jj=1:nbrokerlong
                HBROKER(brokermarklong:brokermarklong+ibrokerlength-1,...
                    brokermarkwide:brokermarkwide+ibrokerwidth-1)=...
                    ii*nbrokerlong-(nbrokerlong-jj);
                
                if jj==nbrokerlong && extralengthb > 0
                    HBROKER(ibrokerlength*jj+1:ibrokerlength*jj+extralengthb,...
                        brokermarkwide:brokermarkwide+ibrokerwidth-1)=...
                        ii*nbrokerlong-(nbrokerlong-jj);
                end
                if ii==nbrokerwide && extrawidthb > 0
                    HBROKER(brokermarklong:brokermarklong+ibrokerlength-1,...
                        ibrokerwidth*ii+1:ibrokerwidth*ii+extrawidthb)=...
                        ii*nbrokerlong-(nbrokerlong-jj);
                end
                brokermarklong=mod(ibrokerlength*jj+1,ibrokerlength*nbrokerlong);
            end
            if jj==nbrokerlong && ii==nbrokerwide && extralengthb > 0
                HBROKER(ibrokerlength*jj+1:ibrokerlength*jj+extralengthb,...
                    ibrokerwidth*ii+1:ibrokerwidth*ii+extrawidthb)=...
                    ii*nbrokerlong-(nbrokerlong-jj);
            end
            brokermarkwide=mod(ibrokerwidth*ii+1,ibrokerwidth*nbrokerwide);
        end
        
        minibmap=reshape(unique(HBROKER),nbrokerlong,nbrokerwide);
        for ib=1:Nbrokers
            BROKER{ib,1}=find(HBROKER==ib);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%% Landscape Template %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % stream.Substream=86;
        rndstr.State=repeatstate{1,erun};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%    Housing Layer    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Nlots=zeros(1,TMAX);
        inurbarea=288;
        globeind=[(1:NLENGTH*NWIDTH)' reshape(Sdist2cbd.dist2cbd,NLENGTH*NWIDTH,1)];
        globeind=globeind(ismember(globeind(:,1),iscapelist),:);
        iurbprox=sortrows(globeind,2);
        iurblist=iurbprox(1:inurbarea,1);
        iagrlist=globeind(~ismember(globeind(:,1),iurblist),1);
        BASELAYER(iurblist)=ones(length(iurblist),1);
        BASELAYER(iagrlist)=zeros(length(iagrlist),1);
        BASELAYER(icoast)=zeros(length(icoast),1);
        indevedge=find(BASELAYER==1,1,'last');
        
        
        %<><><><><><><> Non-random, 1 acre distribution <><><><><><><><><><><><><>
        % Nlots(1:TSTART+1)=length(iurblist);
        % for i=1:length(iurblist)
        %     LOTS(iurblist(i))=i;
        % end
        % t=1;
        % ilots=LOTS(iurb);
        % LOTS(iurb)=1:length(iurblist);
        
        %<><><><><><><><><><><> Randomized Lot Allocation <><><><><><><><><><><><>
        Nlots(1)=170;
        % Nlots(1)=225;
        openlots=length(iurblist)-Nlots(1);
        Nlotstry=round(openlots*1.2);
        bridgelayer=zeros(size(BASELAYER));
        bridge=iurblist;
        openstartpos=iurblist(min(ceil(length(iurblist)*rand(Nlotstry,1)),length(iurblist)));
        while length(unique(openstartpos)) <= openlots
            for k=1:length(openstartpos)
                check=(openstartpos(k)==openstartpos);
                if length(find(check==1)) > 1
                    openstartpos(check)=iurblist(min(ceil(length(iurblist)*rand(length(find(check==1)),1)),length(iurblist)));
                end
            end
        end
        iopenunique=unique(openstartpos);
        iopenunique=iopenunique(1:openlots);
        bridgelayer(iopenunique)=1;
        ibridge=(bridgelayer(iurblist)==0);
        lotstartpos=iurblist(ibridge);
        LOTS(lotstartpos)=(1:Nlots(1));
        LOTS(~ismember((1:NCELLS),iurblist))=Nlots(1)+1000;
        % LOTS(SCAPE==0)=Nlots(1)+1000;
        
        %%% Assign 2 acre lots
        adjcells=diff(find(LOTS==0));
        if isempty(find(adjcells==1,1))==0
            for f=1:Nlots(1)
                mylots=find(LOTS==f);
                for ml=1:length(mylots)
                    %             zoneid=ceil(find(cat(2,ZONES{:,1})==mylots(ml))/length(ZONES{1,1}));
                    %             zoneid=ZONES{mylots(ml)};
                    growthlim=2-length(mylots);
                    if growthlim == 0
                        continue
                    else
                        [frows,fcols]=ind2sub([NLENGTH,NWIDTH],mylots(ml));
                        rnei=min((fcols+1),NWIDTH)*NLENGTH-(NLENGTH-frows);
                        lnei=max((fcols-1),0)*NLENGTH-(NLENGTH-frows);
                        upnei=fcols*NLENGTH-(NLENGTH-max((frows-1),0));
                        dnnei=fcols*NLENGTH-(NLENGTH-min((frows+1),NLENGTH));
                        distcells=[rnei;lnei;upnei;dnnei];
                        edgecells=find(LOTS(distcells)==0);
                        if isempty(edgecells)==1
                            continue
                        else
                            LOTS(distcells(edgecells(1:growthlim)))=f;
                        end
                    end
                end
            end
            adjcells=diff(find(LOTS==0));
        end
        ilotfill=find(LOTS==0 & SCAPE==1);
        LOTS(ilotfill)=Nlots(1)+(1:length(ilotfill));
        %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
        
        housesize=[1500 2500]';
        house2cells=unique(LOTS((LOTS < 250)));
        
        %%%% Amenity Level
        % amlevel=[100 300 200 100 300 200 400]';
        % amlevel=[100 200 150 100 200 150 400]';
        
        for ih=1:length(house2cells)
            %     HOUSECHAR(ih)=num2cell([housesize(ceil(length(housesize)*rand(1))) ...
            %         mean(coastprox{isamecell})
            isamecell=find(LOTS == house2cells(ih));
            HOUSESIZE(isamecell)=housesize(ceil(length(housesize)*rand(1)));
            %     AMLEVEL(isamecell)=amlevel(ceil(length(housesize)*rand(1)));
            %%% Starting with all conventional subdivision types %%%
            %     HOUSESIZE(isamecell)=housesize(1);
            %     AMLEVEL(isamecell)=amlevel(1);
            %     AMLEVEL(isamecell)=max(max(coastdist))-mean(coastdist(isamecell));
            AMLEVEL(isamecell)=mean(cat(1,coastprox{isamecell}));
        end
        
        z=[reshape(repmat(zs,1,HTYPE)',HT,1) repmat(housesize,LTYPE,1)];
        
        % devpref=[6; 4.67; 5.75; 4.67; 4.58; 4.33; 4.11]./7;
        % devpref=[0.851; 0.646; 0.823; 0.620; 0.604; 0.601; 0.641];
        devpref=ones(HT,1);
        
        %%% Define and locate similar lots
        for lt=1:HT
            if z(lt,1) < 1
                simlotrange(lt,1)=1;
                simlotrange(lt,2)=find(z(:,1) < 1,1,'last');
            elseif z(lt,1) >=1 && z(lt,1) < 2
                simlotrange(lt,1)=find(z(:,1) >= 1,1,'first');
                simlotrange(lt,2)=find(z(:,1) < 2,1,'last');
            elseif z(lt,1) >= 2
                simlotrange(lt,1)=find(z(:,1) >= 2,1,'first');
                simlotrange(lt,2)=find(z(:,1) >= 2,1,'last');
            end
        end
        
        %%% Add lots smaller than 1 acre %%%
        subLotinfo=[LOTS(iurblist) iurblist HOUSESIZE(iurblist) AMLEVEL(iurblist)];
        subLotinfo=sortrows(subLotinfo,1);
        Lotinfo=zeros([],4);
        singletag=zeros(length(subLotinfo),1);
        for iii=1:length(subLotinfo)
            lotlength=length(find(subLotinfo(:,1)==subLotinfo(iii,1)));
            if lotlength==1
                singletag(iii)=1;
            end
        end
        ismall=find(singletag==1);
        
        for iss=1:length(singletag)
            if singletag(iss) == 0
                Lotinfo(length(Lotinfo(:,1))+1,:)=subLotinfo(iss,:);
            elseif singletag(iss) == 1
                randsize=ceil(length(find(zs <= 1))*rand(1));
                numcell=house2acre(randsize);
                Lotinfo(length(Lotinfo(:,1))+1:length(Lotinfo(:,1))+numcell,:)=...
                    [(subLotinfo(iss,1)+(0:numcell-1))' ...
                    ones(numcell,1).*subLotinfo(iss,2) ...
                    ones(numcell,1).*subLotinfo(iss,3) ...
                    ones(numcell,1).*subLotinfo(iss,4)];
                subLotinfo(iss+1:length(subLotinfo),1)=...
                    subLotinfo(iss+1:length(subLotinfo),1)+(numcell-1);
            end
        end
        % ########################################################################
        % Procedure with real input layer will use the centroids of lots for
        % location
        %[lotid,location index,lotsize,housesize,ltype,ccost,occ/vac,amlevel,travelcost]
        meantc=zeros([],1);
        for nl=1:length(unique(Lotinfo(:,1)))
            ilot=find(Lotinfo(:,1)==nl);
            lotloc=Lotinfo(Lotinfo(:,1) == nl,2);
            Lottype{nl,1}=nl;
            lotchoice{nl,1}=nl;
            Lottype{nl,2}=lotloc;
            lotchoice{nl,2}=lotloc(1);
            ilotid=find(Lotinfo(:,1) == nl);
            ilotsize=Lotinfo(Lotinfo(:,1)==nl,2);
            Lottype{nl,3}=length(Lotinfo(ilot,1))/...
                length(find(Lotinfo(:,2)==Lotinfo(ilot(1),2)));
            Lottype{nl,4}=HOUSESIZE(lotloc(1));
            Lottype{nl,5}=find(z(:,1)==Lottype{nl,3} & z(:,2)==Lottype{nl,4});
            lotchoice(nl,3)=Lottype(nl,5);
            lotchoice{nl,4}=0;
            lotchoice{nl,5}=0;
            Lottype{nl,6}=ccost(Lottype{nl,5});
            Lottype{nl,7}=AMLEVEL(lotloc(1));
            Lottype{nl,8}=mean(cat(1,travelcost{cat(1,Lottype{nl,2})}));
            Lottype{nl,9}=TSTART;
            
            meantc(nl)=mean(Lottype{nl,9});
            %     Lottype(nl,1,1:TSTART)=num2cell(nl*ones(1,TSTART),[1,TSTART]);
            %     Lottype(nl,2,1:TSTART)=num2cell(repmat(lotloc,1,TSTART),[1,TSTART]);
            %     ilotid=find(Lotinfo(:,1) == nl);
            %     ilotsize=length(find(ismember(Lotinfo(:,2),Lotinfo(Lotinfo(:,1)==nl,2))==1));
            %     Lottype(nl,3,1:TSTART)=num2cell(repmat(length(ilotid)/length(ilotsize),1,TSTART),[1,TSTART]);
            %     Lottype(nl,4,1:TSTART)=num2cell(repmat(HOUSESIZE(lotloc(1)),1,TSTART),[1,TSTART]);
            %     Lottype(nl,5,1:TSTART)=num2cell(repmat(find(z(:,1)==Lottype{nl,3,1} & ...
            %         z(:,2)==Lottype{nl,4,1}),1,TSTART),[1,TSTART]);
            %     Lottype(nl,6,1:TSTART)=num2cell(repmat(ccost(Lottype{nl,5,1}),1,TSTART),[1,TSTART]);
            %     Lottype(nl,7,1:TSTART)=num2cell(zeros(1,TSTART),[1,TSTART]);
            %     Lottype(nl,8,1:TSTART)=num2cell(AMLEVEL(lotloc(1)),[1,TSTART]);
            %     Lottype(nl,9,1:TSTART)=num2cell(travelcost(Lottype{nl,2,1}),[1,TSTART]);
            %     meantc(nl)=mean(Lottype{nl,9,1})
        end
        Nlots(1:TSTART+1)=length(Lottype(:,1));
        BIDLEVEL=num2cell(ones(Nlots(TSTART),1));
        AVGUTIL=num2cell(ones(Nlots(TSTART),1));
        BASELAYER(cat(1,Lottype{:,2}))=1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%  Agricultural Layer   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % stream.Substream=35;
        rndstr.State=repeatstate{2,erun};
        farmmapfname='C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\data_files\FARMMAP_grid.mat';
        Sfarmmap=load_farmmap(farmmapfname);
        %         load FARMMAP_grid
        
        
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@    AGENTS    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        
        % <><><><><><><><><><><><><><>    FARMERS    <><><><><><><><><><><><><><><>
        LANDINFO(1,1:TMAX)=mat2cell(repmat(reshape(Sfarmmap.AGLAYER,NCELLS,1),...
            1,TMAX),NCELLS,ones(1,TMAX));
        sublandvalue=zeros(NCELLS,1);
        subpland=zeros(NCELLS,1);
        for nf=1:Nfarmers
            farmacres=find(LANDINFO{1,1}==nf);
            %     farmprod=ones(length(farmacres),1)*FARMPROD+PRODSTD*randn(1);
            %     farmcost=ones(length(farmacres),1)*FARMCOST+COSTSTD*randn(1);
            %     farmret=ones(length(farmacres),1)*normrnd(AVGFARMRETURN,STDFARMRETURN,1,1);
            
            %   %%%  AVGFARMRETURN=2486.3;
            %     [farmrow,farmcol]=ind2sub([NLENGTH NWIDTH],farmacres);
            farmmindist=min(coastdist(farmacres))*cell2mile;
            %     farmmindist=(min(farmcol)-5)*cell2mile;
            farmprod=ones(length(farmacres),1)*FARMPROD;
            farmcost=ones(length(farmacres),1)*FARMCOST;
            if farmmindist < 0.3
                farmret=ones(length(farmacres),1)*coastvalue(erun)*AVGFARMRETURN(erun)-...
                    cat(1,travelcost{farmacres});
            elseif farmmindist >=0.3 && farmmindist < 1
                farmret=ones(length(farmacres),1)*midvalue(erun)*AVGFARMRETURN(erun)-...
                    cat(1,travelcost{farmacres});
            elseif farmmindist >= 1
                farmret=ones(length(farmacres),1)*inlandvalue(erun)*AVGFARMRETURN(erun)-...
                    cat(1,travelcost{farmacres});
            end
            sublandvalue(farmacres)=farmret;
            subpland(farmacres)=farmret;
            if length(farmacres) < 3
                farmcntr=farmacres(1);
            else
                farmcntr=median(farmacres);
                if isempty(find(ismember(farmacres,farmcntr),1))==1
                    farmcntr=farmacres(find(farmacres > farmcntr,1,'first'));
                end
            end
            Farminfo{nf,1}=farmcntr;
            Farminfo{nf,2}=[farmacres farmprod farmcost farmret];
            %     Plandinfo(nf,1:TSTART+1)=mat2cell(reshape(repmat(reshape([farmret*ones(1,3) ...
            %         zeros(length(farmret),2)],length(farmret)*5,1),1,TSTART+1),...
            %         length(farmret),5,TSTART+1),length(farmret),5,ones(1,TSTART+1));
            farmretinfo(nf)=mean(farmret);
        end
        %[acres prod_costs value_acre]
        LANDINFO(2,1:TSTART+1)=mat2cell(repmat(sublandvalue,1,TSTART+1),...
            NCELLS,ones(1,TSTART+1));
        LANDINFO(3,1:TSTART+1)=mat2cell(repmat(sublandvalue,1,TSTART+1),...
            NCELLS,ones(1,TSTART+1));
        iNfarmers=unique(LANDINFO{1,TSTART});
        iNfarmers=iNfarmers(iNfarmers~=0);
        % subland=cat(2,LANDINFO{:,1});
        
        %@@@@@@@ Farmer Projection Models @@@@@@@@@@@@@@@@@@
        
        %Distance-Discounting models for farmers
        distcoeff=mincoeff+(maxcoeff-mincoeff)*rand(Nfarmers,NUMMODELDIST);
        %distance coefficient is in $1000/acre_distance
        
        landmodel = ceil(FARMNUMCLASS*rand(Nfarmers,NUMMODEL));
        %         for i = 1:FARMNUMCLASS
        %             strl = sprintf('landclass%d = find(landmodel == %d);',i,i);
        % %             eval(strl);
        %         end
        landclass1=find(landmodel == 1);
        landclass2=find(landmodel == 2);
        landclass3=find(landmodel == 3);
        landclass4=find(landmodel == 4);
        landclass5=find(landmodel == 5);
        landclass6=find(landmodel == 6);
        aa = zeros(Nfarmers,NUMMODEL);    %land models
        
        % mirror model
        aa(landclass1) = rand(1); % fraction that pred is away from 1/2 from mirror image
        
        % mean model
        aa(landclass2) = ceil(MAXMEANMODEL*rand(length(landclass2),1));
        
        %cycle model
        aa(landclass3) = ceil(MAXCYCLEMODEL*rand(length(landclass3),1));
        
        % projection model
        aa(landclass4) = ceil(2+((MAXPROJECT-1)-2)*rand(length(landclass4),1));
        
        % rescale model
        aa(landclass5) = 2*rand(length(landclass5),1);
        
        %local or regional trends
        aa(landclass6) = round(rand(length(landclass6),1));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%    FARMERS LEARNING MODULE    %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% FarmerModule_Coast_base
        % Price Projection models
        meanPLAND=zeros([],1);
        Planddist=zeros(Nfarmers,[]);
        farmerid=LANDINFO{1,TSTART};
        PLAND=LANDINFO{3,TSTART};
        subtc=cell2mat(travelcost);
        ifarmtrans=unique(farmerid(iurblist));      %co-opted from land market
        meanPLAND(1)=mean(PLAND(iurblist));
        FARMPROJ=zeros(Nfarmers,1);
        for rf=1:length(iNfarmers)
            ifarmcalc=find(farmerid == iNfarmers(rf));
            FARMPROJ(iNfarmers(rf))=mean(meanPLAND(1)-subtc(ifarmcalc));
            wtpland(iNfarmers(rf),1)=FARMPROJ(iNfarmers(rf));
        end
        learnwtaland(iNfarmers,1)=wtaland(iNfarmers,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        learnwtaland(iNfarmers,2)=learnwtaland(iNfarmers,1);
        successflag=1;
        successcount=0;
        successmark=zeros(1,[]);
        tlearn=1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%    Farmer Price Prediction    %%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while successcount <= 2
            tlearn=tlearn+1;
            %%% Linear, Random price signal
            PLAND(iurblist)=PLAND(iurblist)+200*rand(1);
            
            meanPLAND(tlearn)=mean(PLAND(iurblist));
            for rf=1:length(iNfarmers)
                ifarmcalc=find(farmerid == iNfarmers(rf));
                FARMPROJ(iNfarmers(rf))=mean(meanPLAND(tlearn)-subtc(ifarmcalc));
                wtpland(iNfarmers(rf),tlearn)=FARMPROJ(iNfarmers(rf));
            end
            
            %%%%% LAND MARKET %%%%%%%%
            Paskland(iNfarmers,tlearn)=learnwtaland(iNfarmers,tlearn);
            Plandproj(iNfarmers,tlearn)=mean([learnwtaland(iNfarmers,tlearn) ...
                wtpland(iNfarmers,tlearn)],2);
            
            % Calculate Distances to Transactions
            ntrans=length(ifarmtrans);   %number of land market trans actions, farms sold. In this module, for farmer learning, initial city farms sold held constatnt
            indtransdist=zeros(NLENGTH,NWIDTH,ntrans);
            transdist=zeros(NLENGTH,NWIDTH);
            Planddistproj=zeros(length(iNfarmers),ntrans);
            
            if ntrans > 0
                for nt=1:ntrans
                    transland=find(farmerid==ifarmtrans(nt));
                    for nf=1:length(iNfarmers)
                        rc=zeros([],1);
                        [transrow,transcol]=ind2sub([NLENGTH NWIDTH],Farminfo{ifarmtrans(nt),1});
                        [nfrow,nfcol]=ind2sub([NLENGTH NWIDTH],Farminfo{iNfarmers(nf),1});
                        avgtransdist=sqrt((nfrow-transrow)^2+(nfcol-transcol)^2);
                        coeffmark=(Paskland(iNfarmers(nf),tlearn)-mean(PLAND(transland)))/...
                            avgtransdist;
                        fitness(iNfarmers(nf),:,tlearn)=fitness(iNfarmers(nf),:,tlearn-1)+...
                            abs(distcoeff(iNfarmers(nf),:)-coeffmark);
                        fitsort=sort(fitness(iNfarmers(nf),:,tlearn),'ascend');
                        stratcount=1;
                        isurvive=zeros(length(nextgen),1);
                        for x=1:nextgen
                            numstrat=find(fitness(iNfarmers(nf),:,tlearn)==fitsort(x));
                            if length(numstrat) > 1
                                isurvive(x,1)=numstrat(stratcount);
                                stratcount=stratcount+1;
                            else
                                isurvive(x,1)=numstrat(stratcount);
                            end
                            if stratcount > length(numstrat)
                                stratcount=1;
                            end
                        end
                        distcoeff(iNfarmers(nf),1:nextgen)=distcoeff(iNfarmers(nf),isurvive);
                        for xx=1:nextgen
                            irc=xx*recombo/nextgen;
                            rc(length(rc)+1:irc,1)=distcoeff(iNfarmers(nf),xx)+rand(1,4);
                        end
                        distcoeff(iNfarmers(nf),nextgen+1:nextgen+recombo)=rc';
                        distcoeff(iNfarmers(nf),nextgen+recombo+1:NUMMODELDIST)=...
                            mincoeff+(maxcoeff-mincoeff)*rand(1,NUMMODELDIST-...
                            (nextgen+recombo));
                    end
                    Planddistproj(iNfarmers,nt)=distcoeff(iNfarmers,1).*avgtransdist+...
                        mean(PLAND(transland));
                end
                for nnw=1:NWIDTH
                    for nnl=1:NLENGTH
                        transdist(nnl,nnw)=min(indtransdist(nnl,nnw,1:ntrans));
                    end
                end
                Planddist(iNfarmers,tlearn)=mean(Planddistproj(iNfarmers,:),2);
            else
                Planddist(iNfarmers,tlearn)=mean([learnwtaland(iNfarmers,tlearn) wtpland(iNfarmers,tlearn)],2);
            end
            %%% Land Price projections for tlearn+1
            for nf=1:length(iNfarmers)
                
                ilandclass1=find(landmodel(iNfarmers(nf),:)==1);
                ilandclass2=find(landmodel(iNfarmers(nf),:)==2);
                ilandclass3=find(landmodel(iNfarmers(nf),:)==3);
                ilandclass4=find(landmodel(iNfarmers(nf),:)==4);
                ilandclass5=find(landmodel(iNfarmers(nf),:)==5);
                ilandclass6=find(landmodel(iNfarmers(nf),:)==6);
                
                % mimic models
                learnlandproj(iNfarmers(nf),ilandclass1) = Plandproj(iNfarmers(nf),tlearn)+(1-aa...
                    (iNfarmers(nf),ilandclass1)).*(0.5*Plandproj(iNfarmers(nf),tlearn)-...
                    (Plandproj(iNfarmers(nf),tlearn)-Plandproj(iNfarmers(nf),tlearn-1)));
                
                % mean model
                for jl = 1:length(ilandclass2)
                    learnlandproj(iNfarmers(nf),ilandclass2(jl)) = mean(Plandproj(iNfarmers(nf),...
                        tlearn:-1:max((tlearn-aa(iNfarmers(nf),ilandclass2(jl))),1)));
                end
                
                
                %cycle model
                learnlandproj(iNfarmers(nf),ilandclass3) = Plandproj(iNfarmers(nf),...
                    max(tlearn-round(aa(iNfarmers(nf),ilandclass3)),1));
                
                % projection model
                for jl = 1:length(ilandclass4)
                    %Nonlinear Forecast
                    warning off
                    indata=Plandproj(iNfarmers(nf),tlearn-min(tlearn-1,(1+aa(iNfarmers(nf),ilandclass4(jl)))):tlearn);
                    %                     pcoef=polyfit(1:length(indata),indata,2);
                    %                     pline=pcoef(1).*(1:length(indata)+1).^2+pcoef(2).*(1:...
                    %                         length(indata)+1)+pcoef(3);
                    pcoef=polyfit(1:length(indata),indata,1);
                    pline=pcoef(1).*(1:length(indata)+1)+pcoef(2);
                    learnlandproj(iNfarmers(nf),ilandclass4(jl))=pline(length(pline));
                    warning on
                end
                
                
                % rescale model
                learnlandproj(iNfarmers(nf),ilandclass5) = aa(iNfarmers(nf),...
                    ilandclass5)*Plandproj(iNfarmers(nf),tlearn);
                
                % local(1) or regional(0) trends
                ilandreg=(aa(iNfarmers(nf),ilandclass6)==0);
                ilandlocal=(aa(iNfarmers(nf),ilandclass6)==1);
                % Local: just spatially discounted
                learnlandproj(iNfarmers(nf),ilandclass6(ilandlocal)) = Planddist...
                    (iNfarmers(nf),tlearn);
                if isempty(iNfarmers)==1
                    break
                end
                % Regional: density-dependent
                learnlandproj(iNfarmers(nf),ilandclass6(ilandreg)) = Planddist...
                    (iNfarmers(nf),tlearn).*(1+1/length(iNfarmers));
                
                if tlearn > 2
                    learnlanderror(iNfarmers(nf),:) = (1-learnDELTA)*learnlanderror(iNfarmers(nf),:)+...
                        learnDELTA*abs(Plandproj(iNfarmers(nf),tlearn)-learnlandprojSAVE(iNfarmers(nf),:));
                    % Use model that predicted this period's price best
                    [landbest,ilandbest] = min(learnlanderror(iNfarmers(nf),:),[],2);
                else
                    [landbest,ilandbest] = min(learnlanderror(iNfarmers(nf),:),[],2);
                end
                subfarminfo=Farminfo{iNfarmers(nf),2};
                learnwtaland(iNfarmers(nf),tlearn+1)=max(learnlandproj(iNfarmers(nf),ilandbest),...
                    mean(subfarminfo(:,4)));
                %         learnwtaland(iNfarmers(nf),tlearn+1)=max(learnlandproj(iNfarmers(nf),ilandbest),...
                %             Farmstats(iNfarmers(nf),5,1));
                
                learnlandprojSAVE(iNfarmers(nf),:)=learnlandproj(iNfarmers(nf),:);
                learnlandbestSAVE(iNfarmers(nf),tlearn) = landbest;
                ilearnlandbestSAVE(iNfarmers(nf),tlearn) = ilandbest;
                learnlandprojbestSAVE(iNfarmers(nf),tlearn+1) = learnlandproj(iNfarmers(nf),ilandbest);
                learnlandmodelbestSAVE(iNfarmers(nf),tlearn) = landmodel(iNfarmers(nf),ilandbest);
            end
            
            %%%% End tlearn loop
            
            %         figure(1)
            %         hist(landmodelbestSAVE(:,tlearn),1:6)
            
            landprojdiff(:,tlearn)=Plandproj(:,tlearn)-learnlandprojbestSAVE(:,tlearn);
            landpctdiff(:,tlearn)=abs(landprojdiff(:,tlearn))./Plandproj(:,tlearn);
            landpricevar(:,tlearn)=abs(diff(Plandproj(:,tlearn-1:tlearn),1,2)./Plandproj(:,tlearn));
            if tlearn > TSTART+1
                %             successflag=(pctdiff(:,tlearn) >= pricevar(:,tlearn));
                %             successflag=(sum(pctdiff(:,tlearn)) >= sum(pricevar(:,tlearn)))+...
                %                 (length(find(pctdiff(:,tlearn) <= pricevar(:,tlearn))==1) <= 0.9*Nfarmers);
                successflag=(length(find(landpctdiff(:,tlearn) < 0.10)==1) <= 0.9*Nfarmers);
                if successflag==1
                    successmark(length(successmark)+1)=0;
                    successcount=0;
                elseif successflag==0
                    successmark(length(successmark)+1)=1;
                    successcount=successcount+1;
                end
            end
        end
        
        % Feed model success information to TSTART
        
        fitnessSAVE=fitness(:,:,tlearn);
        
        
        %%%%%%%%%%%%%%%   ADD LEARNED INFO   %%%%%%%%%%%%%%%
        
        % %Implicitly added:
        % landerror
        % distcoeff
        fitness(:,:,TSTART)=fitnessSAVE;
        landerror=learnlanderror;
        wtaland(iNfarmers,1:TSTART+1)=max(learnwtaland(iNfarmers,tlearn-TSTART:tlearn),...
            repmat(farmretinfo(iNfarmers),1,length(tlearn-TSTART:tlearn)));
        Paskland(iNfarmers,1:TSTART+1)=Paskland(iNfarmers,tlearn-TSTART:tlearn);
        
        landbestSAVE(iNfarmers,1:TSTART)=learnlandbestSAVE(iNfarmers,tlearn-TSTART:tlearn-1);
        ilandbestSAVE(iNfarmers,1:TSTART)=ilearnlandbestSAVE(iNfarmers,tlearn-TSTART:tlearn-1);
        landprojbestSAVE(iNfarmers,1:TSTART)=learnlandprojbestSAVE(iNfarmers,tlearn-TSTART:...
            tlearn-1);
        landmodelbestSAVE(iNfarmers,1:TSTART)=learnlandmodelbestSAVE(iNfarmers,...
            tlearn-TSTART:tlearn-1);
        Plandproj(iNfarmers,1:TSTART)=learnlandprojbestSAVE(iNfarmers,tlearn-TSTART:tlearn-1);
        
        % transacres=cat(1,Farminfo{ifarmtrans,2});
        subfarminfo=LANDINFO{1,TSTART};
        subfarminfo(BASELAYER == 1)=0;
        LANDINFO{1,TSTART}=subfarminfo;
        % transacres=find(ismember(LANDINFO{1,TSTART},ifarmtrans)==1);
        % ivacland=ismember(transacres,find(BASELAYER==0));
        % vacland(1,1:TSTART)=mat2cell(repmat(transacres(ivacland,1),1,TSTART),...
        %     length(find(ivacland==1)),ones(1,TSTART));
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%    Developers    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % stream.Substream=4;
        rndstr.State=repeatstate{3,erun};
        %%%% Developer's Population Prediction Models %%%%%%%
        
        classagentmodel = ceil(POPNUMCLASS*rand(Ndevelopers,NUMMODEL));
        
%         for i = 1:POPNUMCLASS
%             str = sprintf('indclass%d = find(classagentmodel == %d);',i,i);
%             eval(str);
%         end
        indclass1=find(classagentmodel == 1);
        indclass2=find(classagentmodel == 2);
        indclass3=find(classagentmodel == 3);
        indclass4=find(classagentmodel == 4);
        indclass5=find(classagentmodel == 5);
        dd = zeros(Ndevelopers,NUMMODEL);
        
        % mirror model
        dd(indclass1) = 0.40+(0.60-0.40)*rand(1,length(indclass1)); % fraction that pred is away from 1/2 from mirror image
        
        % mean model
        dd(indclass2) = ceil(MAXMEANMODEL*rand(length(indclass2),1));
        
        %cycle model
        dd(indclass3) = ceil(MAXCYCLEMODEL*rand(length(indclass3),1));
        
        % projection model
        dd(indclass4) = 1+ceil((MAXPROJECT-1)*rand(length(indclass4),1));
        
        % rescale model
        dd(indclass5) = 2*rand(length(indclass5),1);
        
        nproj = zeros(Ndevelopers,NUMMODEL);
        errorsq = zeros(Ndevelopers,NUMMODEL);
        errorsq(1:Ndevelopers,1:NUMMODEL) = rand(Ndevelopers,NUMMODEL);
        
        if POPNUMCLASS >= 4
            xxx = zeros(length(indclass4),MAXPROJECT);
            yyy = zeros(length(indclass4),MAXPROJECT);
            for j = 1:length(indclass4)
                xxx(j,1:dd(indclass4(j))) = ((-dd(indclass4(j))+1):0);
            end
            sumx = sum(xxx,2);
            sumx2 = sum(xxx.^2,2);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%@@@@@@@    Broker Projection Models    @@@@@@@@@@@@@@@@@@%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        bbfull = zeros(HT,NUMMODEL,Nbrokers);    %broker models
        brokermodel = ceil(BROKERNUMCLASS*rand(HT,NUMMODEL,Nbrokers));
        
        for nb=1:Nbrokers
            %             for i = 1:BROKERNUMCLASS
            %                 strb = sprintf('brokerclass%d = find(brokermodel(:,:,nb) == %d);',i,i);
            %                 eval(strb);
            %             end
            bb=zeros(HT,NUMMODEL);
            
            %         bb=bbfull(:,:,nb);
            brokerclass1=(brokermodel(:,:,nb) == 1);
            brokerclass2=find(brokermodel(:,:,nb) == 2);
            brokerclass3=find(brokermodel(:,:,nb) == 3);
            brokerclass4=find(brokermodel(:,:,nb) == 4);
            brokerclass5=find(brokermodel(:,:,nb) == 5);
            brokerclass6=find(brokermodel(:,:,nb) == 6);
            
            % mirror model
            bb(brokerclass1) = rand(1); % fraction that pred is away from 1/2 from mirror image
            
            % mean model
            bb(brokerclass2) = ceil(MAXMEANMODEL*rand(length(brokerclass2),1));
            
            %cycle model
            bb(brokerclass3) = ceil(MAXCYCLEMODEL*rand(length(brokerclass3),1));
            
            % projection model
            bb(brokerclass4) = ceil(2+((MAXPROJECT-1)-2)*rand(length(brokerclass4),1));
            
            % rescale model
            bb(brokerclass5) = 2*rand(length(brokerclass5),1);
            
            %local or regional trends
            bb(brokerclass6) = ceil(rand(length(brokerclass6),1));
        end
        
        bbfull(:,:,nb)=bb;
    
%         clearvars brokerclass1 brokerclass2 brokerclass3 brokerclass4 brokerclass5 brokerclass6 bb
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%    Broker's Learning Module    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% BrokersModule_Coast_0v3
        for lt=1:HT
            pcoeffs(lt,:)=polyfit(testtime,testrents(lt,:),1);
        end
        
        % rentinfo=pcoeffs(1:12,1)*(1:TMAX)+pcoeffs(1:12,2)*ones(1,TMAX)+(500*randn(12,TMAX));
        rentinfo=pcoeffs(1:HT,1)*(1:TMAX)+(pcoeffs(1:HT,2)*...
            ones(1,TMAX).*(0.85+(1.15-0.85)*rand(HT,TMAX)));
        
        successflag=1;
        successcount=0;
        successmark=zeros(1,[]);
        tlearn=TSTART;
        
        while successcount <= 2
            tlearn=tlearn+1;
            
            for nb=1:Nbrokers
                
                bb=bbfull(:,:,nb);
                for lt=1:HT
%                     for i = 1:BROKERNUMCLASS
%                         %         strb = sprintf('brokerclass%d = find(brokermodel(:,:,nb) == %d);',i,i);
%                         strb = sprintf('brokerclass%d = find(brokermodel(lt,:,nb) == %d);',i,i);
%                         eval(strb);
%                     end
                    brokerclass1=find(brokermodel(lt,:,nb)==1);
                    brokerclass2=find(brokermodel(lt,:,nb)==2);
                    brokerclass3=find(brokermodel(lt,:,nb)==3);
                    brokerclass4=find(brokermodel(lt,:,nb)==4);
                    brokerclass5=find(brokermodel(lt,:,nb)==5);
                    brokerclass6=find(brokermodel(lt,:,nb)==6);
                    
                    % mimic models
                    learnbproj(lt,brokerclass1) = rentinfo(lt,tlearn)+(1-bb...
                        (lt,brokerclass1)).*(0.5*rentinfo(lt,tlearn)-...
                        (rentinfo(lt,tlearn)-rentinfo(lt,tlearn-1)));
                    
                    % mean model
                    for jl = 1:length(brokerclass2)
                        learnbproj(lt,brokerclass2(jl)) = mean(rentinfo(lt,...
                            tlearn:-1:(tlearn-bb(lt,brokerclass2(jl)))));
                    end
                    
                    %cycle model
                    learnbproj(lt,brokerclass3) = rentinfo(lt,tlearn-...
                        round(max(1,bb(lt,brokerclass3))));
                    
                    % projection model
                    for jl = 1:length(brokerclass4)
                        %Nonlinear Forecast
                        indata=rentinfo(lt,tlearn-(1+bb(lt,brokerclass4(jl))):tlearn);
                        subindata=reshape(indata,1,length(indata));
                        pcoef=polyfit(1:length(indata),subindata,1);
                        pline=pcoef(1).*(1:length(indata)+1)+pcoef(2);
                        %                         pcoef=polyfit(1:length(indata),subindata,2);
                        %                         pline=pcoef(1).*(1:length(indata)+1).^2+pcoef(2).*(1:...
                        %                             length(indata)+1)+pcoef(3);
                        learnbproj(lt,brokerclass4(jl))=pline(length(pline));
                    end
                    
                    % rescale model
                    learnbproj(lt,brokerclass5) = bb(lt,brokerclass5)*rentinfo(lt,tlearn);
                    
                    [brows,bcols]=ind2sub([nbrokerlong nbrokerwide],nb);
                    brnei=(bcols+1)*nbrokerlong-(nbrokerwide-brows);
                    blnei=(bcols-1)*nbrokerlong-(nbrokerwide-brows);
                    bupnei=bcols*nbrokerlong-(nbrokerwide-(brows-1));
                    bdnnei=bcols*nbrokerlong-(nbrokerwide-(brows+1));
                    ibnei=[brnei blnei bupnei bdnnei];
                    realbnei=find(minibmap==brnei | minibmap==blnei | ...
                        minibmap==bupnei | minibmap==bdnnei);
                    learnbproj(lt,brokerclass6) = (bb(lt,brokerclass6)+...
                        rand(1,length(brokerclass6)))*mean(rentinfo(lt,tlearn));
                    
                    
                    learnbrokererror(lt,:,nb,tlearn) = (1-learnDELTA)*learnbrokererror(lt,:,nb,tlearn-1)+...
                        learnDELTA*abs(rentinfo(lt,tlearn)-learnbproj(lt,:));
                    learnabserror(lt,:,nb,tlearn)=rentinfo(lt,tlearn)-learnbproj(lt,:);
                    [brokerbest,ibrokerbest] = min(learnbrokererror(lt,:,nb,tlearn),[],2);
                    if tlearn > TSTART+1
                        difflearnerror(lt,:,nb,tlearn)=learnbrokererror(lt,:,nb,tlearn)-...
                            learnbrokererror(lt,:,nb,tlearn-1);
                        learnbestdiffSAVE(nb,lt,tlearn)=difflearnerror(lt,ibrokerbest,nb,tlearn);
                        learnbestabsSAVE(nb,lt,tlearn)=learnabserror(lt,ibrokerbest,nb,tlearn);
                    else
                        difflearnerror(lt,:,nb,tlearn)=0;
                    end
                    learnbrokerbestSAVE(nb,lt,tlearn) = brokerbest';
                    ilearnbrokerbestSAVE(nb,lt,tlearn) = ibrokerbest';
                    learnbrokerprojSAVE(nb,lt,tlearn) = learnbproj(lt,ibrokerbest);
                    learnbrokermodelSAVE(nb,lt,tlearn) = brokermodel(lt,ibrokerbest,nb);
                    
                end
                brokerproj(:,:,nb)=learnbproj;
                
                brkrprojdiff(:,nb,tlearn)=rentinfo(:,tlearn)-learnbrokerprojSAVE(nb,:,tlearn)';
                brkrpctdiff(:,nb,tlearn)=abs(brkrprojdiff(:,nb,tlearn))./rentinfo(:,tlearn);
                brkrpricevar(:,nb,tlearn)=abs(diff(rentinfo(:,tlearn-1:tlearn),1,2)./...
                    rentinfo(:,tlearn));
            end
            if tlearn > TSTART+TSTART+1
                meansuccess=mean(brkrpctdiff(:,:,tlearn));
                successflag=(length(find(meansuccess < 0.10)==1) <= 0.9*Nbrokers);
                if successflag==1
                    successmark(length(successmark)+1)=0;
                    successcount=0;
                elseif successflag==0
                    successmark(length(successmark)+1)=1;
                    successcount=successcount+1;
                end
            end
        end
        %%%%%%%%%%%%%%%   ADD LEARNED INFO   %%%%%%%%%%%%%%%
        brokererror(:,:,:,1:TSTART)=learnbrokererror(:,:,:,tlearn-9:tlearn);
        brokerbestdiffSAVE(:,:,1:TSTART)=learnbestdiffSAVE(:,:,max(12,tlearn-9):tlearn);
        brokerabserror(:,:,:,1:TSTART)=learnabserror(:,:,:,tlearn-9:tlearn);
        brokerbestabsSAVE(:,:,1:TSTART)=learnbestabsSAVE(:,:,max(12,tlearn-9):tlearn);
        
        distfname='C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\data_files\master_dist.mat';
        Sdist=load_distmat(distfname);
%         load master_dist

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%    Model Runs    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% CHALMS_Coast_batch
        %%% initial vacant land
        ivac=(VACLAND ~= 0);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%    Consumers    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assign random income and houseing preferences to consumers
        housepref=zeros(Nconstart,1);
        subincome=round(exp(lognrnd(parmhat(1),parmhat(2),Nconstart,1)));
        % %CONINFO=[income,searchtime,consumer_good,housesize,lotsize,proximity,subrisk,occ/vac,utility]
        CONINFO(:,1)=num2cell(min(max(subincome,minwage),maxwage));
        CONINFO(:,2)=num2cell(TSTART+ceil(searchtimemin+...
            (searchtimemax-searchtimemin)*rand(Nconstart,1)));
        income=cat(1,CONINFO{:,1});
        income3=(income >= lowminwage & income <= lowmaxwage);     %income 1 2 3 = hi mid low wages
        income2=(income >= midminwage & income <= midmaxwage);
        income1=(income >= himinwage & income <= himaxwage); % Different proportions of income spent on housing depending on income level
        housepref(income1)=HIBETA(1)+(HIBETA(2)-HIBETA(1))*rand(length(find(income1==1)),1);
        housepref(income2)=MIDBETA(1)+(MIDBETA(2)-MIDBETA(1))*rand(length(find(income2==1)),1);
        housepref(income3)=LOWBETA(1)+(LOWBETA(2)-LOWBETA(1))*rand(length(find(income3==1)),1);
        
        CONINFO(:,3)=num2cell(1-housepref);
        CONINFO(:,6)=num2cell((ampref_min(erun)+(ampref_max(erun)-ampref_min(erun))*...
            rand(length(housepref),1)).*housepref);
        CONINFO(:,4)=num2cell(rand(length(housepref),1).*...
            (housepref-cat(1,CONINFO{:,6}))/2);
        CONINFO(:,5)=num2cell(housepref-(cat(1,CONINFO{:,4})+cat(1,CONINFO{:,6})));
        
        % CONINFO(:,4)=num2cell((0.1+(0.9-0.1)*rand(length(housepref),1)).*housepref);
        % CONINFO(:,5)=num2cell(rand(length(housepref),1).*...
        %     (housepref-cat(1,CONINFO{:,4}))/2);
        % CONINFO(:,6)=num2cell(housepref-(cat(1,CONINFO{:,4})+cat(1,CONINFO{:,5})));
        CONINFO(:,7)=num2cell(ones(length(CONINFO(:,1)),1));
        CONINFO(:,8)=num2cell(zeros(length(CONINFO(:,1)),1));
        
        %%% Subjective risk perception
        for tt=1:TSTART
            damcoef(:,tt)=mat2cell(ones(Nconstart,Nconstart),Nconstart,ones(Nconstart,1));
        end
        % stream.Substream=mrun;
        rndstr.State=savedState;
        
        %%% Spin-up housing market, developer learns pricing
        rentmodelfname='C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\data_files\rentmodel.mat';
        Srent=load_rentmodel(rentmodelfname);       
%         load rentmodel

        
        %%% Initialize based on calibrated model
        Paskhouse(1:Nlots(1),1)=predict(Srent.rentmdl,[cat(1,Lottype{:,3}) cat(1,Lottype{:,8})...
            cat(1,Lottype{:,7}) cat(1,[CONINFO{1:length(CONINFO),1} ...
            CONINFO{1:Nlots(1)-length(CONINFO),1}])']);
        
        deltadiff=zeros(1,[]);
        utildiff=zeros(1,TMAX);
        meandiffprice=5000;
        utildiff(1:TSTART)=1000;
        diffcheck=0;
        stillvaccheck=0;
        meanprice=zeros(HT,[]);
        varprice=zeros(HT,TMAX);
        iterhmc=zeros(1,[]);
        diffvac=zeros(1,[]);
        diffprice=zeros(1,[]);
        iter=0;
        deltaiter=1;
        deltafac=1;
        diffflag=0;
        killflag=0;
        vaccheck=Nconstart;
        %%%%%% Learning Period %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        while abs(deltafac) > 0.001
            initialdiff=meandiffprice;
            initialutildiff=utildiff(TSTART);
            
            
            %%% HouseMarketInitial_coast_batch
            %%%%%%%%%%%%%% House Market Initial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            wtpcon=zeros(length(CONINFO(:,1)),Nlots(TSTART));
            wtpconstar=zeros(length(CONINFO(:,1)),Nlots(TSTART));
            wtbcon=zeros(length(CONINFO(:,1)),Nlots(TSTART));
            Rn=zeros(length(CONINFO(:,1)),Nlots(TSTART));
            Phousebid=zeros(length(CONINFO(:,1)),Nlots(TSTART));
            subnhouselook=zeros(length(CONINFO(:,1)),Nlots(TSTART));
            notherbuyers=zeros(length(CONINFO(:,1)),Nlots(TSTART));
            iotherbuyers=zeros(length(CONINFO(:,1)),Nlots(TSTART));
            housemp=zeros(length(CONINFO(:,1)),1);
            EUmit=zeros(length(CONINFO(:,1)),Nlots(TSTART));
            EUnomit=zeros(length(CONINFO(:,1)),Nlots(TSTART));
            U=zeros(length(CONINFO(:,1)),Nlots(TSTART));
            Unorm=zeros(length(CONINFO(:,1)),Nlots(TSTART));
            maxUset=zeros(length(CONINFO(:,1)),LTYPE*HTYPE);
            exptCmit=zeros(length(CONINFO(:,1)),Nlots(TSTART));
            exptCnomit=zeros(length(CONINFO(:,1)),Nlots(TSTART));
            mitchoice=zeros(length(CONINFO(:,1)),Nlots(TSTART));
            exptcost=zeros(2,Nlots(TSTART),length(CONINFO(:,1)));
            % %Lottype=[id,location index,lotsize,housesize,ltype,ccost,amlevel,travelcost,buildtime,brokerid]
            % %lotchoice=[id,location index,ltype,occ/vac,consumer id,residence time,sell price,mitchoice]
            % %CONINFO=[income,searchtime,consumer_good,housesize,lotsize,proximity,subrisk,occ/vac,utility]
            
            % average damage across storm categories
            Cdam(1:Nlots(TSTART),TSTART)=num2cell(Paskhouse.*...
                (housedam(cat(1,lotchoice{:,2})))./100);
%             Cdam(1:Nlots(TSTART),TSTART)=num2cell(mean((Paskhouse*housedam).*...
%                 repmat(mean(Psevere,1),Nlots(TSTART),1).*...
%                 repmat(cat(1,Pflood{cat(1,lotchoice{:,2})}),1,5),2));
            for c=1:length(CONINFO(:,1))
                wtpcon(c,:)=(CONINFO{c,1}-cat(1,travelcost{cat(1,lotchoice{:,2})})-...
                    cat(1,Cdam{1:Nlots(TSTART),TSTART}))*(CONINFO{c,4}+...
                    CONINFO{c,5}+CONINFO{c,6});
                EUnomit(c,:)=sum(Psevere).*...
                    (((CONINFO{c,1}-cat(1,travelcost{cat(1,lotchoice{:,2})})-...
                    Paskhouse-cat(1,Cdam{1:Nlots(TSTART),TSTART})).^CONINFO{c,3}).*...
                    (cat(1,Lottype{:,4}).^CONINFO{c,4}).*(cat(1,Lottype{:,3}).^...
                    CONINFO{c,5}).*(cat(1,Lottype{:,7}).^CONINFO{c,6}))+...
                    (1-sum(Psevere)).*...
                    (((CONINFO{c,1}-cat(1,travelcost{cat(1,lotchoice{:,2})})-...
                    Paskhouse).^CONINFO{c,3}).*...
                    (cat(1,Lottype{:,4}).^CONINFO{c,4}).*(cat(1,Lottype{:,3}).^...
                    CONINFO{c,5}).*(cat(1,Lottype{:,7}).^CONINFO{c,6}));
%                 U(c,:)=((CONINFO{c,1}-cat(1,travelcost{cat(1,lotchoice{:,2})})-...
%                     Paskhouse).^CONINFO{c,3}).*(cat(1,Lottype{:,4}).^CONINFO{c,4}).*...
%                     (cat(1,Lottype{:,3}).^CONINFO{c,5}).*(cat(1,Lottype{:,7}).^CONINFO{c,6});
                %%%%%%%% Expected Costs %%%%%%%%%
                %     subdamcoef=damcoef{c,TSTART};
                %     exptCmit(c,:)=subdamcoef.*cat(1,Pflood{cat(1,lotchoice{:,2})}).*...
                %         (Cmit(erun)+miteff(erun)*Cdam(erun)*Paskhouse)+(1-subdamcoef).*...
                %         (1-cat(1,Pflood{cat(1,lotchoice{:,2})}))*Cmit(erun);
                %     exptCnomit(c,:)=subdamcoef.*cat(1,Pflood{cat(1,lotchoice{:,2})}).*(Cdam(erun)*Paskhouse);
                %     exptcost(:,:,c)=[exptCmit(c,:); exptCnomit(c,:)];
                %     wtpcon(c,:)=(CONINFO{c,1}-cat(1,travelcost{cat(1,lotchoice{:,2})}))*(CONINFO{c,4}+...
                %         CONINFO{c,5}+CONINFO{c,6});
                %     EUmit(c,:)=(max(CONINFO{c,1}-cat(1,travelcost{cat(1,lotchoice{:,2})})-...
                %         Paskhouse-exptCmit(c,:)',0).^CONINFO{c,3}).*...
                %         (cat(1,Lottype{:,4}).^CONINFO{c,4}).*(cat(1,Lottype{:,3}).^CONINFO{c,5}).*...
                %         (cat(1,Lottype{:,7}).^CONINFO{c,6});
                %     EUnomit(c,:)=(max(CONINFO{c,1}-cat(1,travelcost{cat(1,lotchoice{:,2})})-...
                %         Paskhouse-exptCnomit(c,:)',0).^CONINFO{c,3}).*...
                %         (cat(1,Lottype{:,4}).^CONINFO{c,4}).*(cat(1,Lottype{:,3}).^CONINFO{c,5}).*...
                %         (cat(1,Lottype{:,7}).^CONINFO{c,6});
                %     [Umax,mitcheck]=max([EUmit(c,:); EUnomit(c,:)]);
                    U(c,:)=EUnomit(c,:);
                %     mitchoice(c,:)=2*ones(1,length(EUnomit(c,:)));
                % %     U(c,:)=Umax;
                % %     mitchoice(c,:)=mitcheck;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %     ihousein=(cat(1,BIDLEVEL{:}).*Paskhouse < wtpcon(c,:)' & U(c,:)' > ...
                %         CONINFO{c,3}*max(U(c,:)));
                ihousein=cat(1,BIDLEVEL{:}).*Paskhouse < wtpcon(c,:)';
                ihouseout=find(ihousein == 0);
                
                Unorm(c,ihousein)=U(c,ihousein)./max(U(c,ihousein));
                Rn(c,ihousein)=Paskhouse(ihousein).*Unorm(c,ihousein)';
                wtbcon(c,ihousein)=min(wtpcon(c,ihousein)-(Paskhouse(ihousein)'-...
                    Rn(c,ihousein)),wtpcon(c,ihousein));
                subnhouselook(c,:)=(ihousein==1);
            end
            nhouselook=(subnhouselook == 1);
            for nl=1:Nlots(TSTART)
                notherbuyers(nhouselook(:,nl),nl)=find(nhouselook(:,nl)==1);
            end
            
            for c=1:length(CONINFO(:,1))
                if isempty(find(nhouselook(c,:),1))==1
                    continue
                else
                    nhouses=length(find(nhouselook(c,:)==1));
                    subbuyers=unique(notherbuyers(:,nhouselook(c,:)));
                    subbuyers=subbuyers(subbuyers~=0);
                    nbuyers=length(subbuyers);
                    housemp(c)=0.5*(nbuyers-nhouses)/(nbuyers+nhouses);
                    
                    if housemp(c) >= 0
                        Phousebid(c,nhouselook(c,:))=min(max(Rn(c,nhouselook(c,:))+...
                            (wtbcon(c,nhouselook(c,:))-Rn(c,nhouselook(c,:)))*...
                            housemp(c),Rn(c,nhouselook(c,:))),wtbcon(c,nhouselook(c,:)));
                    elseif housemp(c) < 0
                        Phousebid(c,nhouselook(c,:))=min(Rn(c,nhouselook(c,:))+...
                            Rn(c,nhouselook(c,:)).*(1./(wtbcon(c,nhouselook(c,:))-...
                            Rn(c,nhouselook(c,:))))*housemp(c),wtbcon(c,nhouselook(c,:)));
                    end
                end
            end
            
            avghousemp(1:TSTART)=mean(housemp);
            openhouse=Nlots(TSTART);
            % con2lot=zeros(Nlots(TSTART),4);     %[WinBid conid lotid restime]
            subPhousebid=Phousebid;
            subU=U;
            for nl=1:Nlots(TSTART)
                iunderbid=(subPhousebid(:,nl) <= 0);
                subPhousebid(iunderbid,nl)=0;
            end
            while openhouse > 0
                if isempty(find(subPhousebid > 0,1))==1
                    break
                end
                istillopen=find(cat(1,lotchoice{:,4})==0);
                wincon=zeros(1,Nlots(TSTART));
                
                [maxbid,imaxbid]=max(subPhousebid,[],1);
                
                for nl=1:Nlots(TSTART)
                    if maxbid(nl) <= 0
                        continue
                    end
                    %check for multiple consumers with same, highest bid
                    iwincon=find(subPhousebid(:,nl)==maxbid(nl));
                    if length(iwincon) > 1
                        icon=ceil(length(iwincon)*rand(1));
                        wincon(nl)=iwincon(icon);
                    else
                        wincon(nl)=iwincon;
                    end
                end
                conset=unique(wincon);      %Highest bidders at the moment
                conset=conset(conset~=0);
                randorder=randperm(length(conset));
                conset=conset(randorder);
                
                for cs=1:length(conset)
                    ilotmatch=find(wincon==conset(cs));
                    %        subexptcost=exptcost(:,:,conset(cs));
                    %        iexptcost=sub2ind(size(subexptcost),mitchoice(conset(cs),ilotmatch),ilotmatch);
                    %        uset=(CONINFO{conset(cs),1}-cat(1,travelcost{cat(1,lotchoice{ilotmatch,2})})-...
                    %            subPhousebid(conset(cs),ilotmatch)'-subexptcost(iexptcost)'.^CONINFO{conset(cs),3}).*...
                    %            (cat(1,Lottype{ilotmatch,4}).^CONINFO{conset(cs),4}).*...
                    %            (cat(1,Lottype{ilotmatch,3}).^CONINFO{conset(cs),5}).*...
                    %            (cat(1,Lottype{ilotmatch,7}).^CONINFO{conset(cs),6});
                    uset=((CONINFO{conset(cs),1}-cat(1,travelcost{cat(1,lotchoice{ilotmatch,2})})-...
                        subPhousebid(conset(cs),ilotmatch)').^CONINFO{conset(cs),3}).*...
                        (cat(1,Lottype{ilotmatch,4}).^CONINFO{conset(cs),4}).*...
                        (cat(1,Lottype{ilotmatch,3}).^CONINFO{conset(cs),5}).*...
                        (cat(1,Lottype{ilotmatch,7}).^CONINFO{conset(cs),6});
                    ilotid=find(uset==max(uset));
                    
                    %Match to lot with highest bid
                    
                    if length(ilotid) > 1
                        ipick=ceil(length(ilotid)*rand(1));
                        lotid=ilotmatch(ilotid(ipick));
                    else
                        ipick=1;
                        lotid=ilotmatch(ilotid);
                    end
                    conid=conset(cs);
                    CONINFO{conid,8}=1;
                    CONINFO{conid,9}=uset(ilotid(ipick));
                    %        CONINFO{conid,9}=uset(conid,lotid);
                    lotchoice{lotid,4}=1;
                    lotchoice{lotid,5}=conid;
                    lotchoice{lotid,6}=max(ceil(TSTART+avgrestime/2+normrnd(avgrestime/...
                        2,stdrestime/2,1,1)),TSTART+1);
                    lotchoice{lotid,7}=subPhousebid(conid,lotid);
                    %        lotchoice{lotid,8}=mitchoice(conid,lotid);
                    lotchoice{lotid,8}=1;
                    
                    %        con2lot(lotid,1)=subPhousebid(conid,lotid);
                    %        con2lot(lotid,2)=conid;
                    %        con2lot(lotid,3)=lotid;
                    %        lotchoice(lotid,7)=1;
                    %        Lottype(Lottype(:,1)==lotid,7)=1;
                    %        cellinfo(:,7)=Lottype(icells,7);
                    %        con2lot(lotid,4)=max(ceil(TSTART+avgrestime/2+normrnd(avgrestime/2,stdrestime/2,1,1)),TSTART+1);
                    %        RESTIME(Lottype(Lottype(:,1)==lotid,2))=con2lot(lotid,4);
                    subPhousebid(conid,:)=0;
                    subPhousebid(:,lotid)=0;
                    subU(conid,:)=0;
                    subU(:,lotid)=0;
                    MITIGATE(cat(1,Lottype{lotid,2}))=num2cell(rem(lotchoice{lotid,8},2));
                    openhouse=openhouse-1;
                end
            end
            
            % maxbids=zeros(length(Nlots(TSTART)),1);
            conlist=(1:length(CONINFO(:,1)))';
            ifilled=find(cat(1,lotchoice{:,4})==1);
            istillvac=find(cat(1,lotchoice{:,4})==0);
            popin=cat(1,lotchoice{ifilled,5});
            popout=conlist(~ismember(conlist,popin));
            lotchoice(istillvac,7)=num2cell(Paskhouse(istillvac)./(1+discount));
            lotchoice(istillvac,6)=num2cell(ones(length(istillvac),1)*TSTART+1);
            lotchoice(istillvac,8)=num2cell(zeros(length(istillvac),1));
            
            
            %%%%%%%% Utility check %%%%%%%%%%%
            % realulot=zeros(Nconsumers,3);
            % inhouselist=conlist(popin);
            % maxulot=zeros(1,[]);
            % imaxulot=zeros(1,[]);
            % for c=1:Nconsumers
            %     [maxulot(c),imaxulot(c)]=max(U(c,:),[],2);
            % end
            % for c=1:length(inhouselist)
            %     ireallot=find(cat(1,lotchoice{:,5})==inhouselist(c));
            %     realulot(inhouselist(c),1:3)=[ireallot lotchoice(ireallot,5) U(inhouselist(c),ireallot)];
            % end
            % maxuset=[imaxulot' lotchoice(imaxulot,5) maxulot'];
            % fullset=[conlist realulot maxuset];
            % utildiff(TSTART)=mean(fullset(fullset(:,4)~=0,7)-fullset(fullset(:,4)~=0,4));
            % pctutildiff(TSTART)=mean(fullset(fullset(:,4)~=0,4)./fullset(fullset(:,4)~=0,7));
            % incomediff=mean(Income(popin));
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if killflag==1
                break
            end
            
            %iterate until difference between observed and expected prices is below
            %threshold, consumer utility is maximized, and vacancies are minimized
            iter=iter+1;
            deltadiff(deltaiter)=meandiffprice;
            meandiffprice=mean(abs(Paskhouse-cat(1,lotchoice{:,7})));
            diffprice(1:length(Paskhouse),iter)=cat(1,lotchoice{:,7})-Paskhouse;
            diffvac(iter)=length(istillvac);
            
            deltaiter=deltaiter+1;
            
            for lt=1:HT
                if isempty(find(cat(1,lotchoice{:,3})==lt,1))==1
                    continue
                end
                meanprice(lt,iter)=mean(cat(1,lotchoice{cat(1,lotchoice{:,3})==lt,7}));
            end
            iterhmc(length(iterhmc)+1)=avghousemp(TSTART);
            deltadiff(deltaiter)=meandiffprice;
            deltafac=(deltadiff(deltaiter-1)-deltadiff(deltaiter))/deltadiff(deltaiter-1);
            diffcheck=initialdiff-meandiffprice;
            utildiffcheck=initialutildiff-utildiff(TSTART);
            if iter > 2
                vaccheck=diffvac(iter-1)-diffvac(iter);
            end
            
            if diffcheck < 0 && utildiffcheck < 0 && abs(deltafac) < 0.01
                Paskhouse=savePaskhouse;
                killflag=1;
                continue
            elseif diffcheck < 0 && utildiffcheck < 0 && vaccheck < 0
                Paskhouse=savePaskhouse;
                killflag=1;
                continue
            end
            
            savePaskhouse=Paskhouse;
            Paskhouse=Paskhouse+abs(meandiffprice/mean(meanprice(meanprice(:,iter)~=0,iter))).*diffprice(:,iter);
            
            RENT(cat(1,lotchoice{:,2}),1)=cat(1,lotchoice{:,7});
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%    Population    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        popinhouse=find(cat(1,CONINFO{:,8})==1);
        ioldpop=find(cat(1,CONINFO{:,8})==0);
        POP(1:TSTART)=round(length(CONINFO(:,1))./(1+POPGROW).^(TSTART:-1:1));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%    Broker Info    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %@@@@@@@ Broker Projection Models @@@@@@@@@@@@@@@@@@
        subplandinfo=cat(2,LANDINFO{3,TSTART});
        lotlocate=zeros([],2);
        brkrlocate=zeros([],2);
        for il=1:length(Lottype(:,1))
            %lotlocate=[lotid lot_index]
            lotind=cat(1,Lottype{il,2});
            lotlocate(length(lotlocate(:,1))+1:length(lotlocate(:,1))+length(lotind),:)=...
                [ones(length(lotind),1)*il lotind];
            %brkrlocate=[brokerid index]
            brokerassign=unique(HBROKER(Lottype{il,2}));
            brokerassign=brokerassign(brokerassign~=0);
            Lottype{il,10}=brokerassign;
            brkrlocate(length(brkrlocate(:,1))+1:length(brkrlocate(:,1))+length(brokerassign),:)=...
                [ones(length(brokerassign),1)*il brokerassign];
        end
        warning('off','all');
        for ibr=1:Nbrokers
            bhood=find(HBROKER==ibr);
            brokerind=unique(brkrlocate(cat(1,Lottype{:,10})==ibr,1));
            iblots=cell2mat(lotchoice(brokerind,[1 3 4 5]));
            
            if isempty(iblots) == 1
                continue
            else
                % Calculate bids per lot
                nbids=zeros([],1);
                for b=1:length(iblots(:,1))
                    bidover=(Phousebid(:,iblots(b,1)) > 0);
                    nbids(b,1)=length(find(bidover==1));
                end
                
                % Calculate average consumer utilities, preferences for rent
                % projection
                isnotvac=(iblots(:,3)==1);
                if isempty(find(isnotvac,1))==1
                    conintue
                else
                    brkravgstats(ibr,:)=[median(cat(1,CONINFO{iblots(isnotvac,4),1})) ...
                        median(cat(1,CONINFO{iblots(isnotvac,4),3})) ...
                        mean(cat(1,CONINFO{iblots(isnotvac,4),4})) ...
                        mean(cat(1,CONINFO{iblots(isnotvac,4),5})) ...
                        mean(cat(1,CONINFO{iblots(isnotvac,4),6}))];
                end
                AVGUTIL(brokerind(isnotvac))=CONINFO(iblots(isnotvac,4),9);
                BIDLEVEL(brokerind(isnotvac))=num2cell(cat(1,lotchoice{brokerind(isnotvac),7})./...
                    (cat(1,Lottype{brokerind(isnotvac),6})+...
                    discount*subplandinfo(cat(1,lotchoice{brokerind(isnotvac),2})).*...
                    cat(1,Lottype{brokerind(isnotvac),3})));
                sampleinfo=zeros(length(brokerind),6);
                sampleinfo(:,1)=cat(1,lotchoice{brokerind,7});
                sampleinfo(:,2)=cat(1,lotchoice{brokerind,3});
                sampleinfo(:,3)=cat(1,lotchoice{brokerind,4});
                sampleinfo(:,4)=nbids;
                sampleinfo(:,5)=cat(1,lotchoice{brokerind,7})./(cat(1,Lottype{brokerind,6})+...
                    discount*subplandinfo(cat(1,lotchoice{brokerind,2})).*cat(1,Lottype{brokerind,3}));
                sampleinfo(:,6)=cat(1,lotchoice{brokerind,8});
                subexpthouse=zeros(HT,1);
                for lt=1:HT
                    ils=(iblots(:,2)==lt);
                    houseinfo(lt,2,ibr,TSTART)=z(lt,1);
                    houseinfo(lt,3,ibr,TSTART)=z(lt,2);
                    if isempty(find(ils,1))==1
                        subexpthouse(lt)=0;
                        houseinfo(lt,[1 4:7],ibr,TSTART)=0;
                    else
                        subexpthouse(lt)=mean(cat(1,lotchoice{iblots(ils,1),7}));
                        houseinfo(lt,1,ibr,TSTART)=mean(sampleinfo(ils,1));
                        houseinfo(lt,4,ibr,TSTART)=mean(sampleinfo(ils,4));
                        houseinfo(lt,5,ibr,TSTART)=min(sampleinfo(ils,5));
                        houseinfo(lt,6,ibr,TSTART)=length(find(ils==1));
                        houseinfo(lt,7,ibr,TSTART)=mean(cat(1,coastprox{cat(1,Lottype{iblots(ils,1),2})}));
                    end
                end
                ilotlocate=ismember(lotlocate(:,1),iblots(:,1));
                % Rent expectations per parcel for each housing type
                EXPTHOUSE(lotlocate(ilotlocate,2),TSTART+1)=subexpthouse(...
                    cat(1,lotchoice{lotlocate(ilotlocate,1),3}));
            end
        end
        
        warning('on','all');
        for lt=1:HT
            ilt=(cat(1,lotchoice{:,3})==lt);
            % average bid level relative to base housing cost
            bidlevel(lt,1:TSTART+1)=mean(cat(1,lotchoice{ilt,7})./(cat(1,Lottype{ilt,6})+...
                discount*subplandinfo(cat(1,lotchoice{ilt,2})).*...
                cat(1,Lottype{ilt,3})));
        end
        numlt(:,1:TSTART)=repmat(histc(cat(1,Lottype{:,5}),1:HT),1,TSTART);
        
        bcheck(:,1:Nbrokers)=houseinfo(:,1,1:Nbrokers,TSTART);
        ibcheck=(bcheck' ~= 0);
        htexist=ismember(1:HT,cat(1,Lottype{:,5}));
        hset=htset(htexist);
        for lt=1:HT
            ihtexist=(ibcheck(:,lt)==1);
            if isempty(find(ihtexist,1))==1
                continue
            else
                avgbrokervar(lt,TSTART)=mean(var(brokerbestabsSAVE(ihtexist,lt,1:TSTART),0,3));
                probloss(lt,TSTART)=mean(sum((brokerbestabsSAVE(ihtexist,lt,1:TSTART)>0),3)./...
                    length(1:TSTART));
                abserror=brokerbestabsSAVE(ihtexist,lt,1:TSTART);
                abserror=reshape(abserror,length(abserror(:,1,1))*TSTART,1);
                [mu,sigma]=normfit(abserror);
                phat(lt,:)=[mu sigma];
                probeven(lt,TSTART)=cdf('norm',0,phat(lt,1),phat(lt,2));
                probover(lt,TSTART)=length(find(abserror > 0))/length(abserror);
                probunder(lt,TSTART)=length(find(abserror < 0))/length(abserror);
                overvalue(lt,TSTART)=icdf('norm',max(min(probeven(lt,TSTART)+...
                    (1-probeven(lt,TSTART))*probover(lt,TSTART),0.99),0.01),phat(lt,1),phat(lt,2));
                undervalue(lt,TSTART)=icdf('norm',max(min(probeven(lt,TSTART)*...
                    (1-probunder(lt,TSTART)),0.99),0.01),phat(lt,1),phat(lt,2));
                maxvalue(lt,TSTART)=icdf('norm',probeven(lt,TSTART)+...
                    (1-probeven(lt,TSTART))*0.99,phat(lt,1),phat(lt,2));
                minvalue(lt,TSTART)=icdf('norm',probeven(lt,TSTART)*...
                    (1-0.99),phat(lt,1),phat(lt,2));
            end
            
        end
        ihtnexist=(htexist==0);
        isimvar=hset(ismember(htset(htexist),min(simlotrange(htset(ihtnexist),1)):...
            max(simlotrange(htset(ihtnexist),2))));
        avgbrokervar(ihtnexist,TSTART)=max(avgbrokervar(isimvar,TSTART));
        probloss(ihtnexist,TSTART)=alpha_gain/(alpha_gain+alpha_loss);
        overvalue(ihtnexist,TSTART)=mean(maxvalue(isimvar,TSTART));
        undervalue(ihtnexist,TSTART)=mean(minvalue(isimvar,TSTART));
        probover(ihtnexist,TSTART)=alpha_gain/(alpha_gain+alpha_loss);
        probunder(ihtnexist,TSTART)=alpha_loss/(alpha_gain+alpha_loss);
        RENTGRAD(lotlocate(:,2))=RENT(lotlocate(:,2),TSTART).*1./cat(1,Lottype{lotlocate(:,1),3});
        
        subfarminfo=cat(2,LANDINFO{1,TSTART});
        subfarminfo(iurblist)=0;
        iNfarmers=unique(subfarminfo);
        iNfarmers=iNfarmers(iNfarmers~=0);
        LANDINFO(1,TSTART+1)=mat2cell(subfarminfo,NCELLS,1);
        landproj=mean(Plandproj(:,TSTART))*(1+min(max(randn(size(landproj)),-0.5),0.5));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%    RESULTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for lt=1:HT
            ilt=find(cat(1,lotchoice{:,3})==lt);
            avgrent(lt,1:TSTART+1)=mean(cat(1,lotchoice{ilt,7}));
            iltocc=ismember(ilt,ifilled);
            iltvac=ismember(ilt,istillvac);
            budget_lt(lt,1:TSTART+1)=sum(cat(1,lotchoice{ilt(iltocc),7})-...
                subplandinfo(cat(1,lotchoice{ilt(iltocc),2})).*...
                cat(1,Lottype{ilt(iltocc),3})*discount-cat(1,Lottype{ilt(iltocc),6}));
            vac_ccost(lt,1:TSTART+1)=sum(discount*cat(1,Lottype{ilt(iltvac),6})+...
                subplandinfo(cat(1,lotchoice{ilt(iltvac),2})).*cat(1,Lottype{ilt(iltvac),3}));
        end
        carrycost(TSTART)=sum(vac_ccost(:,TSTART));
        BUDGET(TSTART)=sum(budget_lt(:,TSTART));
        
        oldincome(1:TSTART)=mean(cat(1,CONINFO{popout,1}));
        devcells(1,TSTART)=length(iurblist);
        devcells(2,TSTART)=devcells(1,TSTART)/NCELLS;
        
        vacantlots(TSTART)=length(istillvac);
        leftoverpop(TSTART)=length(popout);
        
        vacrate(TSTART)=vacantlots(TSTART)/Nlots(TSTART);
        nohouserate(TSTART)=leftoverpop(TSTART)/POP(TSTART);
        
        agrland(TSTART)=length(find(BASELAYER == 0 & reshape(SCAPE,NCELLS,1) == 1));
        
        consumerstats(1,TSTART)=length(CONINFO(:,1));
        consumerstats(4,TSTART)=mean(housemp);
        consumerstats(2,TSTART)=mean(cat(1,lotchoice{ifilled,7}));
        consumerstats(3,TSTART)=mean(cat(1,AVGUTIL{ifilled}));
        
        ifill=ismember(lotlocate(:,1),ifilled);
        BUILDTIME(lotlocate(:,2))=cat(1,Lottype{lotlocate(:,1),9});
        BIDLEVELMAP(lotlocate(ifill,2),1:TSTART)=repmat(cat(1,BIDLEVEL{lotlocate(ifill,1)}),1,TSTART);
        VACLAND(cat(1,vacland{1,1}),1:TSTART)=1;
        AVGRENT(lotlocate(:,2),1:TSTART)=repmat(cat(1,lotchoice{lotlocate(:,1),7}),1,TSTART);
        LOTTYPE(lotlocate(:,2),1:TSTART)=repmat(cat(1,Lottype{lotlocate(:,1),5}),1,TSTART);
        BASEMAP(lotlocate(:,2))=BASELAYER(lotlocate(:,1));
        INCOME(lotlocate(ifill,2),1:TSTART)=repmat(cat(1,CONINFO{cat(1,lotchoice{lotlocate(ifill,1),5}),1}),1,TSTART);
        
        Rpop(1:TSTART)=length(ifilled);
        Rvacrate(1:TSTART)=vacrate(TSTART);
        Rvaclots(1:TSTART)=vacantlots(TSTART);
        Rleftoverpop(1:TSTART)=leftoverpop(TSTART);
        setupmap=Sfarmmap.AGLAYER;
        Ufinset=zeros(length(ifilled),1);
        
        for ires=1:length(ifilled)
            c=lotchoice{ifilled(ires),5};
            hopt=((CONINFO{c,1}-cat(1,travelcost{cat(1,lotchoice{ifilled(ires),2})})-...
                avgrent(:,TSTART)).^CONINFO{c,3}).*(cat(1,Lottype{ifilled(ires),4}).^...
                CONINFO{c,4}).*(cat(1,Lottype{ifilled(ires),3}).^CONINFO{c,5}).*...
                (cat(1,Lottype{ifilled(ires),7}).^CONINFO{c,6});
            
            profopt=(avgrent(:,TSTART)-ones(HT,1)*subplandinfo(cat(1,lotchoice{ifilled(ires),2}))-...
                cat(1,Lottype{ifilled(ires),6}))./cat(1,Lottype{ifilled(ires),3});
            [imaxp,jmaxp]=max(profopt,[],1);
            profset(jmaxp,1:TSTART)=profset(jmaxp,1)+1;
            
            [imaxu,jmaxu]=max(hopt,[],1);
            idealset(jmaxu,1:TSTART)=idealset(jmaxu,1)+1;
        end
        
        % figure(1)
        % surf(reshape(LOTTYPE(:,TSTART),NLENGTH,NWIDTH));
        % axis ij;
        % view(0,90);
        % title('Lot Types, t=TSTART')
        % set(gca,'clim',[1 HT])
        % colorbar
        % MLT(TSTART)=getframe(gcf);
        
        iconleave=0;
        Nconsumers=length(cat(1,CONINFO{:,1}));
        %%
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %@@@@@@@@@@@@@@@@@@@@@@@@@    DYNAMICS    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        for t=TSTART+1:TMAX
            
%             t
            ccost=ccost_base;
            POP(t)=ceil(length(CONINFO(:,1))*(1+POPGROW));
            
            %Existing houses back on market <><><><><><><><><><><><><>
            ileave=find(cat(1,lotchoice{:,6})==t);
            if isempty(find(ileave,1))==0
                returncon=(cat(1,lotchoice{ileave,4})==1);
                CONINFO(cat(1,lotchoice{ileave(returncon),5}),2)=num2cell(t+ceil(searchtimemin+...
                    (searchtimemax-searchtimemin)*rand(length(ileave(returncon)),1)));
                lotchoice(ileave,4)=num2cell(zeros(length(ileave),1));
                iwasocc=(cat(1,lotchoice{ileave,5})~=0);
                CONINFO(cat(1,lotchoice{ileave(iwasocc),5}),8)=...
                    num2cell(zeros(length(cat(1,lotchoice{ileave(iwasocc),5})),1));
                openlothist=histc(cat(1,lotchoice{ileave,3}),1:HT);
                newopenlots(:,t)=max(reshape([openlothist zeros(size(openlothist))],HT,2),[],2);
                iexisthouse=ileave;
            end
            numlt(:,t)=histc(cat(1,Lottype{:,5}),1:HT);
            
            subplandinfo=cat(2,LANDINFO{3,t-1});
            for in=1:length(iNfarmers)
                ifarmland=find(LANDINFO{1,t}==nf);
                subplandinfo(ifarmland)=LOCWGHT*wtaland(iNfarmers(in),t-1)+...
                    REGWGHT*mean(mean(subplandinfo));
                LANDINFO(3,t)=mat2cell(subplandinfo,NCELLS,1);
%                 Dynplandmap(:,:,t)=reshape(cat(2,LANDINFO{3,t}),NLENGTH,NWIDTH);
                Dynplandmap(:,t)=cat(2,LANDINFO{3,t});
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%    Storm Simulator    %%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Random Storm Generation
            %     stormfreq=0.1;
            stormint=round(1/sum(Psevere));
%             stormgen=TSTART:stormint:TMAX;
            stormgen=TMAX+1;
            if isempty(find(stormgen==t,1))==0
                probsurf=cat(1,Pflood{:});
                impactprob=max(cat(1,Pflood{:}))-rand(1);
                iimpact=find(probsurf > impactprob);
                inoimpact=find(probsurf <= impactprob);
                lothit=ismember(lotlocate(:,2),iimpact);
                ihitlist=unique(lotlocate(lothit,1));
                lotmiss=ismember(lotlocate(:,2),inoimpact);
                IMPACT(iimpact,t)=num2cell(ones(length(iimpact),1));
                IMPACT(inoimpact,t)=num2cell(ones(length(inoimpact),1));
                DAMAGE(cat(1,lotchoice{ihitlist,2}))=num2cell(cat(1,DAMAGE{cat(1,lotchoice{ihitlist,2})})+...
                    cat(1,Cdam{ihitlist,t}).*cat(1,lotchoice{ihitlist,7}));
                TSI(iimpact)=num2cell(zeros(length(iimpact),1));
                TSI(inoimpact)=num2cell(cat(1,TSI{inoimpact})+1);
            end
            % %%% Random Storm Generation
            %     stormgen=100*rand(1);
            %     if stormgen <= stormthresh
            %         probsurf=cat(1,Pflood{:});
            %         impactprob=1-rand(1);
            %         iimpact=find(probsurf > impactprob);
            %         inoimpact=find(probsurf <= impactprob);
            %         lothit=ismember(lotlocate(:,2),iimpact);
            %         ihitlist=unique(lotlocate(lothit,1));
            %         lotmiss=ismember(lotlocate(:,2),inoimpact);
            %         IMPACT(iimpact,t)=num2cell(ones(length(iimpact),1));
            %         IMPACT(inoimpact,t)=num2cell(ones(length(inoimpact),1));
            %         DAMAGE(cat(1,lotchoice{ihitlist,2}))=num2cell(cat(1,DAMAGE{cat(1,lotchoice{ihitlist,2})})+...
            %             Cdam.*cat(1,lotchoice{ihitlist,7}));
            %         TSI(iimpact)=num2cell(zeros(length(iimpact),1));
            %         TSI(inoimpact)=num2cell(cat(1,TSI{inoimpact})+1);
            %     end
            %     lotdamage(lotlocate(lothit,1))=cat(1,DAMAGE{cat(1,lotchoice{ihitlist,2})})
            tsilot=zeros(length(Lottype(:,1)),1);
            for nl=1:length(Lottype(:,1))
                tsilot(nl)=min(cat(1,TSI{Lottype{nl,2}}));
            end
            damcoef(1:length(CONINFO(:,1)),t)=num2cell(ones(length(CONINFO(:,1)),1));
            %     for in=1:length(CONINFO(:,1))
            %         damcoef{in,t}=(1-CONINFO{in,7})*damcoef{in,t-1}.*(1./tsilot)+...
            %             CONINFO{in,7}*damcoef{in,t-1};
            %     end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%    Developer's Decisions    %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Demand, by each lot
            reldemand=zeros(Nlots(t),1);
            ltype=zeros(Nlots(t),1);
            extrabid=zeros(Nlots(t),1);
            
            for lt=1:HT
                iocc=(houseinfo(lt,1,:,t-1)~=0);
                if isempty(find(iocc,1))==1
                    bidtot(lt,t)=0;
                else
                    bidtot(lt,t)=mean(houseinfo(lt,4,iocc,t-1),3);
                end
            end
            for iii=1:Nlots(t)
                bidoverx=(Phousebid(:,iii) >= (Paskhouse(iii)));
                extrabid(iii,1)=length(find(bidoverx ==1));
                reldemand(iii,1)=extrabid(iii,1)-newopenlots(lotchoice{iii,3},t);
            end
            % Rank housing types by their proportion of bids received last period
            trackbids=[extrabid reldemand cat(1,Lottype{:,5})];
            demandrank=sortrows(trackbids,-1); %[#of_bids(sorted) relative_demand lotid ind lotsize lottype]
            
            for lt=1:HT
                pctbuildold(lt,t)=sum(demandrank(demandrank(:,3)==lt,1))/...
                    sum(demandrank(:,1));
                pctbuildnew(lt,t)=sum(demandrank(demandrank(:,3)==lt,1))/...
                    sum(demandrank(:,1));
            end
            
            %%% Developer's Population Predictions %%%
            
            % mirror models
            nproj(indclass1) = POP(t-1) - POP(t-2) + (1-dd(indclass1))*(0.5*POP(t-1) - (POP(t-1)-POP(t-2)));
            
            % mean mode
            for j = 1:length(indclass2)
                nproj(indclass2(j)) = mean(POP((t-1):-1:(t-dd(indclass2(j)))));
            end
            
            %cycle model
            nproj(indclass3) = POP(t-max(1,dd(indclass3)));
            
            % projection model
            for j = 1:length(indclass4)
                yyy(j,1:dd(indclass4(j))) = POP((t-dd(indclass4(j))):(t-1));
            end
            sumy = sum(yyy,2);
            sumxy = sum(xxx.*yyy,2);
            slopes = (dd(indclass4)'.*sumxy - sumx.*sumy)./(dd(indclass4)'.*sumx2 - sumx.*sumx);
            intercepts = (sumy - slopes.*sumx)./dd(indclass4)';
            nproj(indclass4) = slopes + intercepts;
            
            % rescale model
            nproj(indclass5) = dd(indclass5)*POP(t-1);
            
            errorsq = (1-DELTA)*errorsq + DELTA*(POP(t)-nproj).^2;
            [best,ibest] = min(errorsq,[],2);
            bestPOPSAVE(1:Ndevelopers,t) = best;
            ibestPOPSAVE(1:Ndevelopers,t) = ibest;
            nprojSAVE(1:Ndevelopers,t) = nproj(ibest);
            
            numnewhouses(t)=round((nproj(ibest)-(length(ifilled)+length(ioldpop)-length(iconleave))));
            numoldhouses(t)=length(ioldpop)-length(istillvac);
            houseset=max(ceil(nproj(ibest)-sum(numlt(:,t))-length(istillvac)),0);
            newhouseset(:,t)=round(houseset.*(bidtot(:,t)./sum(bidtot(:,t))));
            
            %Spatial Rent discounting for all vacant cells
            iurblist=find(BASELAYER == 1);
            ivac=(VACLAND(:,t-1) ~= 0);
            ivaclist=find(ivac==1);
            
            for lt=1:HT
                % identify similar lots on which to base rent projections for
                % housing types that are not in range or do not yet exist
                isimlots=(ismember(cat(1,Lottype{:,5}),simlotrange(lt,1):...
                    simlotrange(lt,2))==1 & cat(1,lotchoice{:,4})~=0);
                isimcells=unique(cat(1,Lottype{ismember(cat(1,Lottype{:,1}),cat(1,lotchoice{isimlots,1})),2}));
                
                %Regional Stats
                if isempty(find(cat(1,lotchoice{:,3})==lt,1))==1
                    simlots_income(lt)=median(cat(1,CONINFO{cat(1,lotchoice{isimlots,5}),1}));
                    %             simlots_util(lt)=median(cat(1,AVGUTIL{isimlots}));
                    simlots_alpha(lt)=mean(cat(1,CONINFO{cat(1,lotchoice{isimlots,5}),3}));
                    simlots_beta(lt)=mean(cat(1,CONINFO{cat(1,lotchoice{isimlots,5}),4}));
                    simlots_gamma(lt)=mean(cat(1,CONINFO{cat(1,lotchoice{isimlots,5}),5}));
                    simlots_ampref(lt)=mean(cat(1,CONINFO{cat(1,lotchoice{isimlots,5}),6}));
                else
                    simlots_income(lt)=median(cat(1,CONINFO{cat(1,lotchoice{ifilled(...
                        cat(1,lotchoice{ifilled,3})==lt),5}),1}));
                    %         simlots_util(lt)=median(cat(1,AVGUTIL{cat(1,lotchoice{ifilled(...
                    %             cat(1,lotchoice{ifilled,3})==lt),5})}));
                    simlots_alpha(lt)=median(cat(1,CONINFO{cat(1,lotchoice{ifilled(...
                        cat(1,lotchoice{ifilled,3})==lt),5}),3}));
                    simlots_beta(lt)=median(cat(1,CONINFO{cat(1,lotchoice{ifilled(...
                        cat(1,lotchoice{ifilled,3})==lt),5}),4}));
                    simlots_gamma(lt)=median(cat(1,CONINFO{cat(1,lotchoice{ifilled(...
                        cat(1,lotchoice{ifilled,3})==lt),5}),5}));
                    simlots_ampref(lt)=median(cat(1,CONINFO{cat(1,lotchoice{ifilled(...
                        cat(1,lotchoice{ifilled,3})==lt),5}),6}));
                end
                
                ireglot=(cat(1,Lottype{:,5})==lt);
                if isempty(find(ireglot,1))==1
                    continue
                else
                    regionaldist(lt,t)=mean(Sdist2cbd.dist2cbd(cat(1,Lottype{ireglot,2})));
                    regionalrent(lt,t)=mean(EXPTHOUSE(cat(1,Lottype{ireglot,2}),t));
                end
            end
            reg_util=median(cat(1,AVGUTIL{ifilled}));
            reg_income=median(cat(1,CONINFO{cat(1,lotchoice{ifilled,5}),1}));
            reg_alpha=mean(mean(cat(1,CONINFO{cat(1,lotchoice{ifilled,5}),3})));
            reg_beta=mean(cat(1,CONINFO{cat(1,lotchoice{ifilled,5}),4}));
            reg_gamma=mean(cat(1,CONINFO{cat(1,lotchoice{isimlots,5}),5}));
            reg_ampref=mean(cat(1,CONINFO{cat(1,lotchoice{isimlots,5}),6}));
            
            %Hedonic housing price estimation
            %     B=mvregress([ones(length(ifilled),1) cat(1,Lottype{ifilled,3}) ...
            %         cat(1,Lottype{ifilled,8}) cat(1,Lottype{ifilled,7})],...
            %         cat(1,lotchoice{ifilled,7}));
            %     rentmdl=fitlm([cat(1,Lottype{ifilled,3}) cat(1,Lottype{ifilled,4})...
            %         cat(1,Lottype{ifilled,8}) cat(1,Lottype{ifilled,7}) ...
            %         cat(1,CONINFO{cat(1,lotchoice{ifilled,5}),1})],...
            %         cat(1,lotchoice{ifilled,7}));
            rentmdl=fitlm([cat(1,Lottype{ifilled,3}) cat(1,Lottype{ifilled,8})...
                cat(1,Lottype{ifilled,7}) cat(1,CONINFO{cat(1,lotchoice{ifilled,5}),1})],...
                cat(1,lotchoice{ifilled,7}));
            
            ddist2hznnei=zeros(NLENGTH,NWIDTH);    %distance to horizontal neighbor from icenter
            ddist2vrtnei=zeros(NLENGTH,NWIDTH);
            potentialbuy=find(BASELAYER==0 & reshape(SCAPE,NCELLS,1) == 1);
            for tl=1:length(potentialbuy)
                %%% Distance calc needs to be revised when applied to irregular grid
                [vacrow,vaccol]=ind2sub([NLENGTH NWIDTH],potentialbuy(tl));
                % pre-loaded distance matrix
                inddist2dev=Sdist.distmat{potentialbuy(tl)};
                
                subdist2dev=[inddist2dev(iurblist) iurblist];
                sortdevdist=sortrows(subdist2dev,1);
                iclosedev=sortdevdist(:,2);
                for lt=1:HT
                    zonecheck=cat(2,ZONES{ceil(find(cat(1,ZONES{:,1})==potentialbuy(tl))/...
                        (NCELLS/length(ZONES(:,1)))),2:3});
                    if z(lt,1) >= zonecheck(1) && z(lt,1) <= zonecheck(2)
                        warning('off','all');
                        likelotcount=round(numlt(lt,t));
                        icountcells=iclosedev(1:round(length(iurblist)*PCTSEARCH));
                        ilocatelot=find(ismember(lotlocate(:,2),icountcells)==1);
                        iutilcells=lotlocate(ilocatelot(cat(1,AVGUTIL{lotlocate(ilocatelot,1)})~=0 & ...
                            cat(1,Lottype{lotlocate(ilocatelot,1),3})>=simlotrange(lt,1) & ...
                            cat(1,Lottype{lotlocate(ilocatelot,1),3})<=simlotrange(lt,2)),2);
                        distutils=inddist2dev(iutilcells);
                        if isempty(distutils)==1
                            termrunflag=1;
                        end
                        warning('on','all');
                        if likelotcount == 0
                            %                     subRENTPROJ(lt)=B(1)+B(2)*z(lt,1)+B(3)*...
                            %                         travelcost{potentialbuy(tl)}+B(4)*coastprox{potentialbuy(tl)};
                            %                     subRENTPROJ(lt)=predict(rentmdl,[z(lt,1) z(lt,2) ...
                            %                         travelcost{potentialbuy(tl)} coastprox{potentialbuy(tl)} ...
                            %                         simlots_income(lt)]);
                            subRENTPROJ(lt)=predict(rentmdl,[z(lt,1) ...
                                travelcost{potentialbuy(tl)} coastprox{potentialbuy(tl)} ...
                                simlots_income(lt)]);
                        else
                            % find close lots upon which to base rent projections
                            nclosecells=max(likelotcount*min(z(lt,1),1),round(length(iurblist)*PCTSEARCH));
                            iclosecells=find(ismember(cat(1,Lottype{:,2}),iclosedev(1:nclosecells))==1);
                            icloselot=find(cat(1,Lottype{lotlocate(iclosecells,1),5}) == lt,nclosecells,'first');
                            devrents=EXPTHOUSE(lotlocate(iclosecells(icloselot),2),t);
                            distrents=inddist2dev(lotlocate(iclosecells(icloselot),2));
                            [closecelldist,iclstcell]=min(distrents,[],1);
                            iclosestcell=lotlocate(iclosecells(icloselot(iclstcell)),2);
                            
                            if isempty(icloselot)==1
                                % If similar housing type is not found in search
                                % area ...
                                %                         locrentproj=B(1)+B(2)*z(lt,1)+B(3)*...
                                %                             travelcost{potentialbuy(tl)}+B(4)*coastprox{potentialbuy(tl)};
                                %                         locrentproj=predict(rentmdl,[z(lt,1) z(lt,2) ...
                                %                             travelcost{potentialbuy(tl)} coastprox{potentialbuy(tl)} ...
                                %                             simlots_income(lt)]);
                                locrentproj=predict(rentmdl,[z(lt,1) ...
                                    travelcost{potentialbuy(tl)} coastprox{potentialbuy(tl)} ...
                                    simlots_income(lt)]);
                                regrentproj=regionalrent(lt,t)-margtc*...
                                    (Sdist2cbd.dist2cbd(vacrow,vaccol)-regionaldist(lt,t));
                                subRENTPROJ(lt)=LOCWGHT*locrentproj+...
                                    REGWGHT*regrentproj;
                            else
                                locrentproj=(sum(-distrents.*devrents)/sum(-distrents))-...
                                    margtc*distrents(iclstcell);
                                %                         regrentproj=predict(rentmdl,[z(lt,1) z(lt,2) ...
                                %                             travelcost{potentialbuy(tl)} coastprox{potentialbuy(tl)} ...
                                %                             simlots_income(lt)]);
                                regrentproj=predict(rentmdl,[z(lt,1) ...
                                    travelcost{potentialbuy(tl)} coastprox{potentialbuy(tl)} ...
                                    simlots_income(lt)]);
                                %                         regrentproj=regionalrent(lt,t)-margtc*...
                                %                             (Sdist2cbd.dist2cbd(vacrow,vaccol)-regionaldist(lt,t));
                                subRENTPROJ(lt)=LOCWGHT*locrentproj+...
                                    REGWGHT*regrentproj;
                            end
                        end
                        subRETURN(lt)=(1-discount)*(subRENTPROJ(lt)-ccost(lt))/...
                            z(lt,1)-carrycost(t-1)/ceil(sum(newhouseset(:,t).*z(:,1)));
                        
                        % Esitmate potential high/low return for risk aversion calc
                        highRETURN(lt)=(1-discount)*((overvalue(lt,t-1)+...
                            subRENTPROJ(lt))-ccost(lt))/z(lt,1)-carrycost(t-1)/...
                            ceil(sum(newhouseset(:,t).*z(:,1)));
                        lowRETURN(lt)=(1-discount)*((subRENTPROJ(lt)+...
                            undervalue(lt,t-1))-ccost(lt))/z(lt,1)-carrycost(t-1)/...
                            ceil(sum(newhouseset(:,t).*z(:,1)));
                        
                        % Potential gain /loss relative to zero
                        if highRETURN(lt) > 0
                            potgain(lt)=(highRETURN(lt)/alpha_gain)^(1/alpha_gain);
                        elseif highRETURN(lt) <= 0
                            potgain(lt)=0;
                            potloss(lt)=(abs(highRETURN(lt))/alpha_loss)^(1/alpha_loss);
                        end
                        if lowRETURN(lt) < 0
                            potloss(lt)=potloss(lt)+...
                                (abs(lowRETURN(lt))/alpha_loss)^(1/alpha_loss);
                        elseif lowRETURN(lt) >= 0
                            potloss(lt)=0;
                            potgain(lt)=potgain(lt)+(abs(lowRETURN(lt))/alpha_gain)^(1/alpha_gain);
                        end
                        subEU_dev(lt)=(potgain(lt))/(potgain(lt)+potloss(lt));
                    else
                        subRETURN(lt)=0;
                    end
                end
                RENTPROJ(potentialbuy(tl),t)=mat2cell(subRENTPROJ,HT,1);
                RETURN(potentialbuy(tl),t)=mat2cell(subRETURN,HT,1);
                EU_dev(potentialbuy(tl),t)=mat2cell(subEU_dev,HT,1);
                maxcount(potentialbuy(tl))=length(find(EU_dev{potentialbuy(tl)} == ...
                    max(EU_dev{potentialbuy(tl)})));
            end
            Exptret(:,t)=mean(cat(2,RETURN{potentialbuy,t}),2);
            if isempty(find(newhouseset(:,t),1))==1
                LANDBUDGET(t)=0;
            else
                LANDBUDGET(t)=BUDGET(t-1)/round(sum(newhouseset(:,t).*z(:,1)));
            end
            
            farmprojcells=zeros(length(potentialbuy),3);     %[profitability ltype ind]
            [maxret,retind]=max(cat(2,RETURN{:,t}),[],1);
            [EU_value,EU_ind]=sort(cat(2,EU_dev{:,t}),1,'descend');
            [maxeuval,maxeuind]=max(cat(2,EU_dev{:,t}),[],1);
            MAXRETURN(potentialbuy)=maxret;
            RETIND(potentialbuy)=retind;
            EUIND(:,potentialbuy)=EU_ind;
            MAXEUIND(potentialbuy)=maxeuind;
            MAXEUIND(maxcount > 1)=RETIND(maxcount > 1);
            for tl=1:length(potentialbuy)
                subRETURN=RETURN{potentialbuy(tl),t};
                EUrankret(potentialbuy(tl),:)=subRETURN(EUIND(:,potentialbuy(tl)));
                iposeu=find(EUrankret(potentialbuy(tl),:) > 0);
                poseuset=reshape(EUIND(iposeu,potentialbuy(tl)),length(iposeu),1);
                eubuildset=max(newhouseset(poseuset,t),1);
                subEU_dev=EU_dev{potentialbuy(tl),t};
                MAXRET(potentialbuy(tl),t)=sum(subEU_dev(poseuset).*...
                    subRETURN(EUIND(iposeu,potentialbuy(tl))))/...
                    sum(subEU_dev(poseuset));
                WTPMAP(potentialbuy(tl),t)=MAXRET(potentialbuy(tl),t)/discount;
            end
            farmprojcells(1:length(potentialbuy),:)=[maxret' retind' potentialbuy];
            RETURNPROJ(potentialbuy)=maxret;
            
            Dynmaxretmap(potentialbuy,t)=maxret;
            Dynretltmap(potentialbuy,t)=retind;
            Dyneultmap(potentialbuy,t)=maxeuind;
            
            maxretcells=zeros([],3);
            iundev=find(BASELAYER == 0);
            
            maxretcells(length(maxretcells(:,1))+1:length(maxretcells(:,1))+length(iundev),:)=...
                [RETURNPROJ(iundev) ones(length(iundev),1).*RETIND(iundev) ...
                iundev];
            rankmpcells=sortrows(maxretcells,-1);
            numnewacres(t)=ceil(sum(newhouseset(:,t).*z(:,1)));
            
            subrentproj=cat(2,RENTPROJ{:,t});
            maxRENTPROJ(potentialbuy)=subrentproj(retind);
            CCOST=ccost(retind);
            ZZ=z(retind,1);
            
            %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            %<><><><><><><><><><> Go to Land Market <><><><><><><><><><><><><><><><><><>
            %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            
            landdemand(1,t)=max(sum(numnewacres(t))-max(length(find(ivac==1)),...
                min(histc(cat(2,LANDINFO{1,t}),iNfarmers))),0);
            for rf=1:length(iNfarmers)
                ifarmcalc=find(LANDINFO{1,t}==iNfarmers(rf));
                subWTPMAP=WTPMAP(:,t);
                MAXEUMAP=MAXRET(:,t);
                wtpland(iNfarmers(rf),t)=sum(subWTPMAP(ifarmcalc))/length(ifarmcalc);
            end
            
            %%% LandMarket_coast_base
            %%%%%%%%%%%%%% ABM Land Market %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fdepart=zeros([],1);
            
            for tt=t
                if isempty(iNfarmers)==1 || landdemand(tt) == 0
                    ifarmtrans=[];
                    %         Plandproj(iNfarmers,tt)=mean([wtaland(iNfarmers,tt) wtpland(iNfarmers,tt)],2);
                    continue
                else
                    ipool=find(wtaland(iNfarmers,tt)<=wtpland(iNfarmers,tt));
                    ioutpool=find(wtaland(iNfarmers,tt)>=wtpland(iNfarmers,tt));
                end
                
                if isempty(ipool)==1
                    ifarmtrans=[];
                    %         Plandproj(iNfarmers,tt)=mean([wtaland(iNfarmers,tt) wtpland(iNfarmers,tt)],2);
                    Plandproj(iNfarmers,tt)=0.75*wtaland(iNfarmers,tt)+0.25*wtpland(iNfarmers,tt);
                    Paskland(iNfarmers,tt)=wtaland(iNfarmers,tt);
                    Pdevbid(iNfarmers,tt)=wtpland(iNfarmers,tt);
                    continue
                end
                subpoolinfo=cat(1,Farminfo{iNfarmers(ipool),2});
                epsilon(tt)=zeta*(landdemand(tt)-sum(subpoolinfo(:,1)))/...
                    (landdemand(tt)+sum(subpoolinfo(:,1)));
                %     epsilon(tt)=zeta*(landdemand(tt)-sum(Farmstats(iNfarmers(ipool),2,tt)))/...
                %         (landdemand(tt)+sum(Farmstats(iNfarmers(ipool),2,tt)));
                
                %     Paskland(iNfarmers(ipool),tt)=max(wtaland(iNfarmers(ipool),tt).*...
                %         (1+epsilon(tt)),wtaland(iNfarmers(ipool),tt));
                Paskland(iNfarmers(ipool),tt)=max(wtaland(iNfarmers(ipool),tt).*...
                    (1+epsilon(tt)),farmretinfo(iNfarmers(ipool)));
                Paskland(iNfarmers(ioutpool),tt)=wtaland(iNfarmers(ioutpool),tt);
                Plandproj(iNfarmers(ioutpool),tt)=wtpland(iNfarmers(ioutpool),tt);
                Pdevbid(iNfarmers(ipool),tt)=min(wtpland(iNfarmers(ipool),tt)*...
                    (1+epsilon(tt)),wtpland(iNfarmers(ipool),tt));
                
                
                %%% Farmer Selection Criteria %%%
                isell=find(Paskland(iNfarmers(ipool),tt)<= Pdevbid(iNfarmers(ipool),tt));
                if isempty(isell)==1
                    ifarmtrans=[];
                    Plandproj(iNfarmers(ipool),tt)=Pdevbid(iNfarmers(ipool),tt);
                    %         Plandproj(iNfarmers(ipool),tt)=0.75*Paskland(iNfarmers...
                    %             (ipool),tt)+0.25*Pdevbid(iNfarmers(ipool),tt);
                    continue
                else
                    transprice(iNfarmers(ipool(isell)),tt)=mean([Paskland(iNfarmers...
                        (ipool(isell)),tt) Pdevbid(iNfarmers(ipool(isell)),tt)],2);
                    %         transprice(iNfarmers(ipool(isell)),tt)=Paskland(iNfarmers(ipool(isell)),tt);
                    subfarmerinfo=zeros([],1);
                    for nnf=1:length(isell)
                        ifarmcalc=find(LANDINFO{1,t}==iNfarmers(ipool(isell(nnf))));
                        subfarmerinfo(nnf,1)=length(find(LANDINFO{1,t}==iNfarmers(ipool(isell(nnf)))));
                        %             ifarmcalc=find(AGLAYER == iNfarmers(ipool(isell(nnf))));
                        maxrent(iNfarmers(ipool(isell(nnf))),tt)=sum(maxRENTPROJ(ifarmcalc));
                        maxreturn(iNfarmers(ipool(isell(nnf))),tt)=sum(MAXRET(ifarmcalc,t))-...
                            discount*transprice(iNfarmers(ipool(isell(nnf))),tt).*...
                            subfarmerinfo(nnf,1);
                        %             maxreturn(iNfarmers(ipool(isell(nnf))),tt)=sum(maxret(ifarmcalc))-...
                        %                 discount*transprice(iNfarmers(ipool(isell(nnf))),tt).*...
                        %                 Farmstats(iNfarmers(ipool(isell(nnf))),2,1);
                    end
                    sellinfo=[maxreturn(iNfarmers(ipool(isell)),tt) ipool(isell) ...
                        subfarmerinfo];
                    %         sellinfo=[maxreturn(iNfarmers(ipool(isell)),tt) ipool(isell) ...
                    %             Farmstats(iNfarmers(ipool(isell)),2,1)];
                    lowprice=sortrows(sellinfo,1);
                    %         lowprice=sortrows(sellinfo,1);
                    landsupply=cumsum(lowprice(:,3));
                    ilandbuy=find(landsupply >= landdemand(tt),1,'first');
                    %         ifarmer=lowprice(:,2);
                    %         keyboard
                    if isempty(ilandbuy)==0
                        ibuy=1:ilandbuy;
                    elseif isempty(ilandbuy)==1
                        ibuy=1:length(isell);
                    end
                end
                
                fdepart(length(fdepart)+1:length(fdepart)+length(ibuy),1)=...
                    iNfarmers(lowprice(ibuy,2));
                subfarminfo=cat(2,LANDINFO{1,t});
                sublandinfo=cat(2,LANDINFO{3,t-1});
                for fd=1:length(fdepart)
                    sublandinfo(subfarminfo==fdepart(fd))=transprice(fdepart(fd),tt);
                    %         PLAND(AGLAYER==fdepart(fd))=transprice(fdepart(fd),tt);
                end
                LANDINFO(3,t)=mat2cell(sublandinfo,NCELLS,1);
                sellrecord(fdepart,1)=tt;
                buyrecord(fdepart,1)=transprice(iNfarmers(lowprice(ibuy,2)),tt);
                iselltime=(sellrecord>0);
                
                ifarmtrans=fdepart;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Update acres sold
            transacres=find(ismember(LANDINFO{1,t},ifarmtrans)==1);
            if isempty(find(transacres,1))==1
                vacland(1,t)=vacland(1,t-1);
                VACLAND(vacland{1,t},t)=1;
            else
                vacland{1,t}=cat(1,[vacland{1,t-1}; transacres]);
                VACLAND(vacland{1,t},t)=1;
            end
            
            subfarminfo=cat(2,LANDINFO{1,t});
            subfarminfo(ismember(subfarminfo,ifarmtrans))=0;
            iNfarmers=unique(subfarminfo);
            iNfarmers=iNfarmers(iNfarmers~=0);
            LANDINFO(1,t+1)=mat2cell(subfarminfo,NCELLS,1);
            
            %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            
            ivac=(VACLAND(:,t) ~= 0);
            indvac=find(ivac==1);
            vcells=zeros([],4);
            subplandinfo=cat(2,LANDINFO{3,t});
            for tl=1:length(indvac)
                subrentEU=zeros([],1);
                subRENTPROJ=RENTPROJ{indvac(tl),t};
                EUlandret(indvac(tl),:)=(subRENTPROJ-((subplandinfo(indvac(tl)).*z(:,1).*...
                    discount)+ccost))./z(:,1);
                EUrankret(indvac(tl),:)=EUlandret(indvac(tl),EUIND(:,indvac(tl)));
            end
            subreturn=cat(2,RETURN{:,t});
            subpland=LANDINFO{3,t};
            for ilayer=1:HT
                subEUprofit=EUrankret(:,ilayer);
                subeuind=EUIND(ilayer,:)';
                posret=(subEUprofit(indvac)>=0);
                subprofit=subreturn(ilayer,:)'-subpland(potentialbuy).*discount;
                iltypes=unique(subeuind);
                iltypes=iltypes(iltypes~=0);
                for ip=1:length(iltypes)
                    Exptprofit(iltypes(ip),t)=mean(subprofit(subeuind(potentialbuy)==iltypes(ip)));
                end
                if isempty(find(indvac(posret),1))==1
                    continue
                else
                    vcells(length(vcells(:,1))+1:length(vcells(:,1))+length(indvac(posret)),:)=...
                        [subEUprofit(indvac(posret)) subeuind(indvac(posret)) ...
                        indvac(posret) ilayer*ones(length(indvac(posret)),1)];
                    ilaycells=(vcells(:,4)==ilayer);
                    rankeucells=sortrows(vcells(ilaycells,:),-1);
                    vcells(ilaycells,:)=rankeucells;
                end
            end
            PROFIT(:,t)=MAXRET(:,t)-(LANDINFO{3,t}.*discount);
            subprofit=PROFIT(:,t);
            subprofit(iurblist)=-1;
            
            rankvaccells=vcells(:,1:3);
            
            Dynmaxprofmap(:,t)=PROFIT(:,t);
            Dynprofltmap(potentialbuy,t)=retind;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%% Profit-Driven Construction decisions %%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            newquota=sum(houseset);
            lotidmark=Nlots(t);
            newhousequota=newhouseset(:,t);
            
            while newquota > 0 || isempty(find(newhousequota,1))==0
                iadded=zeros([],1);
                istartcell=find(rankvaccells(:,1),1);
                iadded(length(iadded)+1:length(iadded)+length(istartcell))=rankvaccells(istartcell,3);
                if isempty(istartcell)==1
                    break
                end
                newlotsize=z(rankvaccells(istartcell,2),1);
                newlottype=rankvaccells(istartcell,2);
                
                newhousequota(newlottype)=newhousequota(newlottype)-max(1/newlotsize,1);
                
                isearchcells=istartcell;
                if newlotsize > 1
                    newlottot=1;
                    while newlottot < newlotsize
                        % Look over 4 cardinal neighbors
                        %%% This will need to be changed when using an irregular
                        %%% grid !!!
                        isnei=zeros(length(iadded),4);
                        numnei=zeros(length(iadded),1);
                        for sc=1:length(iadded)
                            [srows,scols]=ind2sub([NLENGTH NWIDTH],iadded(sc));
                            srnei=(scols+1)*NLENGTH-(NWIDTH-srows);
                            slnei=(scols-1)*NLENGTH-(NWIDTH-srows);
                            supnei=scols*NLENGTH-(NWIDTH-(srows-1));
                            sdnnei=scols*NLENGTH-(NWIDTH-(srows+1));
                            isnei(sc,:)=[srnei slnei supnei sdnnei];
                            dblcount=(ismember(isnei(sc,:),iadded)==0);
                            novac=(ismember(isnei(sc,:),cat(1,Lottype{:,2}))==0);
                            owned=(ismember(isnei(sc,:),ivaclist)==1);
                            bound=(isnei(sc,:) > 0 & isnei(sc,:) <= (NLENGTH*NWIDTH));
                            dblcheck=dblcount.*novac.*bound.*owned;
                            isnei(sc,:)=isnei(sc,:).*dblcheck;
                            numnei(sc)=length(find(isnei(sc,:)~=0));
                        end
                        [maxnum,inei]=max(numnei,[],1);
                        if isempty(find(numnei,1))==1
                            rankvaccells(istartcell,:)=0;
                            iadded=zeros([],1);
                            break
                        end
                        realsnei=zeros([],1);
                        ikeep=(isnei ~= 0);
                        ikeepstar=find(ikeep(inei,:)==1);
                        
                        for in=1:numnei(inei)
                            ivacnei=find(rankvaccells(:,3)==isnei(inei,ikeepstar(in)) ...
                                & rankvaccells(:,2)==newlottype);
                            if isempty(find(ivacnei,1))==1
                                realsnei(length(realsnei)+1)=isnei(inei,ikeepstar(in));
                                continue
                            end
                            realsnei(length(realsnei)+1)=rankvaccells(ivacnei,3);
                        end
                        % Randomize selection of new cells
                        addopts=(1:length(realsnei));
                        addopts=circshift(addopts,[0 round(length(addopts)*rand(1))]);
                        numadd=min(length(realsnei)+newlottot,newlotsize);
                        iadd=(1:numadd-newlottot);
                        newbuildcells=realsnei(addopts(iadd));
                        iadded(length(iadded)+1:length(iadded)+length(iadd),1)=...
                            realsnei(addopts(iadd));
                        isearchcells=newbuildcells;
                        
                        doneswitch=newlottot+length(iadd);
                        if doneswitch == newlotsize
                            totnewbuildcells=iadded;
                            BUILDTIME(totnewbuildcells)=t;
                            BASELAYER(totnewbuildcells)=ones(length(totnewbuildcells),1);
                            VACLAND(totnewbuildcells,t)=0;
                            % %Lottype=[id,location index,lotsize,housesize,ltype,ccost,amlevel,travelcost,buildtime,brokerid]
                            % %lotchoice=[id,location index,ltype,occ/vac,consumer id,residence time,sell price,mitchoice]
                            % %CONINFO=[income,searchtime,consumer_good,housesize,lotsize,proximity,subrisk,occ/vac,utility]
                            nl=length(Lottype(:,1))+1;
                            Lottype{nl,1}=nl;
                            lotchoice{nl,1}=nl;
                            Lottype{nl,2}=iadded;
                            lotchoice{nl,2}=iadded(1);
                            Lottype{nl,3}=z(newlottype,1);
                            Lottype{nl,4}=z(newlottype,2);
                            Lottype{nl,5}=newlottype;
                            lotchoice(nl,3)=Lottype(nl,5);
                            lotchoice{nl,4}=0;
                            lotchoice{nl,5}=0;
                            Lottype{nl,6}=ccost(Lottype{nl,5});
                            Lottype{nl,7}=coastprox{iadded(1)};
                            Lottype{nl,8}=mean(cat(1,travelcost{cat(1,Lottype{nl,2})}));
                            lotchoice{nl,8}=0;
                            Lottype{nl,9}=t;
                            Lottype(nl,10)=num2cell(unique(HBROKER(iadded)),1);
                            
                            newquota=newquota-1;
                            lotidmark=lotidmark+1;
                        end
                        newlottot=newlottot+length(iadd);
                    end
                elseif newlotsize < 1
                    totnewbuildcells=iadded;
                    BUILDTIME(totnewbuildcells)=t;
                    BASELAYER(totnewbuildcells)=1;
                    VACLAND(totnewbuildcells,t)=0;
                    namelength=1/newlotsize;
                    nl=length(Lottype(:,1))+1:length(Lottype(:,1))+namelength;
                    Lottype(nl,1)=num2cell(nl);
                    lotchoice(nl,1)=num2cell(nl);
                    Lottype(nl,2)=num2cell(ones(namelength,1)*iadded);
                    lotchoice(nl,2)=num2cell(ones(namelength,1)*iadded);
                    Lottype(nl,3)=num2cell(ones(namelength,1)*z(newlottype,1));
                    Lottype(nl,4)=num2cell(ones(namelength,1)*z(newlottype,2));
                    Lottype(nl,5)=num2cell(ones(namelength,1)*newlottype);
                    lotchoice(nl,3)=Lottype(nl,5);
                    lotchoice(nl,4)=num2cell(zeros(namelength,1));
                    lotchoice(nl,5)=num2cell(zeros(namelength,1));
                    Lottype(nl,6)=num2cell(ones(namelength,1).*ccost(cat(1,Lottype{nl,5})));
                    Lottype(nl,7)=num2cell(ones(namelength,1)*coastprox{iadded});
                    Lottype(nl,8)=num2cell(ones(namelength,1)*...
                        mean(cat(1,travelcost{cat(1,Lottype{nl,2})})));
                    lotchoice(nl,8)=num2cell(zeros(namelength,1));
                    Lottype(nl,9)=num2cell(ones(namelength,1)*t);
                    Lottype(nl,10)=num2cell(ones(namelength,1)*HBROKER(iadded));
                    
                    newquota=newquota-namelength;
                    lotidmark=lotidmark+namelength;
                else
                    totnewbuildcells=iadded;
                    BUILDTIME(totnewbuildcells)=t;
                    BASELAYER(totnewbuildcells)=1;
                    VACLAND(totnewbuildcells,t)=0;
                    nl=length(Lottype(:,1))+1;
                    Lottype{nl,1}=nl;
                    lotchoice{nl,1}=nl;
                    Lottype{nl,2}=iadded;
                    lotchoice{nl,2}=iadded(1);
                    Lottype{nl,3}=z(newlottype,1);
                    Lottype{nl,4}=z(newlottype,2);
                    Lottype{nl,5}=newlottype;
                    lotchoice(nl,3)=Lottype(nl,5);
                    lotchoice{nl,4}=0;
                    lotchoice{nl,5}=0;
                    Lottype{nl,6}=ccost(Lottype{nl,5});
                    Lottype{nl,7}=coastprox{iadded(1)};
                    Lottype{nl,8}=cat(1,travelcost{cat(1,Lottype{nl,2})});
                    lotchoice{nl,8}=0;
                    Lottype{nl,9}=t;
                    Lottype{nl,10}=HBROKER(iadded);
                    
                    newquota=newquota-1;
                    lotidmark=lotidmark+1;
                end
                
                if newhousequota(newlottype)<=0
                    rankvaccells(rankvaccells(:,2)==newlottype,:)=0;
                end
                
                for nls=1:length(iadded)
                    isameind=find(rankvaccells(:,3)==iadded(nls));
                    rankvaccells(isameind,:)=0;
                end
                
                VACLAND(BASELAYER == 1,t)=0;
                iurblist=find(BASELAYER == 1);
                iagrlist=find(BASELAYER == 0);
                ivac=(VACLAND(:,t) ~= 0);
                ivaclist=find(ivac==1);
                iscape=(SCAPE == 1);
                
                if isempty(find(ivac,1))
                    break
                end
                if isempty(find(newhousequota,1))==1
                    break
                end
            end
            subvacland=vacland{1,t};
            subvacland=subvacland(ismember(subvacland,find(BASELAYER==1))==0);
            vacland{1,t}=subvacland;
            
            Nlots(t)=length(Lottype(:,1));
            damcoef(:,t)=mat2cell(ones(length(Lottype(:,1)),length(CONINFO(:,1))),...
                length(Lottype(:,1)),ones(length(CONINFO(:,1)),1));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%    Farmer Price Prediction    %%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Spatially discounting observed land market prices to use in
            %%% prediction models
            % Calculate Distances to Transactions
            ntrans=length(ifarmtrans);   %numer of land market trans actions, farms sold.
            indtransdist=zeros(NLENGTH,NWIDTH,ntrans);
            transdist=zeros(NLENGTH,NWIDTH);
            Planddistproj=zeros(length(iNfarmers),ntrans);
            subfarminfo=LANDINFO{1,t};
            subplandinfo=LANDINFO{3,t};
            if ntrans > 0
                for nt=1:ntrans
                    tdist2hznnei=10000*ones(NLENGTH,NWIDTH);    %distance to horizontal neighbor from icenter
                    tdist2vrtnei=10000*ones(NLENGTH,NWIDTH);
                    transland=Farminfo{ifarmtrans(nt),2};
                    for tl=1:length(transland(:,1))
                        % Will need to be changed when using irregular sptial
                        % extent !!!
                        [itlrow,itlcol]=ind2sub([NLENGTH,NWIDTH],transland(tl,1));
                        for col=1:NWIDTH
                            tdist2hznnei(1:NLENGTH,col)=min(abs(col-itlcol).*...
                                ones(NLENGTH,1),tdist2hznnei(1:NLENGTH,col));
                        end
                        for row=1:NWIDTH
                            tdist2vrtnei(row,1:NWIDTH)=min(abs(row-itlrow).*...
                                ones(1,NLENGTH),tdist2vrtnei(row,1:NWIDTH));
                        end
                        for col=1:NWIDTH
                            for row=1:NLENGTH
                                indtransdist(row,col,nt)=sqrt(tdist2hznnei(row,col)^2+tdist2vrtnei(row,col)^2);
                            end
                        end
                    end
                    subfarminfo(transland(:,1))=0;
                    LANDINFO(1,t)=mat2cell(subfarminfo,length(subfarminfo),1);
                    iremain=(subfarminfo ~= 0);
                    iNfarmers=unique(subfarminfo(iremain));
                    %Calculate land price gradient coefficients using genetic algorithm
                    subtransdist=indtransdist(:,:,nt);
                    distcoeff(ifarmtrans,:)=0;
                    for nf=1:length(iNfarmers)
                        rc=zeros([],1);
                        avgtransdist=mean(subtransdist(subfarminfo==iNfarmers(nf)));
                        coeffmark=(Paskland(iNfarmers(nf),t)-mean(subplandinfo(transland(:,1))))/...
                            avgtransdist;
                        fitness(iNfarmers(nf),:,t)=fitness(iNfarmers(nf),:,t-1)+...
                            abs(distcoeff(iNfarmers(nf),:)-coeffmark);
                        fitsort=sort(fitness(iNfarmers(nf),:,t),'ascend');
                        stratcount=1;
                        for x=1:nextgen
                            numstrat=find(fitness(iNfarmers(nf),:,t)==fitsort(x));
                            if length(numstrat) > 1
                                isurvive(x,1)=numstrat(stratcount);
                                stratcount=stratcount+1;
                            else
                                isurvive(x,1)=numstrat(stratcount);
                            end
                            if stratcount > length(numstrat)
                                stratcount=1;
                            end
                        end
                        distcoeff(iNfarmers(nf),1:nextgen)=distcoeff(iNfarmers(nf),isurvive);
                        for xx=1:nextgen
                            irc=xx*recombo/nextgen;
                            rc(length(rc)+1:irc,1)=distcoeff(iNfarmers(nf),xx)+rand(1,4);
                        end
                        distcoeff(iNfarmers(nf),nextgen+1:nextgen+recombo)=rc';
                        distcoeff(iNfarmers(nf),nextgen+recombo+1:NUMMODELDIST)=...
                            mincoeff+(maxcoeff-mincoeff)*rand(1,NUMMODELDIST-...
                            (nextgen+recombo));
                    end
                    Planddistproj(iNfarmers,nt)=distcoeff(iNfarmers,1).*avgtransdist+...
                        mean(subplandinfo(transland(:,1)));
                end
                for nnw=1:NWIDTH
                    for nnl=1:NLENGTH
                        transdist(nnl,nnw)=min(indtransdist(nnl,nnw,1:ntrans));
                    end
                end
                Plandproj(iNfarmers,t)=mean(Planddistproj(iNfarmers,:),2);
            else
                Plandproj(iNfarmers,t)=(1-DELTA)*Plandproj(iNfarmers,t-1)+DELTA*...
                    mean([wtaland(iNfarmers,t) zeta*wtpland(iNfarmers,t)],2);
            end
            %%% Farmers' Prediction Models %%%
            for nf=1:length(iNfarmers)
                ilandclass1=find(landmodel(iNfarmers(nf),:)==1);
                ilandclass2=find(landmodel(iNfarmers(nf),:)==2);
                ilandclass3=find(landmodel(iNfarmers(nf),:)==3);
                ilandclass4=find(landmodel(iNfarmers(nf),:)==4);
                ilandclass5=find(landmodel(iNfarmers(nf),:)==5);
                ilandclass6=find(landmodel(iNfarmers(nf),:)==6);
                
                landerror(iNfarmers(nf),:) = (1-DELTA)*landerror(iNfarmers(nf),:)+...
                    DELTA*abs(Plandproj(iNfarmers(nf),t)-landproj(iNfarmers(nf),:));
                [landbest,ilandbest] = min(landerror(iNfarmers(nf),:),[],2);
                landbestSAVE(iNfarmers(nf),t) = landbest;
                ilandbestSAVE(iNfarmers(nf),t) = ilandbest;
                landprojSAVE(iNfarmers(nf),t) = landproj(iNfarmers(nf),ilandbest);
                landmodelSAVE(iNfarmers(nf),t) = landmodel(iNfarmers(nf),ilandbest);
                
                for i = 1:FARMNUMCLASS
                    if i == 1
                        % mirror models
                        landproj(iNfarmers(nf),ilandclass1) = Plandproj(iNfarmers(nf),t)+(1-aa...
                            (iNfarmers(nf),ilandclass1)).*(0.5*Plandproj(iNfarmers(nf),t)-...
                            (Plandproj(iNfarmers(nf),t)-Plandproj(iNfarmers(nf),t-1)));
                    elseif i == 2
                        % mean model
                        for jl = 1:length(ilandclass2)
                            landproj(iNfarmers(nf),ilandclass2(jl)) = mean(Plandproj(iNfarmers(nf),...
                                t:-1:(t-aa(iNfarmers(nf),ilandclass2(jl)))));
                        end
                    elseif i == 3
                        %cycle model
                        landproj(iNfarmers(nf),ilandclass3) = Plandproj(iNfarmers(nf),...
                            t-round(max(1,aa(iNfarmers(nf),ilandclass3))));
                    elseif i == 4
                        % projection model
                        for jl = 1:length(ilandclass4)
                            %Nonlinear Forecast
                            indata=Plandproj(iNfarmers(nf),t-(1+aa(iNfarmers(nf),ilandclass4(jl))):t);
                            pcoef=polyfit(1:length(indata),indata,1);
                            pline=pcoef(1).*(1:length(indata)+1)+pcoef(2);
                            landproj(iNfarmers(nf),ilandclass4(jl))=pline(length(pline));
                        end
                    elseif i == 5
                        % rescale model
                        landproj(iNfarmers(nf),ilandclass5) = aa(iNfarmers(nf),...
                            ilandclass5)*Plandproj(iNfarmers(nf),t);
                    elseif i == 6
                        % local(0) or regional(1) trends
                        ilandlocal=(aa(iNfarmers(nf),ilandclass6)==0);
                        ilandreg=(aa(iNfarmers(nf),ilandclass6)==1);
                        if isempty(iNfarmers)==1
                            break
                        end
                        landproj(iNfarmers(nf),ilandclass6(ilandlocal)) = Plandproj...
                            (iNfarmers(nf),t).*(1+1/length(iNfarmers));
                    end
                end
                wtaland(iNfarmers(nf),t+1)=max(landproj(iNfarmers(nf),ilandbest),...
                    farmretinfo(iNfarmers(nf)));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%    Consumers' Choices    %%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%% Add new consumers %%%%%%%%%%
            newpop=ceil(length(CONINFO(:,1))*POPGROW);
            inewpop=length(CONINFO(:,1))+1:length(CONINFO(:,1))+newpop;
            % %CONINFO=[income,searchtime,consumer_good,housesize,lotsize,proximity,subrisk,occ/vac,utility]
            CONINFO(inewpop,1)=num2cell(...
                min(max(normrnd(avgincome,stdincome,newpop,1),minwage),maxwage));
            CONINFO(inewpop,2)=num2cell(...
                t+ceil(searchtimemin+(searchtimemax-searchtimemin)*rand(newpop,1)));
            
            newincome=cat(1,CONINFO{inewpop,1});
            newincome3=(newincome >= lowminwage & newincome <= lowmaxwage);     %income 1 2 3 = hi mid low wages
            newincome2=(newincome >= midminwage & newincome <= midmaxwage);
            newincome1=(newincome >= himinwage & newincome <= himaxwage);
            housepref(inewpop,1)=0;
            housepref(inewpop(newincome1),1)=HIBETA(1)+(HIBETA(2)-HIBETA(1))*rand(length(find(newincome1==1)),1);
            housepref(inewpop(newincome2),1)=MIDBETA(1)+(MIDBETA(2)-MIDBETA(1))*rand(length(find(newincome2==1)),1);
            housepref(inewpop(newincome3),1)=LOWBETA(1)+(LOWBETA(2)-LOWBETA(1))*rand(length(find(newincome3==1)),1);
            
            CONINFO(inewpop,3)=num2cell(1-housepref(inewpop));
            CONINFO(inewpop,6)=num2cell((0.1+(0.9-0.1)*rand(length(housepref(inewpop)),1)).*housepref(inewpop));
            CONINFO(inewpop,4)=num2cell(rand(length(housepref(inewpop)),1).*...)
                (housepref(inewpop)-cat(1,CONINFO{inewpop,6}))/2);
            CONINFO(inewpop,5)=num2cell(housepref(inewpop)-(cat(1,CONINFO{inewpop,4})+cat(1,CONINFO{inewpop,6})));
            %     CONINFO(inewpop,4)=num2cell((0.1+(0.9-0.1)*rand(length(housepref(inewpop)),1)).*housepref(inewpop));
            %     CONINFO(inewpop,5)=num2cell(rand(length(housepref(inewpop)),1).*...)
            %         (housepref(inewpop)-cat(1,CONINFO{inewpop,4}))/2);
            %     CONINFO(inewpop,6)=num2cell(housepref(inewpop)-(cat(1,CONINFO{inewpop,4})+cat(1,CONINFO{inewpop,5})));
            CONINFO(inewpop,7)=num2cell(ones(newpop,1));
            CONINFO(inewpop,8)=num2cell(zeros(newpop,1));
            CONINFO(inewpop,9)=num2cell(zeros(newpop,1));
            damcoef(inewpop,t)=mat2cell(ones(length(Lottype(:,1)),newpop),length(Lottype(:,1)),ones(newpop,1));
            
            %%%%%%%%%%%%%
            
            minprof=zeros(Nlots(t),1);
            vachouse=find(cat(1,lotchoice{:,4})==0);
            moveouts=find(cat(1,lotchoice{:,6})==t);
            inewlots=unique([vachouse; moveouts]);
            inewcon=find(cat(1,CONINFO{:,8})==0);
            for il=1:length(inewlots)
                if Lottype{inewlots(il),9} < t
                    continue
                end
                lotind=cat(1,Lottype{inewlots(il),2});
                lotlocate(length(lotlocate(:,1))+1:length(lotlocate(:,1))+length(lotind),:)=...
                    [ones(length(lotind),1)*inewlots(il) lotind];
                brokerassign=Lottype{inewlots(il),10};
                brkrlocate(length(brkrlocate(:,1))+1:length(brkrlocate(:,1))+length(brokerassign),:)=...
                    [ones(length(brokerassign),1)*inewlots(il) brokerassign];
            end
            
            bt=cat(1,Lottype{inewlots,9});
            Paskhouse(inewlots(bt < t))=cat(1,lotchoice{inewlots(bt < t),7});
            if isempty(find(bt == t,1)) == 0
                newrentinfo=cat(2,RENTPROJ{cat(1,lotchoice{inewlots(bt == t),2}),t});
                newrentind=sub2ind(size(newrentinfo),cat(1,lotchoice{inewlots(bt==t),3}),...
                    (1:length(newrentinfo(1,:)))');
                Paskhouse(inewlots(bt == t))=newrentinfo(newrentind);
            end
            BIDLEVEL(inewlots(bt ==t))=num2cell(zeros(length(inewlots(bt == t)),1));
            AVGUTIL(inewlots(bt ==t))=num2cell(zeros(length(inewlots(bt == t)),1));
            lotchoice(inewlots,7)=num2cell(Paskhouse(inewlots));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%   House Market   %%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% HouseMarketDynamic_coast_batch
            wtpcon=zeros(length(CONINFO(:,1)),Nlots(t));
            wtpconstar=zeros(length(CONINFO(:,1)),Nlots(t));
            wtbcon=zeros(length(CONINFO(:,1)),Nlots(t));
            Rn=zeros(length(CONINFO(:,1)),Nlots(t));
            Phousebid=zeros(length(CONINFO(:,1)),Nlots(t));
            subnhouselook=zeros(length(CONINFO(:,1)),Nlots(t));
            notherbuyers=zeros(length(CONINFO(:,1)),Nlots(t));
            iotherbuyers=zeros(length(CONINFO(:,1)),Nlots(t));
            housemp=zeros(length(CONINFO(:,1)),1);
            EUmit=zeros(length(CONINFO(:,1)),Nlots(t));
            EUnomit=zeros(length(CONINFO(:,1)),Nlots(t));
            U=zeros(length(CONINFO(:,1)),Nlots(t));
            Unorm=zeros(length(CONINFO(:,1)),Nlots(t));
            maxUset=zeros(length(CONINFO(:,1)),LTYPE*HTYPE);
            exptCmit=zeros(length(CONINFO(:,1)),Nlots(t));
            exptCnomit=zeros(length(CONINFO(:,1)),Nlots(t));
            mitchoice=zeros(length(CONINFO(:,1)),Nlots(t));
            exptcost=zeros(2,Nlots(t),length(CONINFO(:,1)));
            % %Lottype=[id,location index,lotsize,housesize,ltype,ccost,amlevel,travelcost,buildtime,brokerid]
            % %lotchoice=[id,location index,ltype,occ/vac,consumer id,residence time,sell price,mitchoice]
            % %CONINFO=[income,searchtime,consumer_good,housesize,lotsize,proximity,subrisk,occ/vac,utility]
%             Cdam(inewlots,t)=mat2cell(sum(Pimpact(cat(1,lotchoice{inewlots,2})).*...
%                 repmat(housedam,length(inewlots),1).*repmat(Paskhouse(inewlots),1,5),2));
            if isempty(find(inewlots,1))==0
                Cdam(inewlots,t)=num2cell(Paskhouse(inewlots).*...
                    (housedam(cat(1,lotchoice{inewlots,2})))./100);
%                 Cdam(inewlots,t)=num2cell(mean((Paskhouse(inewlots)*housedam).*...
%                     repmat(mean(Psevere,1),length(inewlots),1).*...
%                     repmat(cat(1,Pflood{cat(1,lotchoice{inewlots,2})}),1,5),2));
            end
            for c=1:length(inewcon)
                if isempty(find(inewlots,1))==1
                    break
                end
                wtpcon(inewcon(c),inewlots)=(CONINFO{inewcon(c),1}-...
                    cat(1,travelcost{cat(1,lotchoice{inewlots,2})})-cat(1,Cdam{inewlots,t}))*...
                    (CONINFO{inewcon(c),4}+CONINFO{inewcon(c),5}+CONINFO{inewcon(c),6});
                EUnomit(inewcon(c),inewlots)=sum(Psevere)*...
                    (((CONINFO{inewcon(c),1}-cat(1,travelcost{cat(1,lotchoice{inewlots,2})})-...
                    Paskhouse(inewlots)-cat(1,Cdam{inewlots,t})).^CONINFO{inewcon(c),3}).*(cat(1,Lottype{inewlots,4}).^...
                    CONINFO{inewcon(c),4}).*(cat(1,Lottype{inewlots,3}).^CONINFO{inewcon(c),5}).*...
                    (cat(1,Lottype{inewlots,7}).^CONINFO{inewcon(c),6}))+...
                    (1-sum(Psevere))*...
                    (((CONINFO{inewcon(c),1}-cat(1,travelcost{cat(1,lotchoice{inewlots,2})})-...
                    Paskhouse(inewlots)).^CONINFO{inewcon(c),3}).*(cat(1,Lottype{inewlots,4}).^...
                    CONINFO{inewcon(c),4}).*(cat(1,Lottype{inewlots,3}).^CONINFO{inewcon(c),5}).*...
                    (cat(1,Lottype{inewlots,7}).^CONINFO{inewcon(c),6}));
%                 U(inewcon(c),inewlots)=((CONINFO{inewcon(c),1}-cat(1,travelcost{cat(1,lotchoice{inewlots,2})})-...
%                     Paskhouse(inewlots)).^CONINFO{inewcon(c),3}).*(cat(1,Lottype{inewlots,4}).^...
%                     CONINFO{inewcon(c),4}).*(cat(1,Lottype{inewlots,3}).^CONINFO{inewcon(c),5}).*...
%                     (cat(1,Lottype{inewlots,7}).^CONINFO{inewcon(c),6});
                %%%%%%%%%% Expected Costs %%%%%%%%%%%%
                %     subdamcoef=damcoef{inewcon(c),t};
                %     exptCmit(inewcon(c),inewlots)=subdamcoef(inewlots).*cat(1,Pflood{cat(1,lotchoice{inewlots,2})}).*...
                %         (Cmit(erun)+miteff(erun)*Cdam(erun)*Paskhouse(inewlots))+(1-subdamcoef(inewlots)).*...
                %         (1-cat(1,Pflood{cat(1,lotchoice{inewlots,2})}))*Cmit(erun);
                %     exptCnomit(inewcon(c),inewlots)=subdamcoef(inewlots).*cat(1,Pflood{cat(1,lotchoice{inewlots,2})}).*...
                %         (Cdam(erun)*Paskhouse(inewlots));
                %     exptcost(:,inewlots,inewcon(c))=[exptCmit(inewcon(c),inewlots); exptCnomit(inewcon(c),inewlots)];
                %     wtpcon(inewcon(c),inewlots)=(CONINFO{inewcon(c),1}-cat(1,travelcost{cat(1,lotchoice{inewlots,2})}))*...
                %         (CONINFO{inewcon(c),4}+CONINFO{inewcon(c),5}+CONINFO{inewcon(c),6});
                %     EUmit(inewcon(c),inewlots)=(max(CONINFO{inewcon(c),1}-cat(1,travelcost{cat(1,lotchoice{inewlots,2})})-...
                %         Paskhouse(inewlots)-exptCmit(inewcon(c),inewlots)',0).^CONINFO{inewcon(c),3}).*...
                %         (cat(1,Lottype{inewlots,4}).^CONINFO{inewcon(c),4}).*(cat(1,Lottype{inewlots,3}).^...
                %         CONINFO{inewcon(c),5}).*(cat(1,Lottype{inewlots,7}).^CONINFO{inewcon(c),6});
                %     EUnomit(inewcon(c),inewlots)=(max(CONINFO{inewcon(c),1}-cat(1,travelcost{cat(1,lotchoice{inewlots,2})})-...
                %         Paskhouse(inewlots)-exptCnomit(inewcon(c),inewlots)',0).^CONINFO{inewcon(c),3}).*...
                %         (cat(1,Lottype{inewlots,4}).^CONINFO{inewcon(c),4}).*(cat(1,Lottype{inewlots,3}).^...
                %         CONINFO{inewcon(c),5}).*(cat(1,Lottype{inewlots,7}).^CONINFO{inewcon(c),6});
                %     [Umax,mitcheck]=max([EUmit(c,:); EUnomit(c,:)]);
                    U(inewcon(c),inewlots)=EUnomit(inewcon(c),inewlots);
                %     mitchoice(inewcon(c),inewlots)=2*ones(1,length(EUnomit(inewcon(c),inewlots)));
                % %     U(c,:)=Umax;
                % %     mitchoice(c,:)=mitcheck;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %     ihousein=(cat(1,BIDLEVEL{inewlots}).*Paskhouse(inewlots) < wtpcon(inewcon(c),inewlots)' & ...
                %         U(inewcon(c),inewlots)' > CONINFO{inewcon(c),3}*max(U(inewcon(c),inewlots)));
                ihousein=cat(1,BIDLEVEL{inewlots}).*Paskhouse(inewlots) < wtpcon(inewcon(c),inewlots)';
                ihouseout=find(inewlots(ihousein) == 0);
                
                Unorm(inewcon(c),inewlots(ihousein))=U(inewcon(c),inewlots(ihousein))./...
                    max(U(inewcon(c),inewlots(ihousein)));
                
                if isempty(find(inewlots(ihousein),1))==1
                    Rn(inewcon(c),inewlots)=zeros(length(inewcon(c)),length(inewlots));
                    wtbcon(inewcon(c),inewlots)=zeros(length(inewcon(c)),length(inewlots));
                else
                    Rn(inewcon(c),inewlots(ihousein))=Paskhouse(inewlots(ihousein)).*...
                        Unorm(inewcon(c),inewlots(ihousein))';
                    wtbcon(inewcon(c),inewlots(ihousein))=min(wtpcon(inewcon(c),inewlots(ihousein))-...
                        (Paskhouse(inewlots(ihousein))'-Rn(inewcon(c),inewlots(ihousein))),...
                        wtpcon(inewcon(c),inewlots(ihousein)));
                end
                subnhouselook(inewcon(c),inewlots)=(ihousein==1);
            end
            % if t>=12 && isempty(find(inewlots,1))==0
            %     keyboard
            % end
            nhouselook=(subnhouselook == 1);
            for nl=1:length(inewlots)
                notherbuyers(nhouselook(:,inewlots(nl)),inewlots(nl))=...
                    find(nhouselook(:,inewlots(nl))==1);
            end
            
            for c=1:length(inewcon)
                if isempty(find(nhouselook(inewcon(c),inewlots),1))==1
                    continue
                else
                    nhouses=length(find(nhouselook(inewcon(c),inewlots)==1));
                    subbuyers=unique(notherbuyers(:,nhouselook(inewcon(c),inewlots)));
                    subbuyers=subbuyers(subbuyers~=0);
                    nbuyers=length(subbuyers);
                    housemp(inewcon(c))=(nbuyers-nhouses)/(nbuyers+nhouses);
                    
                    if housemp(inewcon(c)) >= 0
                        Phousebid(inewcon(c),nhouselook(inewcon(c),:))=min(max(Rn(inewcon(c),...
                            nhouselook(inewcon(c),:))+(wtbcon(inewcon(c),nhouselook(inewcon(c),:))-...
                            Rn(inewcon(c),nhouselook(inewcon(c),:)))*housemp(inewcon(c)),...
                            Rn(inewcon(c),nhouselook(inewcon(c),:))),wtbcon(inewcon(c),nhouselook(inewcon(c),:)));
                    elseif housemp(inewcon(c)) < 0
                        Phousebid(inewcon(c),nhouselook(inewcon(c),:))=min(Rn(inewcon(c),...
                            nhouselook(inewcon(c),:))+Rn(inewcon(c),nhouselook(inewcon(c),:)).*...
                            (1./(wtbcon(inewcon(c),nhouselook(inewcon(c),:))-Rn(inewcon(c),...
                            nhouselook(inewcon(c),:))))*housemp(inewcon(c)),...
                            wtbcon(inewcon(c),nhouselook(inewcon(c),:)));
                    end
                end
            end
            
            avghousemp(t)=mean(housemp);
            openhouse=length(inewlots);
            subPhousebid=Phousebid;
            subU=U;
            iunderbid=zeros(size(Phousebid));
            sublandinfo=cat(2,LANDINFO{3,t});
            newbidlevel=cat(1,lotchoice{inewlots,7})./((cat(1,Lottype{inewlots,6})+...
                discount*sublandinfo(cat(1,lotchoice{inewlots,2})).*...
                cat(1,Lottype{inewlots,3})).^(1+(t-cat(1,Lottype{inewlots,9}))./TMAX));
            for nl=1:length(inewlots)
                if Lottype{inewlots(nl),9} == t
                    iunderbid(:,inewlots(nl))=(subPhousebid(:,inewlots(nl)) < Paskhouse...
                        (inewlots(nl))/(1+discount));
                elseif Lottype{inewlots(nl),9} < t
                    iunderbid(:,inewlots(nl))=(subPhousebid(:,inewlots(nl)) < Paskhouse...
                        (inewlots(nl))*newbidlevel(nl));
                end
                subPhousebid((iunderbid(:,inewlots(nl))==1),inewlots(nl))=0;
                %     if BUILDTIME(lotchoice(inewlots(nl),2)) == t
                %         iunderbid(:,inewlots(nl))=(subPhousebid(:,inewlots(nl)) < Paskhouse...
                %             (inewlots(nl))/(1+discount));
                %     elseif BUILDTIME(lotchoice(inewlots(nl),2)) < t
                %         iunderbid(:,inewlots(nl))=(subPhousebid(:,inewlots(nl)) < Paskhouse...
                %             (inewlots(nl))*BIDLEVEL(lotchoice(inewlots(nl),2)));
                %     end
                %     subPhousebid((iunderbid(:,inewlots(nl))==1),inewlots(nl))=0;
            end
            while openhouse > 0
                if isempty(find(subPhousebid > 0,1))==1
                    break
                end
                wincon=zeros(1,length(inewlots));
                istillopen=find(cat(1,lotchoice{inewlots,4})==0);
                [maxbid,imaxbid]=max(subPhousebid(:,inewlots),[],1);
                for nl=1:length(inewlots)     %find highest bid for each house in this round
                    iwincon=find(subPhousebid(:,inewlots(nl))==maxbid(nl));
                    if maxbid(nl) <= 0
                        continue
                    end
                    if length(unique(iwincon)) > 1
                        icon=ceil(length(iwincon)*rand(1));
                        wincon(nl)=iwincon(icon);
                    else
                        wincon(nl)=iwincon;
                    end
                end
                conset=unique(wincon);      %Highest bidders at the moment
                conset=conset(conset~=0);
                randorder=randperm(length(conset));
                conset=conset(randorder);
                for cs=1:length(conset)
                    ilotmatch=find(wincon==conset(cs));
                    %         subexptcost=exptcost(:,:,conset(cs));
                    %         iexptcost=sub2ind(size(subexptcost),mitchoice(conset(cs),...
                    %             inewlots(ilotmatch)),inewlots(ilotmatch)');
                    %         uset=(CONINFO{conset(cs),1}-cat(1,travelcost{cat(1,lotchoice{inewlots(ilotmatch),2})})-...
                    %            subPhousebid(conset(cs),inewlots(ilotmatch))'-subexptcost(iexptcost)'.^CONINFO{conset(cs),3}).*...
                    %            (cat(1,Lottype{inewlots(ilotmatch),4}).^CONINFO{conset(cs),4}).*...
                    %            (cat(1,Lottype{inewlots(ilotmatch),3}).^CONINFO{conset(cs),5}).*...
                    %            (cat(1,Lottype{inewlots(ilotmatch),7}).^CONINFO{conset(cs),6});
                    uset=((CONINFO{conset(cs),1}-cat(1,travelcost{cat(1,lotchoice{inewlots(ilotmatch),2})})-...
                        subPhousebid(conset(cs),inewlots(ilotmatch))').^CONINFO{conset(cs),3}).*...
                        (cat(1,Lottype{inewlots(ilotmatch),4}).^CONINFO{conset(cs),4}).*...
                        (cat(1,Lottype{inewlots(ilotmatch),3}).^CONINFO{conset(cs),5}).*...
                        (cat(1,Lottype{inewlots(ilotmatch),7}).^CONINFO{conset(cs),6});
                    ilotid=find(uset==max(uset));
                    
                    if length(ilotid) > 1
                        ipick=find(subPhousebid(conset(cs),inewlots(ilotmatch(ilotid)))==...
                            max(subPhousebid(conset(cs),inewlots(ilotmatch(ilotid)))));
                        if length(ipick) > 1
                            ipick=ipick(ceil(length(ipick)*rand(1)));
                        end
                        lotid=inewlots(ilotmatch(ilotid(ipick)));
                    else
                        ipick=1;
                        lotid=inewlots(ilotmatch(ilotid));
                    end
                    
                    conid=conset(cs);
                    CONINFO{conid,8}=1;
                    CONINFO{conid,9}=uset(ilotid(ipick));
                    %         CONINFO{conid,9}=uset(conid,lotid);
                    lotchoice{lotid,4}=1;
                    lotchoice{lotid,5}=conid;
                    lotchoice{lotid,6}=max(ceil(t+avgrestime/2+normrnd(avgrestime/...
                        2,stdrestime/2,1,1)),t+1);
                    lotchoice{lotid,7}=subPhousebid(conid,lotid);
                    %         lotchoice{lotid,8}=mitchoice(conid,lotid);
                    lotchoice{lotid,8}=1;
                    
                    subPhousebid(conid,:)=0;
                    subPhousebid(:,lotid)=0;
                    subU(conid,:)=0;
                    subU(:,lotid)=0;
                    MITIGATE(cat(1,Lottype{lotid,2}))=num2cell(rem(lotchoice{lotid,8},2));
                    openhouse=openhouse-1;
                end
            end
            % VACANT HOUSES
            conlist=(1:length(CONINFO(:,1)))';
            ifilled=find(cat(1,lotchoice{:,4})==1);
            istillvac=find(cat(1,lotchoice{:,4})==0);
            popin=cat(1,lotchoice{ifilled,5});
            popout=conlist(~ismember(conlist,popin));
            lotchoice(istillvac,6)=num2cell(ones(length(istillvac),1)*t+1);
            lotchoice(istillvac,8)=num2cell(zeros(length(istillvac),1));
            
            % %Lottype=[id,location index,lotsize,housesize,ltype,ccost,amlevel,travelcost,buildtime,brokerid]
            % %lotchoice=[id,location index,ltype,occ/vac,consumer id,residence time,sell price,mitchoice]
            % %CONINFO=[income,searchtime,consumer_good,housesize,lotsize,proximity,subrisk,occ/vac,utility]
for sv=1:length(istillvac)
   if Lottype{istillvac(sv),9} == t 
       lotchoice{istillvac(sv),7}=Paskhouse(istillvac(sv))./(1+discount);
   elseif Lottype{istillvac(sv),9} < t
       lotchoice{istillvac(sv),7}=max(Paskhouse(istillvac(sv))-...
           discount*Lottype{istillvac(sv),6}*(t-Lottype{istillvac(sv),9}),0);
   end
end

for kk=1:length(Lottype(:,1))
    RENT(cat(1,Lottype{kk,2}),t)=ones(length(Lottype{kk,2}),1)*cat(1,lotchoice{kk,7});
%    isamecell=ismember(Lottype(:,2),cellinfo(kk,2));
%    if length(find(isamecell==1)) > 1
%        RENT(cellinfo(kk,2))=mean(con2lot(Lottype(isamecell,1),1));
%    elseif length(find(isamecell==1)) == 1
%        RENT(cellinfo(kk,2))=con2lot(cellinfo(kk,1),1);
%    end
end


% %%%%%%%% Utility check %%%%%%%%%%%
% realulot=zeros(length(Income),3);
% conlist=(1:length(Income))';
% inhouselist=conlist(popinhouse);
% maxulot=zeros(1,[]);
% imaxulot=zeros(1,[]);
% for c=1:length(Income)
%     [maxulot(c),imaxulot(c)]=max(U(c,:),[],2);   
% end
% for c=1:length(inhouselist)
%     ireallot=find(con2lot(:,2)==inhouselist(c));
%     realulot(inhouselist(c),1:3)=[ireallot lotchoice(ireallot,5) U(inhouselist(c),ireallot)];
% end
% maxuset=[imaxulot' lotchoice(imaxulot,5) maxulot'];
% fullset=[conlist realulot maxuset];
% utildiff(t)=mean(fullset(fullset(:,4)~=0,7)-fullset(fullset(:,4)~=0,4));
% pctutildiff(t)=mean(fullset(fullset(:,4)~=0,4)./fullset(fullset(:,4)~=0,7));
            
            RENTGRAD(lotlocate(:,2))=RENT(lotlocate(:,2),t).*1./cat(1,Lottype{lotlocate(:,1),3});
            %<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
            popinhouse=(cat(1,CONINFO{:,8})==1);
            
            oldincome(t)=mean(cat(1,CONINFO{popout,1}));
            
            Nlots(t+1)=Nlots(t);
            for lt=1:HT
                vacstats(lt,t)=length(find(cat(1,lotchoice{istillvac,3})==lt));
            end
            devcells(1,t)=length(iurblist);
            devcells(2,t)=devcells(1,t)/NCELLS;
            
            vacantlots(t)=length(istillvac);
            leftoverpop(t)=length(popout);
            
            vacrate(t)=vacantlots(t)/Nlots(t);
            nohouserate(t)=leftoverpop(t)/POP(t);
            
            %Consumers leave area is searchtime is met
            iconleave=cat(1,lotchoice{cat(1,lotchoice{:,6})==t,5});
            OutIncome(length(OutIncome)+1:length(OutIncome)+length(iconleave),1)=...
                cat(1,CONINFO{iconleave,1});
            CONINFO(iconleave,1)=num2cell(zeros(length(iconleave),1));
            Nconsumers=length(find(cat(1,CONINFO{:,1})~=0));
            subrealret=zeros(NCELLS,1);
            subexptrent=zeros(NCELLS,1);
            subnewlt=zeros(NCELLS,1);
            subnewbid=zeros(NCELLS,1);
            subrealexptret=zeros(NCELLS,1);
            subplandinfo=LANDINFO{3,t};
            
            if isempty(find(inewlots,1))==0
                subexptrent(cat(1,lotchoice{inewlots,2}))=cat(1,lotchoice{inewlots,7})-...
                    Paskhouse(cat(1,lotchoice{inewlots,1}));
                subrealret(cat(1,lotchoice{inewlots,2}))=cat(1,lotchoice{inewlots,7})-...
                    (subplandinfo(cat(1,lotchoice{inewlots,2})).*z(cat(1,lotchoice{inewlots,3}),1)*...
                    discount+ccost(cat(1,lotchoice{inewlots,3})));
                subrealexptret(cat(1,lotchoice{inewlots,2}))=Paskhouse(inewlots)-...
                    (subplandinfo(cat(1,lotchoice{inewlots,2})).*z(cat(1,lotchoice{inewlots,3}),1)*...
                    discount+ccost(cat(1,lotchoice{inewlots,3})));
                subnewbid(cat(1,lotchoice{inewlots,2}))=cat(1,lotchoice{inewlots,7})./...
                    Paskhouse(inewlots);
            end
            
            Exptrentdiff(:,t)=subexptrent;
            Realreturn(:,t)=subrealret;
            Realexptret(:,t)=subrealexptret;
            vac_land(t)=sum(subplandinfo(ivac)*discount);
            for lt=1:HT
                avgrent(lt,t)=mean(RENT(cat(1,Lottype{cat(1,Lottype{:,5})==lt,2}),t));
                ifindlt=(cat(1,Lottype{inewlots,5})==lt);
                Avgexptdiff(lt,t)=mean(subexptrent(cat(1,lotchoice{inewlots(ifindlt),2})));
                Realavgret(lt,t)=mean(subrealret(cat(1,lotchoice{inewlots(ifindlt),2})));
                Realavgexptret(lt,t)=mean(subrealexptret(cat(1,lotchoice{inewlots(ifindlt),2})));
                Avgnewbid(lt,t)=mean(subnewbid(cat(1,lotchoice{inewlots(ifindlt),2})));
                
                ilt=find(cat(1,lotchoice{istillvac,5})==lt);
                if isempty(find(ilt,1))==1
                    vac_ccost(lt,t)=0;
                else
                    vac_ccost(lt,t)=sum(discount*(cat(1,Lottype{istillvac(ilt),6})+...
                        subplandinfo(cat(1,Lottype{istillvac(ilt),2})).*...
                        z(cat(1,Lottype{istillvac(ilt),5}),1)));
                end
                vac_rent(lt,t)=sum(Paskhouse(istillvac(ilt)));
                
                ilt_t=find(cat(1,Lottype{:,5}) == lt & cat(1,Lottype{:,9}) == t);
                if isempty(find(ilt_t,1))==1
                    profits(lt,t)=0;
                    budget_lt(lt,t)=0;
                else
                    profits(lt,t)=sum(cat(1,lotchoice{ilt_t,7})-(subplandinfo(cat(1,lotchoice{ilt_t,2})).*...
                        z(cat(1,lotchoice{ilt_t,3}),1)*discount+cat(1,Lottype{ilt_t,6})));
                    budget_lt(lt,t)=profits(lt,t)-vac_ccost(lt,t);
                end
            end
            Newbidlevel(:,t)=subnewbid;
            carrycost(t)=sum(vac_ccost(:,t))+vac_land(t);
            
            BUDGET(t)=BUDGET(t-1)+sum(budget_lt(:,t))-vac_land(t);
            
            %     figure(1)
            %     surf(reshape(LOTTYPE(:,t),NLENGTH,NWIDTH));
            %     axis ij;
            %     view(0,90);
            %     title(sprintf('Lot Types, t=%d',t))
            %     set(gca,'clim',[1 HT])
            %     colorbar
            %     MLT(t)=getframe(gcf);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%    Brokers' Price Projections    %%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            warning('off','all');
            for ibr=1:Nbrokers
                bhood=find(HBROKER==ibr);
                brokerind=unique(brkrlocate(cat(1,Lottype{:,10})==ibr,1));
                iblots=cell2mat(lotchoice(brokerind,[1 3 4 5]));
                
                if isempty(iblots) == 1
                    continue
                else
                    nbids=zeros([],1);
                    for b=1:length(iblots(:,1))
                        bidover=(Phousebid(:,iblots(b,1)) > 0);
                        nbids(b,1)=length(find(bidover==1));
                    end
                    
                    % Calculate average utilities for rent projection of unknown
                    isnotvac=(iblots(:,3)==1);
                    if isempty(find(isnotvac,1))==1
                        continue
                    else
                        brkravgstats(ibr,:)=[median(cat(1,CONINFO{iblots(isnotvac,4),1})) ...
                            median(cat(1,CONINFO{iblots(isnotvac,4),3})) ...
                            mean(cat(1,CONINFO{iblots(isnotvac,4),4})) ...
                            mean(cat(1,CONINFO{iblots(isnotvac,4),5})) ...
                            mean(cat(1,CONINFO{iblots(isnotvac,4),6}))];
                        AVGUTIL(brokerind(isnotvac))=CONINFO(iblots(isnotvac,4),9);
                    end
                    BIDLEVEL(brokerind)=num2cell(cat(1,lotchoice{brokerind,7})./(cat(1,Lottype{brokerind,6})+...
                        discount*subplandinfo(cat(1,lotchoice{brokerind,2})).*...
                        cat(1,Lottype{brokerind,3})));
                    sampleinfo=zeros(length(brokerind),6);
                    sampleinfo(:,1)=cat(1,lotchoice{brokerind,7});
                    sampleinfo(:,2)=cat(1,lotchoice{brokerind,3});
                    sampleinfo(:,3)=cat(1,lotchoice{brokerind,4});
                    sampleinfo(:,4)=nbids;
                    sampleinfo(:,5)=cat(1,lotchoice{brokerind,7})./(cat(1,Lottype{brokerind,6})+...
                        discount*subplandinfo(cat(1,lotchoice{brokerind,2})).*cat(1,Lottype{brokerind,3}));
                    sampleinfo(:,6)=cat(1,lotchoice{brokerind,8});
                    
                    subexpthouse=zeros(HT,1);
                    for lt=1:HT
                        ils=(iblots(:,2)==lt);
                        houseinfo(lt,2,ibr,t)=z(lt,1);
                        houseinfo(lt,3,ibr,t)=z(lt,2);
                        if isempty(find(ils,1))==1
                            subexpthouse(lt)=0;
                            houseinfo(lt,[1 4:7],ibr,t)=0;
                        else
                            subexpthouse(lt)=mean(cat(1,lotchoice{iblots(ils,1),7}));
                            houseinfo(lt,1,ibr,t)=mean(sampleinfo(ils,1));
                            houseinfo(lt,4,ibr,t)=mean(sampleinfo(ils,4));
                            houseinfo(lt,5,ibr,t)=min(sampleinfo(ils,5));
                            houseinfo(lt,6,ibr,t)=length(find(ils==1));
                            houseinfo(lt,7,ibr,t)=mean(cat(1,coastprox{cat(1,Lottype{iblots(ils,1),2})}));
                        end
                    end
                end
                bb=bbfull(:,:,ibr);
                subexpthouse=zeros(HT,1);
                for lt=1:HT
                    brokererror(lt,:,ibr,t) = (1-DELTA)*brokererror(lt,:,ibr,t-1)+...
                        DELTA*abs(houseinfo(lt,1,ibr,t)-brokerproj(lt,:,ibr));
                    brokerabserror(lt,:,ibr,t)=houseinfo(lt,1,ibr,t)-brokerproj(lt,:,ibr);
                    [brokerbest,ibrokerbest] = min(brokererror(lt,:,ibr,t),[],2);
                    diffbrokererror(lt,:,ibr,t)=brokererror(lt,:,ibr,t)-brokererror(lt,:,ibr,t-1);
                    brokerbestabsSAVE(ibr,lt,t)=brokerabserror(lt,ibrokerbest,ibr,t);
                    brokerbestdiffSAVE(ibr,lt,t)=diffbrokererror(lt,ibrokerbest,ibr,t);
                    brokerbestSAVE(ibr,lt,t) = brokerbest';
                    ibrokerbestSAVE(ibr,lt,t) = ibrokerbest';
                    brokerprojSAVE(ibr,lt,t) = brokerproj(lt,ibrokerbest,ibr);
                    brokermodelSAVE(ibr,lt,t) = brokermodel(lt,ibrokerbest,ibr);
                    %                     for i = 1:BROKERNUMCLASS
                    %                         strb = sprintf('brokerclass%d = find(brokermodel(lt,:,ibr) == %d);',i,i);
                    %                         eval(strb);
                    %                     end
                    brokerclass1=find(brokermodel(lt,:,ibr)==1);
                    brokerclass2=find(brokermodel(lt,:,ibr)==2);
                    brokerclass3=find(brokermodel(lt,:,ibr)==3);
                    brokerclass4=find(brokermodel(lt,:,ibr)==4);
                    brokerclass5=find(brokermodel(lt,:,ibr)==5);
                    brokerclass6=find(brokermodel(lt,:,ibr)==6);
                    
                    if houseinfo(lt,1,ibr,t) == 0
                        bproj(lt,:)=0;
                    else
                        
                        % mirror models
                        bproj(lt,brokerclass1) = houseinfo(lt,1,ibr,t)+(1-bb...
                            (lt,brokerclass1)).*(0.5*houseinfo(lt,1,ibr,t)-...
                            (houseinfo(lt,1,ibr,t)-houseinfo(lt,1,ibr,t-1)));
                        
                        % mean model
                        for jl = 1:length(brokerclass2)
                            bproj(lt,brokerclass2(jl)) = mean(houseinfo(lt,1,ibr,...
                                t:-1:(t-bb(lt,brokerclass2(jl)))));
                        end
                        
                        %cycle model
                        bproj(lt,brokerclass3) = houseinfo(lt,1,ibr,t-...
                            round(max(1,bb(lt,brokerclass3))));
                        
                        % projection model
                        for jl = 1:length(brokerclass4)
                            %Nonlinear Forecast
                            indata=houseinfo(lt,1,ibr,t-(1+bb(lt,brokerclass4(jl))):t);
                            subindata=reshape(indata,1,length(indata));
                            pcoef=polyfit(1:length(indata),subindata,1);
                            pline=pcoef(1).*(1:length(indata)+1)+pcoef(2);
                            bproj(lt,brokerclass4(jl))=pline(length(pline));
                        end
                        
                        % rescale model
                        bproj(lt,brokerclass5) = bb(lt,brokerclass5)*houseinfo(lt,1,ibr,t);
                        
                        [brows,bcols]=ind2sub([nbrokerlong nbrokerwide],ibr);
                        brnei=(bcols+1)*nbrokerlong-(nbrokerwide-brows);
                        blnei=(bcols-1)*nbrokerlong-(nbrokerwide-brows);
                        bupnei=bcols*nbrokerlong-(nbrokerwide-(brows-1));
                        bdnnei=bcols*nbrokerlong-(nbrokerwide-(brows+1));
                        ibnei=[brnei blnei bupnei bdnnei];
                        realbnei=find(minibmap==brnei | minibmap==blnei | ...
                            minibmap==bupnei | minibmap==bdnnei);
                        bproj(lt,brokerclass6) = bb(lt,brokerclass6)*mean(houseinfo(lt,1,realbnei,t));
                        
                        
                    end
                    subexpthouse(lt)=bproj(lt,ibrokerbest);
                end
                ilotlocate=ismember(lotlocate(:,1),iblots(:,1));
                EXPTHOUSE(lotlocate(ilotlocate,2),t+1)=subexpthouse(...
                    cat(1,lotchoice{lotlocate(ilotlocate,1),3}));
                brokerproj(:,:,ibr)=bproj;
            end
            numlt(:,t)=histc(cat(1,Lottype{:,5}),1:HT);
            
            bcheck(:,1:Nbrokers)=houseinfo(:,1,1:Nbrokers,t);
            ibcheck=(bcheck' ~= 0);
            htexist=ismember(1:HT,cat(1,Lottype{:,5}));
            hset=htset(htexist);
            for lt=1:HT
                ilt=(cat(1,lotchoice{:,3})==lt);
                bidlevel(lt,t)=mean(cat(1,lotchoice{ilt,7})./(cat(1,Lottype{ilt,6})+...
                    discount*subplandinfo(cat(1,lotchoice{ilt,2})).*...
                    cat(1,Lottype{ilt,3})));
                ihtexist=(ibcheck(:,lt)==1);
                if isempty(find(ihtexist,1))==1
                    continue
                else
                    avgbrokervar(lt,t)=mean(var(brokerbestabsSAVE(ihtexist,lt,1:t),0,3));
                    probloss(lt,t)=mean(sum((brokerbestabsSAVE(ihtexist,lt,1:t)>0),3)./...
                        length(1:t));
                    abserror=brokerbestabsSAVE(ihtexist,lt,1:t);
                    abserror=reshape(abserror,length(abserror(:,1,1))*t,1);
                    [mu,sigma]=normfit(abserror);
                    phat(lt,:)=[mu sigma];
                    probeven(lt,t)=cdf('norm',0,phat(lt,1),phat(lt,2));
                    probover(lt,t)=length(find(abserror > 0))/length(abserror);
                    probunder(lt,t)=length(find(abserror < 0))/length(abserror);
                    overvalue(lt,t)=icdf('norm',max(min(probeven(lt,t)+...
                        (1-probeven(lt,t))*probover(lt,t),0.99),0.01),phat(lt,1),phat(lt,2));
                    undervalue(lt,t)=icdf('norm',max(min(probeven(lt,t)*...
                        (1-probunder(lt,t)),0.99),0.01),phat(lt,1),phat(lt,2));
                    maxvalue(lt,t)=icdf('norm',probeven(lt,t)+...
                        (1-probeven(lt,t))*0.99,phat(lt,1),phat(lt,2));
                    minvalue(lt,t)=icdf('norm',probeven(lt,t)*...
                        (1-0.99),phat(lt,1),phat(lt,2));
                end
            end
            ihtnexist=(htexist==0);
            isimvar=hset(ismember(htset(htexist),min(simlotrange(htset(ihtnexist),1)):...
                max(simlotrange(htset(ihtnexist),2))));
            avgbrokervar(ihtnexist,t)=max(avgbrokervar(isimvar,t));
            probloss(ihtnexist,t)=alpha_gain/(alpha_gain+alpha_loss);
            overvalue(ihtnexist,t)=mean(maxvalue(isimvar,t));
            undervalue(ihtnexist,t)=mean(minvalue(isimvar,t));
            probover(ihtnexist,t)=alpha_gain/(alpha_gain+alpha_loss);
            probunder(ihtnexist,t)=alpha_loss/(alpha_gain+alpha_loss);
            
            %%% RESULTS %%%
            numtotbids(:,t)=sum(houseinfo(:,4,:,t),3);
            for lt=1:HT
                iocc=(cat(1,lotchoice{:,4})~=0 & cat(1,lotchoice{:,3})==lt);
                htincome(lt,t)=mean(cat(1,CONINFO{cat(1,lotchoice{iocc,5}),1}));
            end
            
            consumerstats(1,t)=length(CONINFO(:,1));
            consumerstats(4,t)=mean(housemp);
            consumerstats(2,t)=mean(cat(1,lotchoice{ifilled,7}));
            consumerstats(3,t)=mean(cat(1,AVGUTIL{ifilled}));
            agrland(t)=length(find(BASELAYER == 0 & reshape(SCAPE,NCELLS,1) == 1));
            
            BIDLEVELMAP(lotlocate(:,2),t)=cat(1,BIDLEVEL{lotlocate(:,1)});
            AVGRENT(lotlocate(:,2),t)=cat(1,lotchoice{lotlocate(:,1),7});
            LOTTYPE(lotlocate(:,2),t)=cat(1,Lottype{lotlocate(:,1),5});
            BASEMAP(lotlocate(:,2))=BASELAYER(lotlocate(:,1));
            
            for ires=1:length(ifilled)
                c=lotchoice{ifilled(ires),5};
                hopt=((CONINFO{c,1}-cat(1,travelcost{cat(1,lotchoice{ifilled(ires),2})})-...
                    avgrent(:,t)).^CONINFO{c,3}).*(cat(1,Lottype{ifilled(ires),4}).^...
                    CONINFO{c,4}).*(cat(1,Lottype{ifilled(ires),3}).^CONINFO{c,5}).*...
                    (cat(1,Lottype{ifilled(ires),7}).^CONINFO{c,6});
                
                profopt=(avgrent(:,t)-ones(HT,1)*subplandinfo(cat(1,lotchoice{ifilled(ires),2}))-...
                    cat(1,Lottype{ifilled(ires),6}))./cat(1,Lottype{ifilled(ires),3});
                [imaxp,jmaxp]=max(profopt,[],1);
                profset(jmaxp,t)=profset(jmaxp,1)+1;
                
                [imaxu,jmaxu]=max(hopt,[],1);
                idealset(jmaxu,t)=idealset(jmaxu,1)+1;
                
            end
            ifill=ismember(lotlocate(:,1),ifilled);
            INCOME(lotlocate(ifill,2),t)=cat(1,CONINFO{cat(1,lotchoice{lotlocate(ifill,1),5}),1});
            PREFMAP(lotlocate(ifill,2),t)=cat(1,CONINFO{cat(1,lotchoice{lotlocate(ifill,1),5}),6});
            SUBRISKMAP(lotlocate(ifill,2),t)=cat(1,CONINFO{cat(1,lotchoice{lotlocate(ifill,1),5}),7});
            
            Rpop(t)=length(ifilled);
            Rvacrate(t)=vacrate(t);
            Rvaclots(t)=vacantlots(t);
            Rleftoverpop(t)=leftoverpop(t);
            newbt=cat(1,Lottype{cat(1,Lottype{:,9})==t,5});
            if isempty(find(newbt,1)) == 1
                htperyear(:,t)=zeros(HT,1);
            else
                htperyear(:,t)=histc(newbt,1:HT);
            end
            
            subfarminfo=LANDINFO{1,t};
            for nf=1:length(iNfarmers)
                Farmdist2dev(iNfarmers(nf),t)=mean(cat(1,Sdist.distmat{subfarminfo==iNfarmers(nf)}));
            end
            subbmodel=zeros(NLENGTH,NWIDTH);
            subbproj=zeros(NLENGTH,NWIDTH);
            for nb=1:Nbrokers
                for lt=1:HT
                    totbrokerrecord(:,lt,t,nb)=[houseinfo(lt,1,nb,t); brokerprojSAVE(nb,lt,t); ...
                        brokermodelSAVE(nb,lt,t)];
                    ibarea=find(HBROKER==nb);
                    idynlt=cat(1,Lottype{cat(1,lotchoice{:,3})==lt,2});
                    itarget=ismember(idynlt,ibarea);
                    subbmodel(idynlt(itarget))=brokermodelSAVE(nb,lt,t);
                    subbproj(idynlt(itarget))=brokerprojSAVE(nb,lt,t);
                    
                end
            end
            Bmodelmap(:,t)=reshape(subbmodel,NCELLS,1);
            Bprojmap(:,t)=reshape(subbproj,NCELLS,1);
        end
        %%% End of time loop
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%    Extract Results Data    %%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % M=M(TSTART:t);
        % MG=MG(TSTART:t);
        % MPL=MPL(TSTART:t);
        % MLT=MLT(TSTART:t);
        % movie(gcf,M,2,1)
        % movie(gcf,MG,2,1)
        % movie(gcf,MPL,2,1)
        % movie(gcf,MLT,2,1)
        
        farmsoldinfo=zeros([],8);
        ifarmsold=find(sellrecord ~= 0);
        avgfarmdist=zeros(length(ifarmsold),1);
        sellepsilon=zeros(length(ifarmsold),1);
        devsellwtp=zeros(length(ifarmsold),1);
        farmsize=zeros(Nfarmers,1);
        for nf=1:Nfarmers
            farmsize(nf)=length(Farminfo{nf,2});
        end
        avgfarmsize=mean(farmsize);
        stdfarmsize=std(farmsize);
        for n=1:length(ifarmsold)
            avgfarmdist(n,1)=mean(Sdist2cbd.dist2cbd(LANDINFO{1,TSTART}==ifarmsold(n)));
            sellepsilon(n,1)=epsilon(sellrecord(ifarmsold(n)));
            devsellwtp(n,1)=wtpland(ifarmsold(n),sellrecord(ifarmsold(n)));
        end
        if isempty(find(ifarmsold,1))==1
            disp('No farm sales')
        else
            pctfarmgain=(buyrecord(ifarmsold)-wtaland(ifarmsold,1))./wtaland(ifarmsold,1);
            farmsoldinfo=[ifarmsold sellrecord(ifarmsold) wtaland(ifarmsold,1) buyrecord(ifarmsold) sellepsilon devsellwtp avgfarmdist pctfarmgain];
            sortfarmsold=sortrows(farmsoldinfo,1);
        end
        Allfarmdist=zeros(Nfarmers,1);
        WTPlandmap=zeros(NCELLS,TMAX);
        WTAlandmap=zeros(NCELLS,TMAX);
        Landprojmap=zeros(NCELLS,TMAX);
        Landmodmap=zeros(NCELLS,TMAX);
        totfarmrecord=mat2cell([wtaland(:,1:TMAX); wtpland; landprojSAVE; landmodelSAVE],ones(4,1)*Nfarmers,TMAX);
        subwtaland=totfarmrecord{1};
        subwtpland=totfarmrecord{2};
        sublandproj=totfarmrecord{3};
        sublandmodel=totfarmrecord{4};
        for ts=11:TMAX
            farmind=unique(LANDINFO{1,ts});
            farmind=farmind(farmind~=0);
            for nf=1:length(farmind)
                WTAlandmap(ismember(LANDINFO{1,ts},farmind(nf)),ts)=subwtaland(farmind(nf),ts);
                WTPlandmap(ismember(LANDINFO{1,ts},farmind),ts)=subwtpland(farmind(nf),ts);
                Landprojmap(ismember(LANDINFO{1,ts},farmind),ts)=sublandproj(farmind(nf),ts);
                Landmodmap(ismember(LANDINFO{1,ts},farmind),ts)=sublandmodel(farmind(nf),ts);
            end
        end
        
        %%%% tolerence levels %%%%
        maxvacrate=0.15;
        maxvacpop=0.20;
        maxpopchange=0.8;
        minhousetypes=3;
        maxstd=0.5;
        maxlotchange=0.20;
        maxdevchange=0.20;
        maxpctriserent=1;
        maxratediff=0.25;
        maxutildiff=0.85;
        
        rentdynstats=zeros(HT,TMAX-TSTART);
        
        avgrentdynms=diff(avgrent(:,TSTART:TMAX),1,2);
        pctbiglots=length(find(LOTTYPE > 9))/length(find(BASELAYER==1));
        
        for lt=1:HT
            rentdynstats(lt,:)=avgrentdynms(lt,:)./avgrent(lt,TSTART:TMAX-1);
        end
        
        % Mean dispersion (Irwin and Bockstael, 2007)
        altprop=zeros(length(iurblist),1);
        subbaselayer=reshape(BASELAYER,NLENGTH,NWIDTH);
        for ic=1:length(iurblist)
            [row,col]=ind2sub([NLENGTH NWIDTH],iurblist(ic));
            updir=max(row-2,1);
            dndir=min(row+2,NLENGTH);
            lfdir=max(col-2,1);
            rtdir=min(col+2,NWIDTH);
            altprop(ic)=length(find(subbaselayer(updir:dndir,lfdir:rtdir)==0))/...
                (length(updir:dndir)*length(lfdir:rtdir));
        end
        meandisp=sum(altprop)/length(iurblist);
        maxdevdist=max(Sdist2cbd.dist2cbd(iurblist));
        
        distzones=ceil(max(max(Sdist2cbd.dist2cbd))/13);
        zonedensity=zeros(distzones,1);
        diststart=0;
        distlim=Sdist2cbd.dist2cbd(indevedge);
        
        for dz=1:distzones
            idistzone=find(Sdist2cbd.dist2cbd > diststart & Sdist2cbd.dist2cbd <= distlim);
            zonedensity(dz)=mean(z(LOTTYPE(idistzone(ismember(idistzone,iurblist))),1));
            diststart=diststart+indevedge;
            distlim=distlim+indevedge;
        end
        
        bidshare=zeros(HT,30);
        buildshare=zeros(HT,30);
        for t=11:30
            subshare=numtotbids(:,t)./htperyear(:,t);
            ibidshare=(isinf(subshare) == 0 & isnan(subshare)==0);
            bidshare(ibidshare,t)=subshare(ibidshare)/sum(subshare(ibidshare),1);
            buildshare(:,t)=htperyear(:,t)/sum(htperyear(:,t),1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %         VARLAYER=VARLAYER+BASELAYER;    %record frequency of development per cell across model runs
        
       
%         ndate=datestr(date,'ddmmyy');
%         savefname=sprintf('storm_%d_%d.mat',erun,mrun);
%         parsave_storm(savefname,...
%             consumerstats,vacstats,BUILDTIME,VACLAND,RENT,RETURN,LOTTYPE,...
%             BASELAYER,Rpop,Rvacrate,Rvaclots,numlt,Rleftoverpop,avgrentdynms,...
%             rentdynstats,farmsoldinfo,avgrent,avgfarmsize,stdfarmsize,DELTA,...
%             survivalrate,LOCWGHT,REGWGHT,PCTSEARCH,zeta,HIBETA,MIDBETA,LOWBETA,...
%             POPGROW,ccost,newhouseset,bidtot,meandisp,maxdevdist,setupmap,....
%             zonedensity,vacrate,Farminfo,oldincome,Realreturn,...
%             Realavgret,Exptrentdiff, Avgexptdiff,htincome,numtotbids,...
%             totfarmrecord,htperyear,Newbidlevel,totbrokerrecord,Farmdist2dev,...
%             Bprojmap,Bmodelmap,WTAlandmap,WTPlandmap,Landprojmap,Landmodmap,...
%             bidshare,buildshare,landdemand,EXPTHOUSE,lotchoice,Exptprofit,...
%             Exptret,Realexptret,Realavgexptret,idealset,profset,avgbrokervar,...
%             carrycost,Lottype,CONINFO,PREFMAP,TSI,IMPACT,DAMAGE,LANDINFO,lotlocate);
        cd X:\model_results\CHALMS_alt_storm_climate

                fname={'results_CHALMS_Coast_singlerun_3'};
        
                save(fname{:},'consumerstats','vacstats','BUILDTIME','VACLAND','RENT','RETURN',...
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
    end
end
toc