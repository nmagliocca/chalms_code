%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%    MASTER CODE    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
tic
stream=RandStream.create('mrg32k3a','Seed',13);
RandStream.setGlobalStream(stream);
stream.Substream=86;
repeatstate1=stream.State;
stream.Substream=35;
repeatstate2=stream.State;
stream.Substream=4;
repeatstate3=stream.State;

MRUNS=30;
rstate=[1 5 18 69 74 20 49 11 56 1009 23 47 58 13 85 52 29 85 46 6 99 216 549 876 316 545 468 736 984 2546];
parmfit=zeros(MRUNS,7);

VARLAYER=zeros(80*80,1);
%%
for mrun=1:MRUNS
    %%
    cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\base-chalms-code
    stream.Substream=mrun;

    disp(mrun)

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
    
    EmpDataInput_coast_base
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%    PARAMETERS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    DELTA=0.0527;
    survivalrate=0.9499;
%     LOCWGHT=0.5;
    LOCWGHT=0.2955;
    REGWGHT=1-LOCWGHT;
    PCTSEARCH=0.5887;
    zeta=0.5;   %dampening on epsilon-- must be >0, controls the rate of price ...
    %increase/preception of market power

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
    celldiag=(cellside*sqrt(2));          %miles diag to neighboring cell's center
    acre2sqft=43560;
    avgdist2nei=mean([cellside celldiag]);
    margtc=1.30*500*cell2mile;        %Assumed: 250 travel days a year, roundtrip
    
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
    
    stormthresh=15;
    Cdam=0.5;
    Cmit=3000;
    miteff=1;
    
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
     
    %%% Zones Layer - for use with zoning settings %%%
    ZONES=cell(NZONES,3);   %[cellid,min zone,max zone]
    Zones=zeros(NLENGTH,NWIDTH);
    zoning=zeros(NZONES,2);     %[(min lotsize) (max lotsize)]
    ZONEMAP=zeros(NLENGTH,NWIDTH);
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
    landerror = zeros(Nfarmers,NUMMODEL);
    landbestSAVE=zeros(Nfarmers,TMAX);
    ilandbestSAVE=zeros(Nfarmers,TMAX);
    landprojSAVE=zeros(Nfarmers,TMAX);
    landprojbestSAVE=zeros(Nfarmers,NUMMODEL);
    landmodelbestSAVE=zeros(Nfarmers,TMAX);
    landmodelSAVE=zeros(Nfarmers,TMAX);
    
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
    FARMPROJ=zeros(NCELLS,1);
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
    AVGUTIL=cell([],1);
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
    
    %%%% RESULTS %%%%
    devcells=zeros(2,TMAX);         %[#dev_cells; pct_dev_land]
    vacstats=zeros(HT,TMAX);
    avgrent=zeros(HT,TMAX);
    vacrate=zeros(1,TMAX);
    nohouse=zeros(1,TMAX);
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
    setupmap=zeros(NLENGTH,NWIDTH);
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
    
    totfarmrecord=zeros(4,TMAX,Nfarmers);   %[wtaland; wtpland; landprojSAVE; landmodelSAVE]
    Allfarmdist=zeros(Nfarmers,1);
    Farmdist2dev=zeros(Nfarmers,TMAX);
    Dynltmap=zeros(NCELLS,TMAX);
    Dynrentmap=zeros(NCELLS,TMAX);
    Dynmaxretmap=zeros(NCELLS,TMAX);
    Dynretltmap=zeros(NCELLS,TMAX);
    Dynmaxprofmap=zeros(NCELLS,TMAX);
    Dynprofltmap=zeros(NCELLS,TMAX);
    Dyneultmap=zeros(NCELLS,TMAX);
    ZONED=zeros(NCELLS,1);
    
    subdynlt=zeros(NCELLS,1);
    subdynrent=zeros(NCELLS,1);
    Bmodelmap=zeros(NCELLS,TMAX);
    Bprojmap=zeros(NCELLS,TMAX);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%    Load Reference Landscape    %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    savedState=stream.State;

    RefLandscape_Coast_baseline

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%    Model Runs    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    
    CHALMS_Coast_baseline
    
    VARLAYER=VARLAYER+BASELAYER;    %record frequency of development per cell across model runs
    
    cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\results
    
    ndate=datestr(date,'ddmmyy');
    fname=sprintf('results_CHALMS_Coast_baseline_%d',mrun);
    
    save(fname,'consumerstats','vacstats','BUILDTIME','VACLAND','RENT','RETURN',...
        'LOTTYPE','BASELAYER','Rpop','Rvacrate','Rvaclots',...
        'numlt','Rleftoverpop','avgrentdynms','rentdynstats',...
        'farmsoldinfo','avgrent','avgfarmsize','stdfarmsize','DELTA',...
        'survivalrate','LOCWGHT','REGWGHT','PCTSEARCH','zeta','HIBETA',...
        'MIDBETA','LOWBETA','POPGROW','ccost','newhouseset','bidtot',...
        'meandisp','maxdevdist','setupmap','zonedensity',...
        'VARLAYER','vacrate','Farminfo','oldincome','Realreturn',...
        'Realavgret','Exptrentdiff','Avgexptdiff','htincome',...
        'numtotbids','totfarmrecord','htperyear','Newbidlevel',...
        'totbrokerrecord','Farmdist2dev','Bprojmap','Bmodelmap',...
        'WTAlandmap','WTPlandmap','Landprojmap','Landmodmap',...
        'bidshare','buildshare','landdemand','EXPTHOUSE','lotchoice',...
        'Exptprofit','Exptret','Realexptret','Realavgexptret','idealset',...
        'profset','avgbrokervar','carrycost','Lottype','CONINFO','PREFMAP',...
        'TSI','IMPACT','DAMAGE','LANDINFO','lotlocate')
end
toc