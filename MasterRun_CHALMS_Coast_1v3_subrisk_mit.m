%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%    MASTER CODE    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% s=RandStream.create('mt19937ar','seed',13);
% RandStream.setDefaultStream(s);
stream=RandStream.create('mrg32k3a','Seed',13);
RandStream.setGlobalStream(stream);
stream.Substream=86;
repeatstate1=stream.State;
stream.Substream=35;
repeatstate2=stream.State;
stream.Substream=4;
repeatstate3=stream.State;

MRUNS=1;
rstate=[1 5 18 69 74 20 49 11 56 1009 23 47 58 13 85 52 29 85 46 6 99 216 549 876 316 545 468 736 984 2546];
% s=RandStream.create('mt19937ar','seed',sum(100*clock));
% RandStream.setDefaultStream(s);
parmfit=zeros(MRUNS,7);

VARLAYER=zeros(80,80);
%%
for mrun=1:MRUNS
    %%
    stream.Substream=mrun;

    
%     reset(s,rstate(mrun));
%     s=RandStream.create('mt19937ar','seed',rstream(mrun));
%     RandStream.setDefaultStream(s);
    mrun
    % load globalvars_103109
    % load masterprange
    % load parmsample_101009
    
%     rand('state',86);
%     s=rand('state');
%     randn('state',86)
%     sn=randn('state');
    
    % rand('state',rndstate(ia))
    % randn('state',rndstate(ia))
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    % M-files needed to run
    % 1. EmpDataInput.m
    % 2. LandMarket2.m
    % 3. HouseMarketInitial4.m
    % 4. HousemarketDynamic3.m
    % 5. ABHLM_V2.m
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    %%% Housing Layer %%%
    HTYPE=3;
    LTYPE=3;
%     AMLEV=4;
    HT=HTYPE*LTYPE;
    
    % global NZONES NCELLS NLENGTH NWIDTH TMAX travelcost distftop LOTSIZE
    
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@    INITIAL CONDITIONS    @@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    EmpDataInput_2v2_dev
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%    PARAMETERS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % DELTA=parmpop(nco,1);
    % survivalrate=parmpop(nco,2);
    % LOCWGHT=parmpop(nco,3);
    % REGWGHT=parmpop(nco,4);
    % PCTSEARCH=parmpop(nco,5);
    % zeta=parmpop(nco,6);
    
    % DELTA=masterparms(2,1);
    % survivalrate=masterparms(2,2);
    % LOCWGHT=masterparms(2,3);
    % REGWGHT=1-LOCWGHT;
    % PCTSEARCH=masterparms(2,5);
    % zeta=masterparms(2,6);
    
%     DELTA=0.2382;
%     survivalrate=0.5958;
%     LOCWGHT=0.4654;
%     REGWGHT=1-LOCWGHT;
%     PCTSEARCH=0.8538;
%     zeta=0.1641;

    DELTA=0.0527;
    survivalrate=0.9499;
    LOCWGHT=0.2955;
    REGWGHT=1-LOCWGHT;
    PCTSEARCH=0.5887;
%     zeta=0.1655;
    zeta=0.5;
    
    
    % DELTA=0.1;
    % survivalrate=0.5738;
    % LOCWGHT=0.6;
    % REGWGHT=1-LOCWGHT;
    % PCTSEARCH=0.5486;
    % zeta=0.3051;

    TMAX=30;
    TSTART=10;
    
    %%% Landscape Layer %%%
    NWIDTH=80;
    NLENGTH=80;
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
    Nfarmers=50;
    
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
    % PCTSEARCH=0.4;
    % LOCWGHT=0.0;
    % REGWGHT=1-LOCWGHT;
    
    %%% Consumers %%%
    % Nconsumers=367;
    
    % RESSPAN=7;
    maxwage=max(incomedata);
    % maxwage=200000;
    minwage=min(incomedata);
    % minwage=20000;
    himaxwage=maxwage;
    himinwage=100000;
    midmaxwage=99999;
    midminwage=60000;
    lowmaxwage=59999;
    lowminwage=minwage;
    WAGECLASS=3;
    HIBETA=[0.16 0.24];     %approximate max and mins for housepref come from calvert county census data
    MIDBETA=[0.25 0.33];
    LOWBETA=[0.34 0.42];
    incomeg=0.05;
    incomesigma=1;
    incomegrow=0.005;
    
    %%% Population %%%
    popg=0.05;
    popsigma=1;
    pop2dem=2.91;   %persons per household, Calvert County quickfacts
    % POPGROW=0.06;   %ditto (19% actually)
    POPGROW=0.10;   %ditto (19% actually)
    popurb=0.7;
    popag=0.3;
    
    % zeta=0.50;   %dampening on epsilon-- must be >0, controls the rate of price ...
    %increase/preception of market power
    
    discount=0.05;
%     discount=0.06;
    tax=0.05;
    thetavac=0.2;
    thetadev=0.2;
    
    %Flags
    epflag=0;
    oldbuildflag=0;
    newbuildflag=0;
    termrunflag=0;
    
    alpha_gain=3;   % skewedness factor (Ligmann-Zielinska, 2009)
    alpha_loss=2.5;
    w_gain=0.5;
    w_loss=0.5;
    proftarget=5000*discount;
    
    Cdam=60000;
    Cmit=3000;
    miteff=0.25;
    
    testtime=[10 15 20 25 30];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%    VARIABLES    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Landscape Layer %%%
    SCAPE=zeros(NLENGTH,NWIDTH);
    BASELAYER=zeros(NLENGTH,NWIDTH);  %developed or undeveloped
    DISTANCE=ones(NLENGTH,NWIDTH);
    dist2hznnei=zeros(NLENGTH,NWIDTH);    %distance to horizontal neighbor from icenter
    dist2vrtnei=zeros(NLENGTH,NWIDTH);
    travelcost=zeros(NLENGTH,NWIDTH);
    
    Pflood=zeros(size(SCAPE));
    TSI=zeros(size(SCAPE));
    IMPACT=zeros(NLENGTH,NWIDTH,TMAX);
    DAMAGE=zeros(NLENGTH,NWIDTH);
    stormstats=zeros(2,TMAX);   %[storm occurance; number of cells impacted;
    MITIGATE=zeros(NLENGTH,NWIDTH);
    
    %%% Zones Layer %%%
    ZONES=zeros(NLENGTH,NWIDTH);
    zoning=zeros(NZONES,2);     %[(min lotsize) (max lotsize)]
    ZONEMAP=zeros(NLENGTH,NWIDTH);
    %%% Housing Layer %%%
    HOUSESIZE=zeros(NLENGTH,NWIDTH);
    AMLEVEL=zeros(NLENGTH,NWIDTH);
    LOTS=zeros(NLENGTH,NWIDTH);
    LOTTYPE=zeros(NLENGTH,NWIDTH);
    
    %%% Broker Layer %%%
    HBROKER=zeros(NLENGTH,NWIDTH);
    UUcell=zeros(NLENGTH,NWIDTH);
    
    %%% Agricultural Layer %%%
    AGLAYER=zeros(NLENGTH,NWIDTH);
    LANDVALUE=zeros(NLENGTH,NWIDTH);
    PLAND=zeros(NLENGTH,NWIDTH);
    
    %%% Farmers %%%
    Farmstats=zeros(Nfarmers,5,TMAX);    %[id acres prod costs value_acre]
    iNfarmers=(1:Nfarmers)';
    wtaland=zeros(Nfarmers,TMAX);
    Paskland=zeros(Nfarmers,TMAX);
    transprice=zeros(Nfarmers,TMAX);
    Plandproj=zeros(Nfarmers,TMAX);
    maxrent=zeros(Nfarmers,TMAX);
    maxreturn=zeros(Nfarmers,TMAX);
    landsold=zeros(NZONES,TMAX);
    sellrecord=zeros(Nfarmers,1);
    pdrrecord=zeros(Nfarmers,1);
    buyrecord=zeros(Nfarmers,1);
    
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
    Paskhouse=zeros([],1);
    Ccost=zeros([],1);
    Pland=zeros([],1);
    Plot=zeros([],1);
    INRENT=zeros(NLENGTH,NWIDTH);
    BUILDTIME=zeros(NLENGTH,NWIDTH);
    PROFIT=zeros(NLENGTH,NWIDTH,TMAX);
    RETURN=zeros(NLENGTH,NWIDTH,HT);
    Uret=zeros(NLENGTH,NWIDTH,HT);
    PCHANGE=zeros(NLENGTH,NWIDTH,HT);
    RETURNPROJ=zeros(NLENGTH,NWIDTH);
    LOTPRICE=zeros(NLENGTH,NWIDTH);
    LOTVALUE=zeros(NLENGTH,NWIDTH);
    INITIALPASK=zeros(NLENGTH,NWIDTH);
    RENT=zeros(NLENGTH,NWIDTH);
    house2cell=zeros(NLENGTH,NWIDTH);
    RENTGRAD=zeros(NLENGTH,NWIDTH);
    XBIDS=zeros(NLENGTH,NWIDTH);
    VACLAND=zeros(NLENGTH,NWIDTH);
    FARMPROJ=zeros(NLENGTH,NWIDTH);
    RENTPROJ=zeros(NLENGTH,NWIDTH);
    RENTPROJLAND=zeros(NLENGTH,NWIDTH);
    maxRENTPROJ=zeros(NLENGTH,NWIDTH);
    CCOST=zeros(NLENGTH,NWIDTH);
    ZZ=zeros(NLENGTH,NWIDTH);
    htset=(1:HT)';
    % subPLAND=zeros(NLENGTH,NWIDTH);
    subRENTPROJ=zeros(NLENGTH,NWIDTH,HT);
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
    WTPMAP=zeros(NLENGTH,NWIDTH,TMAX);
    highRETURN=zeros(NLENGTH,NWIDTH,HT);
    lowRETURN=zeros(NLENGTH,NWIDTH,HT);
    potgain=zeros(NLENGTH,NWIDTH,HT);
    potloss=zeros(NLENGTH,NWIDTH,HT);
    EU_dev=zeros(NLENGTH,NWIDTH,HT);
    MAXRET=zeros(NLENGTH,NWIDTH,TMAX);
%     cvht=zeros(HT,TMAX);
    maxcount=zeros(NLENGTH,NWIDTH);
    
    %%% Consumers %%%
    Income=zeros([],2);     %[Income Search_time]
    OutIncome=zeros([],2);
    damcoef=zeros(NLENGTH,NWIDTH,[],TMAX);
    RESTIME=zeros(NLENGTH,NWIDTH);
    avghousemp=zeros(1,TMAX);
    CONINFO=zeros([],5);    %[consumer_good housesize lotsize proximity subrisk ]
    NCON=zeros(WAGECLASS,TMAX);
    searchtimemin=2;
    searchtimemax=6;
    
    %%% Population %%%
    POP=zeros(1,TMAX);
    
    %%% Brokers %%%
    AVGUTIL=zeros(NLENGTH,NWIDTH);
    AVGINCOME=zeros(NLENGTH,NWIDTH);
    avg_income=zeros(NLENGTH,NWIDTH);
    avg_alpha=zeros(NLENGTH,NWIDTH);
    avg_beta=zeros(NLENGTH,NWIDTH);
    avg_gamma=zeros(NLENGTH,NWIDTH);
    avg_ampref=zeros(NLENGTH,NWIDTH);
    avgstats=zeros(Nbrokers,5);%[avg_income avg_alpha avg_beta avg_gamma avg_ampref]
    houseinfo=zeros(HT,7,Nbrokers,TMAX);      %[avg_price lotsize housesize #of_bidders %above_Pask #lots(lt) amlevel]
    EXPTHOUSE=zeros(NLENGTH,NWIDTH);
    MINBIDLEVEL=zeros(NLENGTH,NWIDTH);
    BIDLEVEL=zeros(NLENGTH,NWIDTH);
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
    EUrankret=zeros(NLENGTH,NWIDTH,HT);
    EUlandret=zeros(NLENGTH,NWIDTH,HT);
    
    %%%% RESULTS %%%%
    devcells=zeros(2,TMAX);         %[#dev_cells; pct_dev_land]
    vacstats=zeros(HT,TMAX);
    avgrent=zeros(HT,TMAX);
    vacrate=zeros(1,TMAX);
    nohouse=zeros(1,TMAX);
    oldincome=zeros(1,TMAX);
    consumerstats=zeros(5,TMAX);     %[Nconsumers; avg_income avg_util housemp income_out]
    vacantlots=zeros(1,TMAX);
    leftoverpop=zeros(1,TMAX);
    agrland=zeros(1,TMAX);
    BT=zeros(NLENGTH,NWIDTH,length(testtime));
    BL=zeros(NLENGTH,NWIDTH,length(testtime));
    Rvacland=zeros(NLENGTH,NWIDTH,length(testtime));
    Rrent=zeros(NLENGTH,NWIDTH,length(testtime));
    subRreturn=zeros(NLENGTH,NWIDTH);
    Rreturn=zeros(NLENGTH,NWIDTH,length(testtime));
    Rlottype=zeros(NLENGTH,NWIDTH,length(testtime));
    Rincome=zeros(NLENGTH,NWIDTH,length(testtime));
    Rbaselayer=zeros(NLENGTH,NWIDTH,length(testtime));
    Rpland=zeros(NLENGTH,NWIDTH,length(testtime));
    PREFMAP=zeros(NLENGTH,NWIDTH,length(testtime));
    SUBRISKMAP=zeros(NLENGTH,NWIDTH,length(testtime));
    Rpop=zeros(1,length(testtime));
    Rvacrate=zeros(1,length(testtime));
    Rvaclots=zeros(1,length(testtime));
    Rnumlot=zeros(1,length(testtime));
    Rleftoverpop=zeros(1,length(testtime));
    idealset=zeros(HT,length(testtime));
    profset=zeros(HT,length(testtime));
    setupmap=zeros(NLENGTH,NWIDTH);
    utilgini=zeros(1,TMAX);
    incgini=zeros(1,TMAX);
    Exptrentdiff=zeros(NLENGTH,NWIDTH,TMAX);
    Avgexptdiff=zeros(HT,TMAX);
    Newlottype=zeros(NLENGTH,NWIDTH,TMAX);
    Realreturn=zeros(NLENGTH,NWIDTH,TMAX);
    Realexptret=zeros(NLENGTH,NWIDTH,TMAX);
    Realavgret=zeros(HT,TMAX);
    Realavgexptret=zeros(HT,TMAX);
    Newbidlevel=zeros(NLENGTH,NWIDTH,TMAX);
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
    Dynltmap=zeros(NLENGTH,NWIDTH,TMAX);
    Dynrentmap=zeros(NLENGTH,NWIDTH,TMAX);
    Dynmaxretmap=zeros(NLENGTH,NWIDTH,TMAX);
    Dynretltmap=zeros(NLENGTH,NWIDTH,TMAX);
    Dynmaxprofmap=zeros(NLENGTH,NWIDTH,TMAX);
    Dynprofltmap=zeros(NLENGTH,NWIDTH,TMAX);
    Dyneultmap=zeros(NLENGTH,NWIDTH,TMAX);
    ZONED=zeros(NLENGTH,NWIDTH);
    
    subdynlt=zeros(NLENGTH,NWIDTH);
    subdynrent=zeros(NLENGTH,NWIDTH);
    Bmodelmap=zeros(NLENGTH,NWIDTH,TMAX);
    Bprojmap=zeros(NLENGTH,NWIDTH,TMAX);
    %     VARLAYER=zeros(NLENGTH,NWIDTH);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%    Load Reference Landscape    %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    savedState=stream.State;

    RefLandscape_Coast_1v3_subrisk
%     s=RandStream.create('mt19937ar','seed',ceil(1000000*rand(1)));
%     RandStream.setDefaultStream(s);
    stream.Substream=mrun;
    stream.State=savedState;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%    Model Runs    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    
    CHALMS_Coast_1v3_subrisk_mit
    
    VARLAYER=VARLAYER+BASELAYER;
    
%     parmfit(mrun,1)=isempty(find(vacrate > maxvacrate,1));
%     parmfit(mrun,2)=isempty(find(popchange > maxpopchange,1));
%     parmfit(mrun,3)=(length(find(lotchange > 0)) >= 3);
%     parmfit(mrun,4)=isempty(find(deltadev > maxdevchange,1));
%     parmfit(mrun,5)=isempty(find(rentdynstats > maxpctriserent,1));
%     parmfit(mrun,6)=isempty(find(ratediff < -maxratediff,1) | ...
%         find(ratediff > maxratediff,1));
%     parmfit(mrun,7)=isempty(find(pctutildiff(TSTART:TMAX) < maxutildiff,1));

    
    ndate=datestr(date,'ddmmyy');
    fname=sprintf('results_CHALMS_Coast_1v3_subrisk_mit_%d',mrun);
    
    save(fname,'consumerstats','vacstats','BT','Rvacland','Rrent','Rreturn',...
        'Rlottype','Rincome','Rbaselayer','Rpop','Rpland','Rvacrate','Rvaclots',...
        'Rnumlot','Rleftoverpop','avgrentdynms','rentdynstats','deltalots',...
        'farmsoldinfo','avgrent','avgfarmsize','stdfarmsize','DELTA',...
        'survivalrate','LOCWGHT','REGWGHT','PCTSEARCH','zeta','HIBETA',...
        'MIDBETA','LOWBETA','POPGROW','ccost','newhouseset','bidtot',...
        'meandisp','maxdevdist','utilgini','incgini','setupmap','zonedensity',...
        'pctutildiff','VARLAYER','vacrate','Farmstats','oldincome','Realreturn',...
        'Realavgret','Exptrentdiff','Avgexptdiff','Newlottype','htincome',...
        'numtotbids','totfarmrecord','htperyear','Newbidlevel','Dynltmap',...
        'Dynrentmap','totbrokerrecord','Farmdist2dev','Bprojmap','Bmodelmap',...
        'Dynplandmap','WTAlandmap','WTPlandmap','Landprojmap','Landmodmap',...
        'bidshare','buildshare','farmacreinfo','totfarmrecord','landdemand',...
        'Exptprofit','Exptret','Realexptret','Realavgexptret','idealset',...
        'profset','avgbrokervar','carrycost','con2lot','CONINFO','PREFMAP',...
        'SUBRISKMAP','TSI','IMPACT','DAMAGE','MITIGATE')
end