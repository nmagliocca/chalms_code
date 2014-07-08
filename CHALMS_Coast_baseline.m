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
income1=(income >= himinwage & income <= himaxwage);
% Different proportions of income spent on housing depending on income
% level
housepref(income1)=HIBETA(1)+(HIBETA(2)-HIBETA(1))*rand(length(find(income1==1)),1);
housepref(income2)=MIDBETA(1)+(MIDBETA(2)-MIDBETA(1))*rand(length(find(income2==1)),1);
housepref(income3)=LOWBETA(1)+(LOWBETA(2)-LOWBETA(1))*rand(length(find(income3==1)),1);

CONINFO(:,3)=num2cell(1-housepref);
CONINFO(:,6)=num2cell((ampref_min+(ampref_max-ampref_min)*rand(length(housepref),1)).*housepref);
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
stream.Substream=mrun;
stream.State=savedState;

%%% Spin-up housing market, developer learns pricing
cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\data_files
load rentmodel
cd C:\Users\nmagliocca\Documents\Matlab_code\CHALMS_coast\base-chalms-code

%%% Initialize based on calibrated model
Paskhouse(1:Nlots(1),1)=predict(rentmdl,[cat(1,Lottype{:,3}) cat(1,Lottype{:,8})...
        cat(1,Lottype{:,7}) cat(1,[CONINFO{1:length(CONINFO),1} ...
        CONINFO{1:Nlots(1)-length(CONINFO),1}])']);
%%% Initialize based on construction and travel costs
% Paskhouse(1:Nlots(1),1)=cat(1,Lottype{:,6})-cat(1,travelcost{cat(1,lotchoice{:,2})});

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
    
    HouseMarketInitial_coast_baseline
    
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

BUILDTIME(lotlocate(:,2))=cat(1,Lottype{lotlocate(:,1),9});
BIDLEVELMAP(lotlocate(:,2),1:TSTART)=repmat(cat(1,BIDLEVEL{lotlocate(:,1)}),1,TSTART);
VACLAND(cat(1,vacland{1,1}),1:TSTART)=1;
AVGRENT(lotlocate(:,2),1:TSTART)=repmat(cat(1,lotchoice{lotlocate(:,1),7}),1,TSTART);
LOTTYPE(lotlocate(:,2),1:TSTART)=repmat(cat(1,Lottype{lotlocate(:,1),5}),1,TSTART);
BASEMAP(lotlocate(:,2))=BASELAYER(lotlocate(:,1));
INCOME(lotlocate(:,2),1:TSTART)=repmat(cat(1,CONINFO{cat(1,lotchoice{lotlocate(:,1),5}),1}),1,TSTART);

Rpop(1:TSTART)=length(ifilled);
Rvacrate(1:TSTART)=vacrate(TSTART);
Rvaclots(1:TSTART)=vacantlots(TSTART);
Rleftoverpop(1:TSTART)=leftoverpop(TSTART);
setupmap=AGLAYER;
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

    t
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
        Dynplandmap(:,:,t)=reshape(cat(2,LANDINFO{3,t}),NLENGTH,NWIDTH);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%    Storm Simulator    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Random Storm Generation
    stormfreq=0.1;
    stormint=1/stormfreq;
    stormgen=TSTART:1/stormfreq:TMAX;
    if isempty(find(stormgen==t,1))==0
        probsurf=cat(1,Pflood{:});
        impactprob=1-rand(1);
        iimpact=find(probsurf > impactprob);
        inoimpact=find(probsurf <= impactprob);
        lothit=ismember(lotlocate(:,2),iimpact);
        ihitlist=unique(lotlocate(lothit,1));
        lotmiss=ismember(lotlocate(:,2),inoimpact);
        IMPACT(iimpact,t)=num2cell(ones(length(iimpact),1));
        IMPACT(inoimpact,t)=num2cell(ones(length(inoimpact),1));
        DAMAGE(cat(1,lotchoice{ihitlist,2}))=num2cell(cat(1,DAMAGE{cat(1,lotchoice{ihitlist,2})})+...
            Cdam.*cat(1,lotchoice{ihitlist,7}));
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

    Ndevelopers=1;

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
    for i = 1:POPNUMCLASS
        if i == 1
            % mirror models
            nproj(indclass1) = POP(t-1) - POP(t-2) + (1-dd(indclass1))*(0.5*POP(t-1) - (POP(t-1)-POP(t-2)));
        elseif i == 2
            % mean mode
            for j = 1:length(indclass2)
                nproj(indclass2(j)) = mean(POP((t-1):-1:(t-dd(indclass2(j)))));
            end
        elseif i == 3
            %cycle model
            nproj(indclass3) = POP(t-max(1,dd(indclass3)));
        elseif i == 4
            % projection model
            for j = 1:length(indclass4)
                yyy(j,1:dd(indclass4(j))) = POP((t-dd(indclass4(j))):(t-1));
            end
            sumy = sum(yyy,2);
            sumxy = sum(xxx.*yyy,2);
            slopes = (dd(indclass4)'.*sumxy - sumx.*sumy)./(dd(indclass4)'.*sumx2 - sumx.*sumx);
            intercepts = (sumy - slopes.*sumx)./dd(indclass4)';
            nproj(indclass4) = slopes + intercepts;
        elseif i == 5
            % rescale model
            nproj(indclass5) = dd(indclass5)*POP(t-1);
        end
    end
    errorsq = (1-DELTA)*errorsq + DELTA*(POP(t)-nproj).^2;
    [best(1:Ndevelopers) ibest(1:Ndevelopers)] = min(errorsq,[],2);
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
            regionaldist(lt,t)=mean(dist2cbd(cat(1,Lottype{ireglot,2})));
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
        inddist2dev=distmat{potentialbuy(tl)};

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
                            (dist2cbd(vacrow,vaccol)-regionaldist(lt,t));
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
%                             (dist2cbd(vacrow,vaccol)-regionaldist(lt,t));
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
    subfarminfo=cat(2,LANDINFO{1,t});

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
    
    LandMarket_coast_base
    
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
        [landbest ilandbest] = min(landerror(iNfarmers(nf),:),[],2);
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
  
    %<><><><><><><><><><> Go to housing market <><><><><><><><><><><><><><><><>
    
    HouseMarketDynamic_coast_baseline

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
            for i = 1:BROKERNUMCLASS
                strb = sprintf('brokerclass%d = find(brokermodel(lt,:,ibr) == %d);',i,i);
                eval(strb);
            end
            if houseinfo(lt,1,ibr,t) == 0
                bproj(lt,:)=0;
            else
                for i = 1:BROKERNUMCLASS
                    if i == 1
                        % mirror models
                        bproj(lt,brokerclass1) = houseinfo(lt,1,ibr,t)+(1-bb...
                            (lt,brokerclass1)).*(0.5*houseinfo(lt,1,ibr,t)-...
                            (houseinfo(lt,1,ibr,t)-houseinfo(lt,1,ibr,t-1)));
                    elseif i == 2
                        % mean model
                        for jl = 1:length(brokerclass2)
                            bproj(lt,brokerclass2(jl)) = mean(houseinfo(lt,1,ibr,...
                                t:-1:(t-bb(lt,brokerclass2(jl)))));
                        end
                    elseif i == 3
                        %cycle model
                        bproj(lt,brokerclass3) = houseinfo(lt,1,ibr,t-...
                            round(max(1,bb(lt,brokerclass3))));
                    elseif i == 4
                        % projection model
                        for jl = 1:length(brokerclass4)
                            %Nonlinear Forecast
                            indata=houseinfo(lt,1,ibr,t-(1+bb(lt,brokerclass4(jl))):t);
                            subindata=reshape(indata,1,length(indata));
                            pcoef=polyfit(1:length(indata),subindata,1);
                            pline=pcoef(1).*(1:length(indata)+1)+pcoef(2);
                            bproj(lt,brokerclass4(jl))=pline(length(pline));
                        end
                    elseif i == 5
                        % rescale model
                        bproj(lt,brokerclass5) = bb(lt,brokerclass5)*houseinfo(lt,1,ibr,t);
                    elseif i == 6
                        [brows bcols]=ind2sub([nbrokerlong nbrokerwide],ibr);
                        brnei=(bcols+1)*nbrokerlong-(nbrokerwide-brows);
                        blnei=(bcols-1)*nbrokerlong-(nbrokerwide-brows);
                        bupnei=bcols*nbrokerlong-(nbrokerwide-(brows-1));
                        bdnnei=bcols*nbrokerlong-(nbrokerwide-(brows+1));
                        ibnei=[brnei blnei bupnei bdnnei];
                        realbnei=find(minibmap==brnei | minibmap==blnei | ...
                            minibmap==bupnei | minibmap==bdnnei);
                        bproj(lt,brokerclass6) = bb(lt,brokerclass6)*mean(houseinfo(lt,1,realbnei,t));
                    end
                end
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
        Farmdist2dev(iNfarmers(nf),t)=mean(cat(1,distmat{subfarminfo==iNfarmers(nf)}));
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
    %%
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
    avgfarmdist(n,1)=mean(dist2cbd(LANDINFO{1,TSTART}==ifarmsold(n)));
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
maxdevdist=max(dist2cbd(iurblist));

distzones=ceil(max(max(dist2cbd))/13);
zonedensity=zeros(distzones,1);
diststart=0;
distlim=dist2cbd(indevedge);

for dz=1:distzones
    idistzone=find(dist2cbd > diststart & dist2cbd <= distlim);
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
