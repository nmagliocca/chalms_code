 %@@@@@@@@@@@@@@@@@ FULL TDR MODEL @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%    FARMERS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for nf=1:Nfarmers
    str= sprintf('ifarmer = find(AGLAYER == idNfarmers(%d));',nf);
    eval(str);
    Farmstats(nf,2,1:TSTART)=length(find(AGLAYER==idNfarmers(nf)));
    LANDVALUE(ifarmer)=Farmstats(iNfarmers(nf),5,1);
end
PLAND=LANDVALUE;

for it=1:TSTART+1
    Farmstats(:,:,it)=Farmstats(:,:,1);
    wtaland(iNfarmers,it)=Farmstats(iNfarmers,5,1);
    Paskland(iNfarmers,it)=Farmstats(iNfarmers,5,1);
    Plandproj(iNfarmers,it)=Paskland(iNfarmers,it);
end

for lv=1:Nlots(1)
    lvind=Lottype((lv == Lottype(:,1)),2);
    LOTVALUE(lvind)=mean(LANDVALUE(lvind),1);
end

%%%%%%%%%%%% Construction Costs %%%%%%%%%%%%%%%

% ccost_base=zeros(HT,1);
% ccost(1:length(z),1)=(infracost+z(:,2).*sqftcost).*discount;
% % ccost(1:3)=ccost(1:3)-(5000*inflate+5000)*discount;   %lower by streets and sewers, respectively
% % ccost(4:6)=ccost(4:6)-(7000*inflate+8000)*discount;
% ccost(1:3)=ccost(1:3)-(7000*inflate+8000)*discount;

ccost_base=([208273.7676
    270773.7676
    333273.7676
    219872.3592
    282372.3592
    344872.3592
    268120.5986
    360620.5986
    423120.5986])*discount;

% ccost_base=ccost_base.*(1+[0; 0.1; 0.05; 0; 0.1; 0.05; 0.15]);
ccost=ccost_base;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i=1:length(Lottype(:,1))
%     Lottype(i,6)=ccost(z(:,1)==Lottype(i,3) & z(:,2)==Lottype(i,4));
% end
Lottype(:,6)=ccost(Lottype(:,5));

lotchoice(:,6)=Lottype(ilots,6);
for nl=1:Nlots(1)
    isamelot=find(Lottype(:,1)==lotchoice(nl,1));
    Rmin(nl,1)=ccost(lotchoice(nl,5));
    INRENT(Lottype(isamelot,2))=Rmin(nl,1);
    INITIALPASK(Lottype(isamelot,2))=Rmin(nl,1)-mean(travelcost(Lottype(isamelot,2)));
    wtahouse(nl,1)=Rmin(nl,1);
end

%%% initial vacant land
ivac=(VACLAND ~= 0);
RENTPROJLAND(iagr)=1;
BUILDTIME(iurb)=TSTART;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%    Consumers    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CONINFO(1:Nconsumers,1)=alpha;
CONINFO(1:Nconsumers,2)=beta;
CONINFO(1:Nconsumers,3)=gamma;
CONINFO(1:Nconsumers,4)=ampref;
CONINFO(1:Nconsumers,5)=subrisk;


%<><><><><><><><><>  Matching consumers with houses   <><><><><><><><><><>
Paskhouse(1:Nlots(1),1)=INITIALPASK(lotchoice(:,2));

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
vaccheck=Nconsumers;
%%%%%% Learning Period %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while abs(deltafac) > 0.001
    Lottype(:,7)=0;
    lotchoice(:,7)=0;
    initialdiff=meandiffprice;
    initialutildiff=utildiff(TSTART);
    
    HouseMarketInitial_coast_v1_subrisk
    
    if killflag==1
        break
    end
    
    iter=iter+1; 
    deltadiff(deltaiter)=meandiffprice;
    meandiffprice=mean(abs(Paskhouse-con2lot(:,1)));
    diffprice(1:length(Paskhouse),iter)=con2lot(:,1)-Paskhouse;
    diffvac(iter)=length(istillvac);
    
    deltaiter=deltaiter+1;
    
    for lt=1:HT
        if isempty(find(lotchoice(:,5)==lt,1))==1
            continue
        end
        meanprice(lt,iter)=mean(con2lot(lotchoice(:,5)==lt,1));
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
%     if iter == 1
%         Paskhouse=(1-deltafac).*Paskhouse+deltafac.*con2lot(:,1)./(1+discount);
%     elseif iter > 1
    Paskhouse=Paskhouse+abs(meandiffprice/mean(meanprice(meanprice(:,iter)~=0,iter))).*diffprice(:,iter);
%     end

    for n=1:length(ifilled)
        HU(con2lot(ifilled(n),2),1)=U(con2lot(ifilled(n),2),con2lot(ifilled(n),3));
        isamecell=find(Lottype(:,1)==con2lot(ifilled(n),3));
        RENT(Lottype(isamecell,2))=mean(con2lot(ifilled(n),1));
    end
    %%% Adjust rent for vacant houses
    RENT(lotchoice(istillvac,2))=con2lot(istillvac,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MINBIDLEVEL(lotchoice(ifilled,2))=con2lot(ifilled,1)./Paskhouse(ifilled);

for n=1:length(ifilled)
    HU(con2lot(ifilled(n),2),1)=U(con2lot(ifilled(n),2),con2lot(ifilled(n),3));
    isamecell=find(Lottype(:,1)==con2lot(ifilled(n),3));
    RENT(Lottype(isamecell,2))=mean(con2lot(ifilled(n),1));
end
%%% Adjust rent for vacant houses
RENT(lotchoice(istillvac,2))=con2lot(istillvac,1);

% popinhouse=ismember(1:Nconsumers,con2lot(:,2))';
% ioldpop=find(popinhouse==0);
popinhouse=ismember(1:length(Income(:,1)),con2lot(:,2))';

ioldpop=find(popinhouse==0 & Income(:,1)~=0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%    Population    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% POP(1)=175;
NCON(1,1:TSTART)=length(Income(income1,1));
NCON(2,1:TSTART)=length(Income(income2,1));
NCON(3,1:TSTART)=length(Income(income3,1));
POP(1:TSTART)=round(sum(NCON(:,TSTART))./(1+POPGROW).^(TSTART:-1:1));

% for tt=1:TSTART-2
%     POP(tt+1)=round(POP(tt)*(1+POPGROW));
% end
% POP(TSTART)=sum(NCON(:,TSTART));
WCGROWTH=[1/3 1/3 1/3]';
FARMERPOS=AGLAYER;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%    Broker Info    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%@@@@@@@ Broker Projection Models @@@@@@@@@@@@@@@@@@

warning('off','all');
for n=1:Nbrokers
    str= sprintf('find(HBROKER == %d);',n);
    ibroker=eval(str);
    ibcellind=ismember(Lottype(:,2),ibroker);
    ibsample=Lottype(ibcellind,[1:2 5 7]);

    if isempty(ibsample) == 1
        continue
    else
        [iblots,ib,jb]=unique(ibsample(:,1));
        nbids=zeros([],1);
        for b=1:length(iblots)  
            bidover=(Phousebid(:,iblots(b)) > 0);
            nbids(b,1)=length(find(bidover==1));
        end
        
        % Calculate average utilities for rent projection of unknown
        % lottypes
        subrent=con2lot(iblots,1);
        subrent(con2lot(iblots,2)==0)=0;
        for sc=1:length(ibroker)
            isamecell=find(ibsample(:,2)==ibroker(sc));
            isnotvac=(ibsample(isamecell,4) ~= 0);
            if isempty(find(ibsample(isamecell(isnotvac),1),1)) == 1
                continue
            else
%             avg_income(ibroker(sc))=mean(Income(con2lot(ibsample(isamecell(isnotvac),1),2)));
%             avg_alpha(ibroker(sc))=mean(CONINFO(con2lot(ibsample(isamecell(isnotvac),1),2),1));
            avg_income(ibroker(sc))=median(Income(con2lot(ibsample(isamecell(isnotvac),1),2)));
            avg_alpha(ibroker(sc))=median(CONINFO(con2lot(ibsample(isamecell(isnotvac),1),2),1));
            avg_beta(ibroker(sc))=mean(CONINFO(con2lot(ibsample(isamecell(isnotvac),1),2),2));
            avg_gamma(ibroker(sc))=mean(CONINFO(con2lot(ibsample(isamecell(isnotvac),1),2),3));
            avg_ampref(ibroker(sc))=mean(CONINFO(con2lot(ibsample(isamecell(isnotvac),1),2),4));
            AVGUTIL(ibroker(sc))=((avg_income(ibroker(sc))-RENT(ibroker(sc))-...
                travelcost(ibroker(sc)))^avg_alpha(ibroker(sc)))*...
                (mean(z(ibsample(isamecell(isnotvac),3),2))^avg_beta(ibroker(sc)))*...
                (mean(z(ibsample(isamecell(isnotvac),3),1))^avg_gamma(ibroker(sc)))*...
                (mean(coastprox(ibroker(sc)))^avg_ampref(ibroker(sc)));
            end
        end
        
        sampleinfo=zeros(length(iblots),6);
        sampleinfo(:,1)=subrent;
        sampleinfo(:,2)=lotchoice(iblots,3);
        sampleinfo(:,3)=lotchoice(iblots,4);
        sampleinfo(:,4)=nbids;
        sampleinfo(:,5)=con2lot(iblots,1)./(ccost(lotchoice(iblots,5))+...
            discount*PLAND(lotchoice(iblots,2)).*z(lotchoice(iblots,5),1));
        sampleinfo(:,6)=lotchoice(iblots,8);
        for lt=1:HT
            ils=(ibsample(ib,3)==lt);                        
            houseinfo(lt,2,n,1)=z(lt,1);
            houseinfo(lt,3,n,1)=z(lt,2);
            
            if isempty(find(ils,1))==1
                houseinfo(lt,[1 4:6],n,1)=0;
            else
                houseinfo(lt,1,n,1)=mean(sampleinfo(ils,1));
                houseinfo(lt,4,n,1)=mean(sampleinfo(ils,4));
                houseinfo(lt,5,n,1)=min(sampleinfo(ils,5));
                houseinfo(lt,6,n,1)=length(find(ils==1));
            end
            isub=(ibsample(:,3) == lt);
            houseinfo(lt,7,n,1)=mean(coastprox(ibsample(isub,2)));
            BIDLEVEL(ibsample(isub,2))=houseinfo(lt,5,n,1);
            EXPTHOUSE(ibsample(isub,2))=houseinfo(lt,1,n,1);
        end    
    end       
end

warning('on','all');
for lt=1:HT
    numlt(lt,1)=length(find(LOTTYPE==lt))/z(lt,1);
    bidlevel(lt,TSTART)=mean(con2lot(lotchoice(:,5)==lt,1)./(ccost(lt)+...
        discount*PLAND(lotchoice(lotchoice(:,5)==lt,2)).*z(lt,1)));
end

for tt=2:TSTART+1
    houseinfo(:,:,:,tt)=houseinfo(:,:,:,1);
    numlt(:,tt)=numlt(:,1);
end

bcheck(:,1:Nbrokers)=houseinfo(:,1,1:Nbrokers,TSTART);
ibcheck=(bcheck' ~= 0);
htexist=ismember(1:HT,LOTTYPE);
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
%         avgover(lt,TSTART)=mean(abserror(abserror >= 0));
%         avgunder(lt,TSTART)=mean(abserror(abserror <= 0));
%         probover(lt,TSTART)=sum(abs(abserror(abserror >= 0)))/sum(sum(abs(abserror)));
%         probunder(lt,TSTART)=sum(abs(abserror(abserror <= 0)))/sum(sum(abs(abserror)));
%         cvht(lt,TSTART)=std(mean(brokerbestabsSAVE(ihtexist,lt,1:TSTART),3))/...
%             mean(con2lot(lotchoice(:,5)==lt,1));

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
% cvht(ihtnexist,TSTART)=max(cvht(isimvar,TSTART));

house2cell(iurb)=1./z(LOTTYPE(iurb),1);
RENTGRAD(iurb)=RENT(iurb).*house2cell(iurb);

FARMERPOS=AGLAYER;
AGLAYER(iurb)=0;
iremain=(AGLAYER ~= 0);
iNfarmers=unique(AGLAYER(iremain));
landproj=mean(Plandproj(:,TSTART))*(1+min(max(randn(size(landproj)),-0.5),0.5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%    RESULTS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for lt=1:HT
    irent=(houseinfo(lt,1,:,10)~=0);
    avgrent(lt,10)=mean(houseinfo(lt,1,irent,10),3);
    ilt=find(lotchoice(:,5)==lt & lotchoice(:,7)==1 & BUILDTIME(lotchoice(:,2))==TSTART);
    budget_lt(lt,TSTART)=sum(con2lot(ilt,1)-PLAND(lotchoice(ilt,2)).*...
        z(lotchoice(ilt,5),1)*discount-lotchoice(ilt,6));
    vac_ccost(lt,TSTART)=sum(discount*(ccost(lt)+...
        PLAND(lotchoice(istillvac(lotchoice(istillvac,5)==lt),2))*z(lt,1)));
%     vac_ccost(lt,TSTART)=sum((discount*(1+cvht(lt,TSTART)))*...
%         (((discount*(1+cvht(lt,TSTART)))/discount)*ccost(lt)+...
%         PLAND(lotchoice(istillvac(lotchoice(istillvac,5)==lt),2))*z(lt,1)));
end
carrycost(TSTART)=sum(vac_ccost(:,TSTART));    
    
BUDGET(TSTART)=sum(budget_lt(:,TSTART));

INMAP=AGLAYER;
INMAP(iurb)=max(max(AGLAYER))+5;

ZONEMAP(iurb)=-5;

oldincome(TSTART)=mean(Income(ioldpop));
devcells(1,TSTART)=length(iurblist);
devcells(2,TSTART)=devcells(1,TSTART)/(NLENGTH*NWIDTH);

vacantlots(TSTART)=length(istillvac);
leftoverpop(TSTART)=length(ioldpop);

vacrate(TSTART)=vacantlots(TSTART)/Nlots(TSTART);
nohouserate(TSTART)=leftoverpop(TSTART)/POP(TSTART);

agrland(TSTART)=length(find(BASELAYER == 0 & SCAPE == 1));

consumerstats(1,TSTART)=Nconsumers;
consumerstats(4,TSTART)=mean(housemp);
consumerstats(2,TSTART)=mean(Income(con2lot(ifilled,2)));
consumerstats(3,TSTART)=mean(AVGUTIL(iurb));

BT(:,:,1)=BUILDTIME;
BL(:,:,1)=BIDLEVEL;
Rvacland(:,:,1)=VACLAND;
Rrent(:,:,1)=RENT;
Rlottype(:,:,1)=LOTTYPE;
Rbaselayer(:,:,1)=BASELAYER;
Rpland(:,:,1)=PLAND;
Rpop(1)=length(ifilled);
Rvacrate(1)=vacrate(TSTART);
Rvaclots(1)=vacantlots(TSTART);
Rnumlot(1)=numlt(TSTART);
Rleftoverpop(1)=leftoverpop(TSTART);
setupmap=AGLAYER;
Newlottype(:,:,TSTART)=LOTTYPE;
Ufinset=zeros(length(ifilled),1);

for ires=1:length(con2lot(ifilled,1))
%     hopt=((Income(con2lot(ifilled(ires),2),1)-travelcost(lotchoice(ifilled(ires),2))-...
%         avgrent(:,TSTART)).^CONINFO(con2lot(ifilled(ires),2),1)).*(z(:,2).^...
%         CONINFO(con2lot(ifilled(ires),2),2)).*(z(:,1).^CONINFO(con2lot(ifilled(ires),2),3)).*...
%         (z(lt,3).^CONINFO(con2lot(ifilled(ires),2),4));
    hopt=((Income(con2lot(ifilled(ires),2),1)-travelcost(lotchoice(ifilled(ires),2))-...
        avgrent(:,TSTART)).^CONINFO(con2lot(ifilled(ires),2),1)).*(z(:,2).^...
        CONINFO(con2lot(ifilled(ires),2),2)).*(z(:,1).^CONINFO(con2lot(ifilled(ires),2),3)).*...
        (coastprox(lotchoice(ifilled(ires),2)).^CONINFO(con2lot(ifilled(ires),2),4));
    [imaxu,jmaxu]=max(hopt,[],1);
    idealset(jmaxu,1)=idealset(jmaxu,1)+1;
end

for i=1:length(ifilled)
    conid=con2lot(ifilled(i),2);
    lotid=con2lot(ifilled(i),3);
    Ufinset(i)=(((Income(conid)-travelcost(lotchoice(ifilled(i),2))-...
        con2lot(ifilled(i),1)).^CONINFO(conid,1)).*...
        (lotchoice(ifilled(i),4).^CONINFO(conid,2)).*(lotchoice(ifilled(i),3).^...
        CONINFO(conid,3)).*(lotchoice(ifilled(i),8).^CONINFO(conid,4)))';
    
    profopt=(avgrent(:,TSTART)-ones(HT,1)*PLAND(lotchoice(ifilled(i),2))-ccost)./z(:,1);
    [imaxp,jmaxp]=max(profopt,[],1);
    profset(jmaxp,1)=profset(jmaxp,1)+1;
end

Ufinset=sort(Ufinset,'descend');
utilgini(TSTART)=(length(Ufinset)+1)/(length(Ufinset)-1)-2/...
    (length(Ufinset)*(length(Ufinset)-1)*mean(Ufinset))*...
    sum((1:length(Ufinset))'.*Ufinset);

Incomeset=sort(Income(con2lot(ifilled,2)),'descend');
incgini(TSTART)=(length(Incomeset)+1)/(length(Incomeset)-1)-2/...
    (length(Incomeset)*(length(Incomeset)-1)*mean(Incomeset))*...
    sum((1:length(Incomeset))'.*Incomeset);

% figure(4)
% surf(RENT);
% axis ij;
% view(0,90);
% title('RENT, t = TSTART')
% colorbar
% M(TSTART)=getframe(gcf);
% 
% figure(5)
% surf(RENTGRAD);
% axis ij;
% view(0,90);
% title('RENT/acre, t=TSTART')
% colorbar
% MG(TSTART)=getframe(gcf);
% keyboard
% figure(6)
% surf(PLAND);
% axis ij;
% view(0,90);
% title('PLAND, t=TSTART')
% colorbar
% MPL(TSTART)=getframe(gcf);
% 
% figure(7)
% surf(LOTTYPE);
% axis ij;
% view(0,90);
% title('Lot Types, t=TSTART')
% set(gca,'clim',[1 18])
% colorbar
% MLT(TSTART)=getframe(gcf);
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@    DYNAMICS    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
incon=find(Income(:,1)~=0);
iconleave=0;
Nconsumers=length(find(Income(:,1) ~= 0));
%%
for t=TSTART+1:TMAX

    t
    ccost=ccost_base;
%     ccost=ccost_base.*((discount*(1+cvht(:,t-1)))./discount);
    lotchoice(:,6)=ccost(lotchoice(:,5));
%     POP(t)=ceil((POP(t-1)-length(iconleave))*(1+POPGROW));
    POP(t)=ceil(Nconsumers*(1+POPGROW));

    %Existing houses back on market <><><><><><><><><><><><><>
    ileave=find(con2lot(:,4)==t);
    lotchoice(ileave,7)=0;
    returncon=(con2lot(ileave,2)~=0);
    Income(con2lot(ileave(returncon),2),2)=t+...
        ceil(searchtimemin+(searchtimemax-searchtimemin)*rand(length(ileave(returncon)),1));
    for lt=1:HT
        leavetype=(lotchoice(ileave,5)==lt);
        newopenlots(lt,t)=length(find(leavetype==1));
    end

    for ilv=1:length(ileave)
        iilots=find(Lottype(:,1)==ileave(ilv));
        lotchoice(ileave(ilv),7)=0;
        Lottype(iilots,7)=0;
    end
    iexisthouse=ileave;
    for lt=1:HT
        numlt(lt,t)=length(find(LOTTYPE==lt))/z(lt,1);
    end

    for in=1:length(iNfarmers)
        ifarmland=(AGLAYER == iNfarmers(in));
%         PLAND(ifarmland)=wtaland(iNfarmers(in),t-1);
        PLAND(ifarmland)=LOCWGHT*wtaland(iNfarmers(in),t-1)+REGWGHT*mean(mean(PLAND));
        Dynplandmap(:,:,t)=PLAND;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%    Storm Simulator    %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Possibility of one storm a year
    probsurf=Pflood+siterisk;
    impactprob=rand(1);
    iimpact=(probsurf > impactprob);
    inoimpact=(probsurf <= impactprob);
    stormstats(1,t)=~isempty(find(iimpact,1));
    stormstats(2,t)=length(find(iimpact==1));
    subimpact=zeros(NLENGTH,NWIDTH);
    subimpact(iimpact)=1;
    IMPACT(:,:,t)=subimpact;
    DAMAGE(iimpact)=DAMAGE(iimpact)+Cdam;
    TSI(iimpact)=0;
    TSI(inoimpact)=TSI(inoimpact)+1;
    
    incon=find(Income(:,1)~=0);
    for in=1:length(incon)
        damcoef(:,:,incon(in),t)=(1-CONINFO(incon(in),5))*damcoef(:,:,incon(in),t-1).*...
            (1./(TSI+1))+CONINFO(incon(in),5)*damcoef(:,:,incon(in),t-1);
    end
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
        isamelot=find(Lottype(:,1) == lotchoice(iii,1));
        bidoverx=(Phousebid(:,iii) >= (Paskhouse(iii)));
        extrabid(iii,1)=length(find(bidoverx ==1));
        reldemand(iii,1)=extrabid(iii,1)-newopenlots(lotchoice(iii,5),t);
        XBIDS(Lottype(isamelot,2))=extrabid(iii,1);
    end
    % Rank housing types by their proportion of bids received last period
    trackbids=[extrabid reldemand lotchoice(:,[1:3 5])];
    demandrank=sortrows(trackbids,-1); %[#of_bids(sorted) relative_demand lotid ind lotsize lottype]
    
    bidrank=demandrank(:,[2 6]);
    
    for lt=1:HT
        pctbuildold(lt,t)=sum(demandrank(demandrank(:,6)==lt,1))/...
            sum(demandrank(:,1));
        pctbuildnew(lt,t)=sum(demandrank(demandrank(:,6)==lt,1))/...
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
    nvac=length(iNfarmers);   %number of vacant plots the developer could buy
    iurb=(BASELAYER == 1);
    iurblist=find(iurb==1);
    ivac=(VACLAND ~= 0);
    ivaclist=find(ivac==1);

    for lt=1:HT
        isimlots=(ismember(lotchoice(:,5),simlotrange(lt,1):...
            simlotrange(lt,2))==1 & lotchoice(:,7)~=0);
        isimcells=unique(Lottype(ismember(Lottype(:,1),lotchoice(isimlots,1)),2));
        
        %Regional Stats
%         simlots_income(lt)=mean(Income(con2lot(isimlots,2)));
%         simlots_util(lt)=mean(mean(AVGUTIL(isimcells)));
        simlots_income(lt)=median(Income(con2lot(isimlots,2)));
        simlots_util(lt)=median(median(AVGUTIL(isimcells)));
        simlots_alpha(lt)=mean(mean(avg_alpha(isimcells)));
        simlots_beta(lt)=mean(mean(avg_beta(isimcells)));
        simlots_gamma(lt)=mean(mean(avg_gamma(isimcells)));
        simlots_ampref(lt)=mean(mean(avg_ampref(isimcells)));
        ireglot=(Lottype(:,5)==lt);
        if isempty(find(ireglot,1))==1
            continue
        else
            regionaldist(lt,t)=mean(DISTANCE(Lottype(ireglot,2)));
            regionalrent(lt,t)=mean(EXPTHOUSE(Lottype(ireglot,2)));
        end
    end

%     reg_util=mean(mean(AVGUTIL(iurb)));
%     reg_income=mean(Income(con2lot(ifilled,2)));
    reg_util=median(median(AVGUTIL(iurb)));
    reg_income=median(Income(con2lot(ifilled,2)));
    reg_alpha=mean(mean(avg_alpha(iurb)));
    reg_beta=mean(mean(avg_beta(iurb)));
    reg_gamma=mean(mean(avg_gamma(iurb)));
    reg_ampref=mean(mean(avg_ampref(iurb)));
    
     B=mvregress([ones(length(ifilled),1) lotchoice(ifilled,3) ...
        travelcost(lotchoice(ifilled,2)) coastprox(lotchoice(ifilled,2))],...
        con2lot(ifilled,1));
%     B=regress(con2lot(ifilled,1),[ones(length(ifilled),1) ...
%         log(Income(con2lot(ifilled,2))) lotchoice(ifilled,3) ...
%         travelcost(lotchoice(ifilled,2))]);

    ddist2hznnei=zeros(NLENGTH,NWIDTH);    %distance to horizontal neighbor from icenter
    ddist2vrtnei=zeros(NLENGTH,NWIDTH);
    potentialbuy=find(BASELAYER==0 & SCAPE == 1);
    for tl=1:length(potentialbuy)
        inddist2dev=10000*ones(NLENGTH,NWIDTH);
        [vacrow,vaccol]=ind2sub([NLENGTH NWIDTH],potentialbuy(tl));

        for col=1:NWIDTH
            ddist2hznnei(1:NLENGTH,col)=abs(col-vaccol).*...
                ones(NLENGTH,1);
            for row=1:NWIDTH
                ddist2vrtnei(row,1:NWIDTH)=abs(row-vacrow).*...
                    ones(1,NLENGTH);
                inddist2dev(row,col)=min(sqrt(ddist2hznnei(row,col)^2+...
                    ddist2vrtnei(row,col)^2),inddist2dev(row,col));
            end
        end
        subdist2dev=[inddist2dev(iurb) iurblist];
        sortdevdist=sortrows(subdist2dev,1);
        iclosedev=sortdevdist(:,2);
        
        for lt=1:HT           
            zonecheck=zoning(ZONES(vacrow,vaccol),:);
            if z(lt,1) >= zonecheck(1) && z(lt,1) <= zonecheck(2)
                warning('off','all');
                likelotcount=round(numlt(lt,t));
                icountcells=iclosedev(1:round(length(iurblist)*PCTSEARCH));
                iutilcells=icountcells(AVGUTIL(icountcells)~=0 & ...
                    (ismember(LOTTYPE(icountcells),simlotrange(lt,1):...
                    simlotrange(lt,2))==1));
                distutils=inddist2dev(iutilcells);
                if isempty(distutils)==1
                    termrunflag=1;
                end
%                 avg_utils=sum(-distutils.*AVGUTIL(iutilcells))/sum(-distutils);
                avg_utils=median(AVGUTIL(iutilcells));
                avgalpha=sum(-distutils.*avg_alpha(iutilcells))/sum(-distutils);
                avgbeta=sum(-distutils.*avg_beta(iutilcells))/sum(-distutils);
                avggamma=sum(-distutils.*avg_gamma(iutilcells))/sum(-distutils);
                avgampref=sum(-distutils.*avg_ampref(iutilcells))/sum(-distutils);
                avgmoney=sum(-distutils.*avg_income(iutilcells))/sum(-distutils);
                warning('on','all');
                if likelotcount == 0
                    
%                     regrentproj=simlots_income(lt)-travelcost(vacrow,vaccol)-...
%                         (simlots_util(lt)/((z(lt,2)^simlots_beta(lt))*...
%                         (z(lt,1)^simlots_gamma(lt))*(z(lt,3)^simlots_ampref(lt))))^...
%                         (1/simlots_alpha(lt));
% %                     locrentproj=avgmoney-travelcost(vacrow,vaccol)-...
% %                         (avg_utils/((z(lt,2)^avgbeta)*(z(lt,1)^avggamma)*...
% %                         (z(lt,3)^avgampref)))^(1/avgalpha);
%                     locrentproj=simlots_income(lt)-travelcost(vacrow,vaccol)-...
%                         (avg_utils/((z(lt,2)^avgbeta)*(z(lt,1)^avggamma)*...
%                         (z(lt,3)^avgampref)))^(1/avgalpha);
% 
%                     subRENTPROJ(vacrow,vaccol,lt)=LOCWGHT*locrentproj+...
%                         REGWGHT*regrentproj;
%                     subRENTPROJ(vacrow,vaccol,lt)=B(1)+...
%                         B(2)*log(simlots_income(lt))+B(3)*z(lt,1)+...
%                         B(4)*travelcost(vacrow,vaccol);
                    subRENTPROJ(vacrow,vaccol,lt)=B(1)+...
                        B(2)*z(lt,1)+B(3)*travelcost(vacrow,vaccol)+...
                        B(4)*coastprox(vacrow,vaccol);
                else
                    nclosecells=max(likelotcount*min(z(lt,1),1),round(length(iurblist)*PCTSEARCH));
                    isamecell=ismember(lotchoice(:,2),iclosedev(1:nclosecells));
                    iclosecells=find(isamecell==1);
                    icloselot=find(lotchoice(iclosecells,5) == lt,nclosecells,'first');
                    devrents=EXPTHOUSE(lotchoice(iclosecells(icloselot),2));
                    distrents=inddist2dev(lotchoice(iclosecells(icloselot),2));
                    [closecelldist,iclstcell]=min(distrents,[],1);
                    iclosestcell=lotchoice(iclosecells(icloselot(iclstcell)),2);
                    
                    if isempty(icloselot)==1
                        locrentproj=B(1)+B(2)*z(lt,1)+B(3)*...
                            travelcost(vacrow,vaccol)+B(4)*coastprox(vacrow,vaccol);
%                         locrentproj=B(1)+...
%                             B(2)*log(simlots_income(lt))+B(3)*z(lt,1)+...
%                             B(4)*travelcost(vacrow,vaccol);
                        regrentproj=regionalrent(lt,t)-margtc*...
                            (DISTANCE(vacrow,vaccol)-regionaldist(lt,t));
                        subRENTPROJ(vacrow,vaccol,lt)=LOCWGHT*locrentproj+...
                            REGWGHT*regrentproj;
%                         
                    else
                        locrentproj=(sum(-distrents.*devrents)/sum(-distrents))-...
                            margtc*distrents(iclstcell);
                        regrentproj=regionalrent(lt,t)-margtc*...
                            (DISTANCE(vacrow,vaccol)-regionaldist(lt,t));
                        subRENTPROJ(vacrow,vaccol,lt)=LOCWGHT*locrentproj+...
                            REGWGHT*regrentproj;
                    end
                end

                RETURN(vacrow,vaccol,lt)=(1-discount)*(subRENTPROJ(vacrow,...
                    vaccol,lt)-ccost(lt))/z(lt,1)-carrycost(t-1)/...
                    ceil(sum(newhouseset(:,t).*z(:,1)));

                highRETURN(vacrow,vaccol,lt)=(1-discount)*((overvalue(lt,t-1)+...
                    subRENTPROJ(vacrow,vaccol,lt))-ccost(lt))/z(lt,1)-carrycost(t-1)/...
                    ceil(sum(newhouseset(:,t).*z(:,1)));
                lowRETURN(vacrow,vaccol,lt)=(1-discount)*((subRENTPROJ(vacrow,vaccol,lt)+...
                    undervalue(lt,t-1))-ccost(lt))/z(lt,1)-carrycost(t-1)/...
                    ceil(sum(newhouseset(:,t).*z(:,1)));
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 % Gain relative to return
%                 if highRETURN(vacrow,vaccol,lt) > RETURN(vacrow,vaccol,lt)
%                     potgain(vacrow,vaccol,lt)=(abs(highRETURN(vacrow,vaccol,lt)...
%                         -RETURN(vacrow,vaccol,lt))/alpha_gain)^(1/alpha_gain);
%                 elseif highRETURN(vacrow,vaccol,lt) <= RETURN(vacrow,vaccol,lt)
%                     potgain(vacrow,vaccol,lt)=0;
% %                     potloss(vacrow,vaccol,lt)=(abs(RETURN(vacrow,vaccol,lt)-...
% %                         highRETURN(vacrow,vaccol,lt))/alpha_loss)^(1/alpha_loss);
%                 end
%                 if lowRETURN(vacrow,vaccol,lt) < RETURN(vacrow,vaccol,lt)
%                     potloss(vacrow,vaccol,lt)=potloss(vacrow,vaccol,lt)+...
%                         (abs(RETURN(vacrow,vaccol,lt)-abs(lowRETURN(vacrow,vaccol,lt)))/...
%                         alpha_loss)^(1/alpha_loss);
%                 elseif lowRETURN(vacrow,vaccol,lt) >= RETURN(vacrow,vaccol,lt)
%                     potloss(vacrow,vaccol,lt)=0;
% %                     potgain(vacrow,vaccol,lt)=potgain(vacrow,vaccol,lt)+...
% %                         (abs(lowRETURN(vacrow,vaccol,lt)-RETURN(vacrow,vaccol,lt))/...
% %                         alpha_gain)^(1/alpha_gain);
%                 end     
                
                % Gain relative to zero
                if highRETURN(vacrow,vaccol,lt) > 0
                    potgain(vacrow,vaccol,lt)=(highRETURN(vacrow,vaccol,lt)/...
                        alpha_gain)^(1/alpha_gain);
                elseif highRETURN(vacrow,vaccol,lt) <= 0
                    potgain(vacrow,vaccol,lt)=0;
                    potloss(vacrow,vaccol,lt)=(abs(highRETURN(vacrow,vaccol,lt))/...
                        alpha_loss)^(1/alpha_loss);
                end
                if lowRETURN(vacrow,vaccol,lt) < 0
                    potloss(vacrow,vaccol,lt)=potloss(vacrow,vaccol,lt)+...
                        (abs(lowRETURN(vacrow,vaccol,lt))/alpha_loss)^(1/alpha_loss);
                elseif lowRETURN(vacrow,vaccol,lt) >= 0
                    potloss(vacrow,vaccol,lt)=0;
                    potgain(vacrow,vaccol,lt)=potgain(vacrow,vaccol,lt)+...
                        (abs(lowRETURN(vacrow,vaccol,lt))/alpha_gain)^(1/alpha_gain);
                end   
                
                EU_dev(vacrow,vaccol,lt)=(potgain(vacrow,vaccol,lt))/...
                    (potgain(vacrow,vaccol,lt)+potloss(vacrow,vaccol,lt));
                
                %                 PCHANGE(vacrow,vaccol,lt)=RETURN(vacrow,vaccol,lt)/...
%                     (risk*avgbrokervar(lt,t-1));
            else
                RETURN(vacrow,vaccol,lt)=0;
            end
            subret=RETURN(:,:,lt);
            Exptret(lt,t)=mean(subret(potentialbuy));
        end
%         [maxeuval,maxretind]=max(EU_dev(vacrow,vaccol,:),[],3);
%         MAXRET(vacrow,vaccol,t)=RETURN(vacrow,vaccol,maxretind);
        maxcount(vacrow,vaccol)=length(find(EU_dev(vacrow,vaccol,:) == ...
            max(EU_dev(vacrow,vaccol,:))));
    end
    if isempty(find(newhouseset(:,t),1))==1
        LANDBUDGET(t)=0;
    else
        LANDBUDGET(t)=BUDGET(t-1)/round(sum(newhouseset(:,t).*z(:,1)));
    end
%     houseset=max(sum(numnewhouses(t)+numoldhouses(t)),0);
    
%     [maxpchg,pchgind]=max(PCHANGE(:,:,1:HT),[],3);
    
    indfarms=find(iagr==1);
    farmprojcells=zeros(length(iagr),3);     %[profitability ltype ind]
    [maxret,retind]=max(RETURN(:,:,1:HT),[],3);
    [EU_value,EU_ind]=sort(EU_dev,3,'descend');
    [maxeuval,maxeuind]=max(EU_dev,[],3);
    maxeuind(maxcount > 1)=retind(maxcount > 1);
    for tl=1:length(potentialbuy)
        [vacrow,vaccol]=ind2sub([NLENGTH NWIDTH],potentialbuy(tl));
%         WTPMAP(vacrow,vaccol,t)=min(RETURN(vacrow,vaccol,retind(vacrow,vaccol))/discount,...
%             LANDBUDGET(vacrow,vaccol,retind(vacrow,vaccol)));
        EUrankret(vacrow,vaccol,:)=RETURN(vacrow,vaccol,EU_ind(vacrow,vaccol,:));
%         if maxcount(vacrow,vaccol) > 1
%             subEUinfo=zeros([],1);
%             subRETinfo=zeros([],1);
%             subEUinfo(1:maxcount(vacrow,vaccol),1)=EU_ind(vacrow,vaccol,1:maxcount(vacrow,vaccol));
%             subRETinfo(1:maxcount(vacrow,vaccol),1)=...
%                 RETURN(vacrow,vaccol,EU_ind(vacrow,vaccol,1:maxcount(vacrow,vaccol)));
%             maxEUinfo=[subEUinfo subRETinfo];
%             sortEUinfo=sortrows(maxEUinfo,-2);
%             EUrankret(vacrow,vaccol,1:maxcount(vacrow,vaccol))=sortEUinfo(:,2);
%             EU_ind(vacrow,vaccol,1:maxcount(vacrow,vaccol))=sortEUinfo(:,1);
%         end
%         iposeu=find(RETURN(vacrow,vaccol,EU_ind(vacrow,vaccol,:)) >= 0);
        iposeu=find(EUrankret(vacrow,vaccol,:) > 0);
        poseuset=reshape(EU_ind(vacrow,vaccol,iposeu),length(iposeu),1);
        eubuildset=max(newhouseset(poseuset,t),1);
%         MAXRET(vacrow,vaccol,t)=sum(z(poseuset,1).*eubuildset.*reshape(RETURN(vacrow,vaccol,...
%             EU_ind(vacrow,vaccol,iposeu)),length(iposeu),1))/sum(z(poseuset,1).*eubuildset);
        MAXRET(vacrow,vaccol,t)=sum(reshape(EU_dev(vacrow,vaccol,poseuset),...
            length(iposeu),1).*reshape(RETURN(vacrow,vaccol,...
            EU_ind(vacrow,vaccol,iposeu)),length(iposeu),1))/...
            sum(reshape(EU_dev(vacrow,vaccol,poseuset),length(iposeu),1));
        %         MAXRET(vacrow,vaccol,t)=RETURN(vacrow,vaccol,maxeuind(vacrow,vaccol));
%         WTPMAP(vacrow,vaccol,t)=min(RETURN(vacrow,vaccol,retind(vacrow,vaccol))/discount,...
%             LANDBUDGET(t));

        WTPMAP(vacrow,vaccol,t)=MAXRET(vacrow,vaccol,t)/discount;
%         WTPMAP(vacrow,vaccol,t)=MAXRET(vacrow,vaccol,t)/...
%             (discount*(1+cvht(maxeuind(vacrow,vaccol),t-1)));
        
    end
    farmprojcells(1:length(indfarms),:)=[maxret(iagr) retind(iagr) indfarms];
    RETURNPROJ=maxret;

    [maxEU_dev,eudevind]=max(EU_dev(:,:,1:HT),[],3);
    
    Dynmaxretmap(:,:,t)=maxret;
    Dynretltmap(:,:,t)=retind;
    Dyneultmap(:,:,t)=maxeuind;
    
    maxretcells=zeros([],3);
    iundev=(BASELAYER == 0);
    indundev=find(iundev==1);
    maxretcells(length(maxretcells(:,1))+1:length(maxretcells(:,1))+length(indundev),:)=...
        [RETURNPROJ(iundev) ones(length(indundev),1).*retind(iundev) ...
        indundev];
    rankmpcells=sortrows(maxretcells,-1);
%     numnewacres(t)=sum(z(rankmpcells(1:min(houseset,length(rankmpcells(:,1))),2),1));
    numnewacres(t)=ceil(sum(newhouseset(:,t).*z(:,1)));
    
    for ww=1:NWIDTH
        for ll=1:NLENGTH
            maxRENTPROJ(ll,ww)=subRENTPROJ(ll,ww,retind(ll,ww));
            CCOST(ll,ww)=ccost(retind(ll,ww));
            ZZ(ll,ww)=z(retind(ll,ww),1);
        end
    end
    
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%<><><><><><><><><><> Go to Land Market <><><><><><><><><><><><><><><><><><>
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
   
%     landdemand(1,t)=max(sum(numnewacres(t))-length(find(ivac==1)),0);
    landdemand(1,t)=max(sum(numnewacres(t))-max(length(find(ivac==1)),...
        min(Farmstats(iNfarmers,2,t))),0);
%     subMAXRET=MAXRET(:,:,t);
    for rf=1:length(iNfarmers)       
        ifarmcalc=find(AGLAYER == iNfarmers(rf));
        subWTPMAP=WTPMAP(:,:,t);
        MAXEUMAP=MAXRET(:,:,t);
        
        wtpland(iNfarmers(rf),t)=sum(subWTPMAP(ifarmcalc))/length(ifarmcalc);
    

    end
    
    LandMarket_6v0
    
    for ft=1:length(ifarmtrans)
        VACLAND(AGLAYER == ifarmtrans(ft))=1;
    end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ivac=(VACLAND ~= 0);
    indvac=find(ivac==1);
    vcells=zeros([],4);
        
    for tl=1:length(indvac)
        subrentEU=zeros([],1);
        [vacrow,vaccol]=ind2sub([NLENGTH NWIDTH],indvac(tl));
        subrentEU(1:HT,1)=subRENTPROJ(vacrow,vaccol,:);
        EUlandret(vacrow,vaccol,:)=(subrentEU-((PLAND(vacrow,vaccol).*z(:,1).*...
            discount)+ccost))./z(:,1);
%         EUlandret(vacrow,vaccol,:)=(subrent-((PLAND(vacrow,vaccol).*z(:,1).*...
%             (discount*(1+cvht(:,t-1))))+ccost))./z(:,1);
        EUrankret(vacrow,vaccol,:)=EUlandret(vacrow,vaccol,EU_ind(vacrow,vaccol,:));
    end
    
    for ilayer=1:HT
        subEUprofit=EUrankret(:,:,ilayer);
        subeuind=EU_ind(:,:,ilayer);
        posret=(subEUprofit(indvac)>=0);
%         vcells(length(vcells(:,1))+1:length(vcells(:,1))+length(indvac(posret)),:)=...
%             [subEUprofit(indvac(posret)) subeuind(indvac(posret)) indvac(posret)];
        subprofit=RETURN(:,:,ilayer)-PLAND.*discount;
%         subprofit=RETURN(:,:,ilayer)-PLAND.*(discount*(1+cvht(ilayer,t-1)));
        iltypes=unique(subeuind);
        for ip=1:length(iltypes)
            Exptprofit(iltypes(ip),t)=mean(subprofit(subeuind==iltypes(ip)));
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
    PROFIT(:,:,t)=MAXRET(:,:,t)-(PLAND.*discount);
%     PROFIT(:,:,t)=MAXRET(:,:,t)-(PLAND.*...
%             (discount*(1+reshape(cvht(maxeuind,t-1),NLENGTH,NWIDTH))));
    subprofit=PROFIT(:,:,t);
    subprofit(iurb)=-1;
    
%     vcells=vcells(vcells(:,1) > 0,:);
%     rankvaccells=sortrows(vcells,-1);
    rankvaccells=vcells(:,1:3);
  
    Dynmaxprofmap(:,:,t)=PROFIT(:,:,t);
    Dynprofltmap(:,:,t)=retind;
    

%     newhouseset(:,t)=round(houseset.*(bidtot(:,t)./sum(bidtot(:,t))));

    %%%%%%%%%%%%%%% Profit-Driven Construction decisions %%%%%%%%%%%%%%%%%%
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
                    novac=(ismember(isnei(sc,:),Lottype(:,2))==0);
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
                    BASELAYER(totnewbuildcells)=1;
                    VACLAND(totnewbuildcells)=0;
%                     Lottype(length(Lottype(:,1))+1:length(Lottype(:,1))+...
%                         length(totnewbuildcells),:)=[ones(length(iadded),1).*(lotidmark+1) ...
%                         totnewbuildcells ones(length(iadded),1)*[z(newlottype,1) ...
%                         z(newlottype,2) newlottype ccost(newlottype)] ...
%                         zeros(length(iadded),1) ones(length(iadded),1)*z(newlottype,3)];                     
                    Lottype(length(Lottype(:,1))+1:length(Lottype(:,1))+...
                        length(totnewbuildcells),:)=[ones(length(iadded),1).*(lotidmark+1) ...
                        totnewbuildcells ones(length(iadded),1)*[z(newlottype,1) ...
                        z(newlottype,2) newlottype ccost(newlottype)] ...
                        zeros(length(iadded),1) coastprox(totnewbuildcells)];
                    newquota=newquota-1;
                    lotidmark=lotidmark+1;
                end
                newlottot=newlottot+length(iadd);
            end
        elseif newlotsize < 1
            totnewbuildcells=rankvaccells(istartcell,3);
            BUILDTIME(totnewbuildcells)=t;
            BASELAYER(totnewbuildcells)=1;
            VACLAND(totnewbuildcells)=0;
            namelength=1/newlotsize;
%             Lottype(length(Lottype(:,1))+1:length(Lottype(:,1))+namelength,:)=...
%                 [(lotidmark+(1:namelength))' (([rankvaccells(istartcell,3) ...
%                 z(rankvaccells(istartcell,2),1) z(rankvaccells(istartcell,2),2) ...
%                 newlottype ccost(rankvaccells(istartcell,2)) 0 ...
%                 z(rankvaccells(istartcell,2),3)])'*ones(1,namelength))'];
            Lottype(length(Lottype(:,1))+1:length(Lottype(:,1))+namelength,:)=...
                [(lotidmark+(1:namelength))' (([rankvaccells(istartcell,3) ...
                z(rankvaccells(istartcell,2),1) z(rankvaccells(istartcell,2),2) ...
                newlottype ccost(rankvaccells(istartcell,2)) 0 ...
                coastprox(totnewbuildcells)])'*ones(1,namelength))'];
            newquota=newquota-namelength;
            lotidmark=lotidmark+namelength;
        else
            totnewbuildcells=rankvaccells(istartcell,3);
            BUILDTIME(totnewbuildcells)=t;
            BASELAYER(totnewbuildcells)=1;
            VACLAND(totnewbuildcells)=0;
%             Lottype(length(Lottype(:,1))+1,:)=[lotidmark+1 rankvaccells(istartcell,3)...
%                 z(rankvaccells(istartcell,2),1) z(rankvaccells(istartcell,2),2)...
%                 newlottype ccost(rankvaccells(istartcell,2)) 0 z(rankvaccells(istartcell,2),3)];
            Lottype(length(Lottype(:,1))+1,:)=[lotidmark+1 rankvaccells(istartcell,3)...
                z(rankvaccells(istartcell,2),1) z(rankvaccells(istartcell,2),2)...
                newlottype ccost(rankvaccells(istartcell,2)) 0 coastprox(totnewbuildcells)];
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
        
        iurb=(BASELAYER == 1);
        iagr=(BASELAYER == 0);
        iurblist=find(iurb==1);
        ivac=(VACLAND ~= 0);
        ivaclist=find(ivac==1);
        iscape=(SCAPE == 1);

        if isempty(find(ivac,1))
            break
        end
        if isempty(find(newhousequota,1))==1
            break
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VACLAND(iurb)=0;
    ivac=(VACLAND ~= 0);
    ivaclist=find(ivac==1);

    lotchoice=zeros(Nlots(1),8);
    cellinfo=zeros(length(iurblist),8);
    [lotids,ilots,jnum]=unique(Lottype(:,1));
    lotchoice(lotids,:)=Lottype(ilots,:);
    Nlots(t)=length(lotchoice(:,1));

    [cellids,icells,jnumcell]=unique(Lottype(:,2));
    cellinfo(1:length(iurblist),:)=Lottype(icells,:);
    LOTTYPE(cellinfo(:,2))=cellinfo(:,5);



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
    if ntrans > 0
        for nt=1:ntrans
            tdist2hznnei=10000*ones(NLENGTH,NWIDTH);    %distance to horizontal neighbor from icenter
            tdist2vrtnei=10000*ones(NLENGTH,NWIDTH);
            transland=find(AGLAYER==ifarmtrans(nt));
            for tl=1:length(transland)
                [itlrow,itlcol]=ind2sub([NLENGTH,NWIDTH],transland(tl));
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
            AGLAYER(transland)=0;
            iremain=(AGLAYER ~= 0);
            iNfarmers=unique(AGLAYER(iremain));
            %Calculate gradient coefficients using genetic algorithm
            subtransdist=indtransdist(:,:,nt);
            distcoeff(ifarmtrans,:)=0;
            for nf=1:length(iNfarmers)
                    rc=zeros([],1);
                    avgtransdist=mean(subtransdist(AGLAYER==iNfarmers(nf)));
                    coeffmark=(Paskland(iNfarmers(nf),t)-mean(PLAND(transland)))/...
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
                mean(PLAND(transland));
        end
        for nnw=1:NWIDTH
            for nnl=1:NLENGTH
                transdist(nnl,nnw)=min(indtransdist(nnl,nnw,1:ntrans));
            end
        end
        Plandproj(iNfarmers,t)=mean(Planddistproj(iNfarmers,:),2);
    else
%         Plandproj(iNfarmers,t)=mean([wtaland(iNfarmers,t) zeta*wtpland(iNfarmers,t)],2);
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
%                     pcoef=polyfit(1:length(indata),indata,2);
%                     pline=pcoef(1).*(1:length(indata)+1).^2+pcoef(2).*(1:...
%                         length(indata)+1)+pcoef(3);
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
            Farmstats(iNfarmers(nf),5,1));
    end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    Consumers' Choices    %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Add new consumers %%%%%%%%%%
    %POP(t)=ceil(Nconsumers*(1+POPGROW));
    newpop=ceil(Nconsumers*POPGROW);
    inewpop=(length(Income)+1:length(Income)+newpop);
    Income(inewpop,1)=min(max(normrnd(avgincome,stdincome,newpop,1),minwage),maxwage);
    Income(inewpop,2)=t+ceil(searchtimemin+(searchtimemax-searchtimemin)*rand(length(inewpop),1));
    
    incon=find(Income(:,1)~=0);

    newincome3=(Income(inewpop,1) >= lowminwage & Income(inewpop,1) <= lowmaxwage);     %income 1 2 3 = hi mid low wages
    newincome2=(Income(inewpop,1) >= midminwage & Income(inewpop,1) <= midmaxwage);
    newincome1=(Income(inewpop,1) >= himinwage & Income(inewpop,1) <= himaxwage);
    housepref(inewpop,1)=0;
    housepref(inewpop(newincome1),1)=0.20+(0.24-0.20)*rand(length(find(newincome1==1)),1);
    housepref(inewpop(newincome2),1)=0.25+(0.29-0.25)*rand(length(find(newincome2==1)),1);
    housepref(inewpop(newincome3),1)=0.30+(0.50-0.30)*rand(length(find(newincome3==1)),1);
    
    Hbeta=housepref(inewpop);
    beta=(0.1+(0.9-0.1)*rand(length(Hbeta),1)).*Hbeta;     %cite Carliner
    gamma=rand(length(Hbeta),1).*(Hbeta-beta);   %cite Carliner
    ampref=Hbeta-(beta+gamma);
    alpha=(1-Hbeta);                    %" " " " " consumer 
%     subrisk=ones(length(Hbeta),1);
    subrisk=rand(length(Hbeta),1);
    CONINFO(inewpop,:)=0;
    CONINFO(inewpop,1)=alpha;
    CONINFO(inewpop,2)=beta;
    CONINFO(inewpop,3)=gamma;
    CONINFO(inewpop,4)=ampref;
    CONINFO(inewpop,5)=subrisk;
    damcoef(:,:,inewpop,t)=1;
%     for in=1:length(incon)
%         damcoef(:,:,incon(in),t)=(1-CONINFO(incon(in),5))*damcoef(:,:,incon(in),t-1).*...
%             (1./(TSI+1))+CONINFO(incon(in),5)*damcoef(:,:,incon(in),t-1);
%     end
    
    income3=(Income(:,1) >= lowminwage & Income(:,1) <= lowmaxwage);     %income 1 2 3 = hi mid low wages
    income2=(Income(:,1) >= midminwage & Income(:,1) <= midmaxwage);
    income1=(Income(:,1) >= himinwage & Income(:,1) <= himaxwage);

%<><><><><><><><>  Matching consumers with houses  <><><><><><><><><><><><>
    minprof=zeros(Nlots(t),1);
    inewlots=find(lotchoice(:,7)==0);
    newcons=find(con2lot(:,2)~=0 & con2lot(:,4)==t);
    oldreturn=(Income(ioldpop,1)~=0);
    inewcon=sort([con2lot(newcons,2)' ioldpop(oldreturn)' inewpop]','ascend');

    for nl=1:length(inewlots)
        isamelot=find(Lottype(:,1)==lotchoice(inewlots(nl),1));
        subprof=PROFIT(:,:,lotchoice(inewlots(nl),5));
        RENTPROJ=subRENTPROJ(:,:,lotchoice(inewlots(nl),5));
        if BUILDTIME(lotchoice(inewlots(nl),2)) < t
            Paskhouse(inewlots(nl))=con2lot(inewlots(nl),1);
        else
            Paskhouse(inewlots(nl))=max(mean(RENT(Lottype(isamelot,2))),...
                mean(RENTPROJ(Lottype(isamelot,2))));
        end
    end

    con2lot(ileave,2)=0;
    con2lot(Nlots(t-1)+1:Nlots(t),:)=0;     %[WinBid conid lotid restime]
    con2lot(Nlots(t-1)+1:Nlots(t),3)=Nlots(t-1)+1:Nlots(t);
    
    %<><><><><><><><><><> Go to housing market <><><><><><><><><><><><><><><><>
    
    HouseMarketDynamic_coast_v1_subrisk

    house2cell(iurb)=1./z(LOTTYPE(iurb),1);
    RENTGRAD(iurb)=RENT(iurb).*house2cell(iurb);
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    popinhouse=ismember(1:length(Income),con2lot(:,2))';
    ioldpop=find(popinhouse==0 & Income(:,1)~=0);
    oldincome(t)=mean(Income(ioldpop));
    RENTPROJ(iurb)=0;
    Nlots(t+1)=Nlots(t);
    for lt=1:HT
        isamelot=(lotchoice(istillvac,5)==lt);
        vacstats(lt,t)=length(find(isamelot==1));
    end
    devcells(1,t)=length(iurblist);
    devcells(2,t)=devcells(1,t)/(NLENGTH*NWIDTH);
    
    vacantlots(t)=length(istillvac);
    leftoverpop(t)=length(ioldpop);
    
    vacrate(t)=vacantlots(t)/Nlots(t);
    nohouserate(t)=leftoverpop(t)/POP(t);
    
    %Consumers leave area
    iconleave=find(Income(:,2) == t);
    OutIncome(length(OutIncome)+1:length(OutIncome)+length(iconleave),1:2)=...
        Income(iconleave,:);
    Income(iconleave,:)=0;
    Nconsumers=length(find(Income(:,1)~=0));
    incon=find(Income(:,1)~=0);
    
    subrealret=zeros(NLENGTH,NWIDTH);
    subexptrent=zeros(NLENGTH,NWIDTH);
    subnewlt=zeros(NLENGTH,NWIDTH);
    subnewbid=zeros(NLENGTH,NWIDTH);
    subrealexptret=zeros(NLENGTH,NWIDTH);

    for nl=1:length(inewlots)
        isamelot=find(Lottype(:,1)==lotchoice(inewlots(nl),1));
        isamecell=find(Lottype(:,2)==Lottype(isamelot(1),2));
        if length(isamecell) > 1
            subexptrent(Lottype(isamelot,2))=mean(con2lot(Lottype(isamecell,1),1)-...
                Paskhouse(Lottype(isamecell,1)));
            subrealret(Lottype(isamelot,2))=mean((con2lot(Lottype(isamecell,1),1)...
                -((PLAND(Lottype(isamecell,2)).*z(Lottype(isamecell,5),1)*discount)+...
                ccost(lotchoice(inewlots(nl),5)))));
            subrealexptret(Lottype(isamelot,2))=mean((Paskhouse(Lottype(isamecell,1))...
                -((PLAND(Lottype(isamecell,2)).*z(Lottype(isamecell,5),1)*discount)+...
                ccost(lotchoice(inewlots(nl),5)))));
            subnewbid(Lottype(isamelot,2))=mean(con2lot(Lottype(isamecell,1),1)./...
                Paskhouse(Lottype(isamecell,1)));
        else
        subexptrent(Lottype(isamelot,2))=con2lot(inewlots(nl),1)-Paskhouse(inewlots(nl));
        subrealret(Lottype(isamelot,2))=(con2lot(inewlots(nl),1)-...
            ((PLAND(lotchoice(inewlots(nl),2)).*z(lotchoice(inewlots(nl),5),1)*discount)+...
            ccost(lotchoice(inewlots(nl),5))));   
        subrealexptret(Lottype(isamelot,2))=(Paskhouse(inewlots(nl))-...
            ((PLAND(lotchoice(inewlots(nl),2)).*z(lotchoice(inewlots(nl),5),1)*discount)+...
            ccost(lotchoice(inewlots(nl),5))));
        subnewlt(Lottype(isamelot,2))=lotchoice(inewlots(nl),5);
        subnewbid(Lottype(isamelot,2))=con2lot(inewlots(nl),1)/Paskhouse(inewlots(nl));
        end
    end
    Exptrentdiff(:,:,t)=subexptrent;
    Realreturn(:,:,t)=subrealret;
    Realexptret(:,:,t)=subrealexptret;
    vac_land(t)=sum(PLAND(ivac)*discount);
    for lt=1:HT
        avgrent(lt,t)=mean(RENT(LOTTYPE==lt));
        ifindlt=(lotchoice(inewlots,5)==lt);
        Avgexptdiff(lt,t)=mean(subexptrent(lotchoice(inewlots(ifindlt),2)));
        Realavgret(lt,t)=mean(subrealret(lotchoice(inewlots(ifindlt),2)));
        Realavgexptret(lt,t)=mean(subrealexptret(lotchoice(inewlots(ifindlt),2)));
        Avgnewbid(lt,t)=mean(subnewbid(lotchoice(inewlots(ifindlt),2)));
        ilt=find(lotchoice(istillvac,5)==lt);

%         vac_ccost(lt,t)=sum(discount*lotchoice(istillvac(ilt),6));
%         vac_ccost(lt,t)=sum(discount*lotchoice(istillvac(ilt),6).*(t-...
%             BUILDTIME(lotchoice(istillvac(ilt),2))));

%         vac_ccost(lt,t)=sum(discount*(lotchoice(istillvac(ilt),6)+...
%             PLAND(lotchoice(istillvac(ilt),2)).*z(lotchoice(istillvac(ilt),5),1)));
        vac_ccost(lt,t)=sum(discount*(lotchoice(istillvac(ilt),6)+...
            PLAND(lotchoice(istillvac(ilt),2)).*z(lotchoice(istillvac(ilt),5),1)));
%         vac_ccost(lt,t)=sum((discount*(1+cvht(lt,t-1)))*(lotchoice(istillvac(ilt),6)+...
%             PLAND(lotchoice(istillvac(ilt),2)).*z(lotchoice(istillvac(ilt),5),1)));

        vac_rent(lt,t)=sum(Paskhouse(istillvac(ilt)));
        
        ilt_t=find(lotchoice(:,5) == lt & BUILDTIME(lotchoice(:,2)) == t);
        profits(lt,t)=sum(con2lot(ilt_t,1)-(PLAND(lotchoice(ilt_t,2)).*...
            z(lotchoice(ilt_t,5),1)*discount+lotchoice(ilt_t,6)));
%         profits(lt,t)=sum(con2lot(ilt_t,1)-(PLAND(lotchoice(ilt_t,2)).*...
%             z(lotchoice(ilt_t,5),1)*(discount*(1+cvht(lt,t-1)))+lotchoice(ilt_t,6)));
%         budget_lt(lt,t)=budget_lt(lt,t-1)+profits(lt,t)-vac_ccost(lt,t)-vac_land(t);
%         budget_lt(lt,t)=budget_lt(lt,t-1)+profits(lt,t)-vac_ccost(lt,t);
        budget_lt(lt,t)=profits(lt,t)-vac_ccost(lt,t);
    end
    Newlottype(:,:,t)=subnewlt;
    Newbidlevel(:,:,t)=subnewbid;
    carrycost(t)=sum(vac_ccost(:,t))+vac_land(t);
    
    BUDGET(t)=BUDGET(t-1)+sum(budget_lt(:,t))-vac_land(t);          
    
%     figure(4)
%     surf(RENT);
%     axis ij;
%     view(0,90);
%     title(sprintf('RENT, t = %d',t))
%     colorbar
%     M(t)=getframe(gcf);
%      
%     figure(5)
%     surf(RENTGRAD);
%     axis ij;
%     view(0,90);
%     title(sprintf('RENT/acre, t = %d',t))
%     colorbar
%     MG(t)=getframe(gcf);
    
%     figure(6)
%     surf(PLAND);
%     axis ij;
%     view(0,90);
%     title(sprintf('PLAND, t = %d',t))
%     colorbar
%     MPL(t)=getframe(gcf);
%     
%     figure(7)
%     surf(LOTTYPE);
%     axis ij;
%     view(0,90);
%     title(sprintf('Lot Types, t=%d',t))
%     set(gca,'clim',[1 18])
%     colorbar
%     MLT(t)=getframe(gcf);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%    Brokers' Price Projections    %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    warning('off','all');
    for nb=1:Nbrokers
        str= sprintf('find(HBROKER == %d);',nb);
        ibroker=eval(str);
        ibcellind=ismember(Lottype(:,2),ibroker);
        ibsample=Lottype(ibcellind,[1:2 5 7]);
        
        if isempty(ibsample) == 1
            continue
        else
            [iblots,ib,jb]=unique(ibsample(:,1));
            nbids=zeros([],1);
            inewbids=find(ismember(inewlots,iblots)==1);
            for b=1:length(iblots)
                bidover=(Phousebid(:,iblots(b)) > 0);
                nbids(b,1)=length(find(bidover==1));
            end

            for sc=1:length(ibroker)
                isamecell=find(ibsample(:,2)==ibroker(sc));
                isnotvac=(ibsample(isamecell,4) ~= 0);
                if isempty(find(ibsample(isamecell(isnotvac),1),1)) == 1
                    continue
                else
%                     avg_income(ibroker(sc))=mean(Income(con2lot(ibsample(isamecell(isnotvac),1),2)));
%                     avg_alpha(ibroker(sc))=mean(CONINFO(con2lot(ibsample(isamecell(isnotvac),1),2),1));
                    avg_income(ibroker(sc))=median(Income(con2lot(ibsample(isamecell(isnotvac),1),2)));
                    avg_alpha(ibroker(sc))=median(CONINFO(con2lot(ibsample(isamecell(isnotvac),1),2),1));                    
                    avg_beta(ibroker(sc))=mean(CONINFO(con2lot(ibsample(isamecell(isnotvac),1),2),2));
                    avg_gamma(ibroker(sc))=mean(CONINFO(con2lot(ibsample(isamecell(isnotvac),1),2),3));
                    avg_ampref(ibroker(sc))=mean(CONINFO(con2lot(ibsample(isamecell(isnotvac),1),2),4));
                    AVGUTIL(ibroker(sc))=((avg_income(ibroker(sc))-RENT(ibroker(sc))-...
                        travelcost(ibroker(sc)))^avg_alpha(ibroker(sc)))*...
                        (mean(z(ibsample(isamecell(isnotvac),3),2))^avg_beta(ibroker(sc)))*...
                        (mean(z(ibsample(isamecell(isnotvac),3),1))^avg_gamma(ibroker(sc)))*...
                        (mean(coastprox(ibroker(sc)))^avg_ampref(ibroker(sc)));
                end
            end
            
            subrent=con2lot(iblots,1);
            subrent(con2lot(iblots,2)==0)=0;
            sampleinfo=zeros(length(iblots),6);
            sampleinfo(:,1)=subrent;
            sampleinfo(:,2)=lotchoice(iblots,3);
            sampleinfo(:,3)=lotchoice(iblots,4);
            sampleinfo(:,4)=nbids;
            sampleinfo(:,5)=con2lot(iblots,1)./(ccost(lotchoice(iblots,5))+...
                discount.*Farmstats(FARMERPOS(lotchoice(iblots,2)),5,1).*...
                z(lotchoice(iblots,5),1));
%             sampleinfo(:,5)=con2lot(iblots,1)./(ccost(lotchoice(iblots,5))+...
%                 (discount*(1+cvht(lotchoice(iblots,5),t-1))).*Farmstats(FARMERPOS(lotchoice(iblots,2)),5,1).*...
%                 z(lotchoice(iblots,5),1));
            sampleinfo(:,6)=lotchoice(iblots,8);
            for lt=1:HT
                ils=(ibsample(ib,3)==lt);
                
                houseinfo(lt,2,nb,t)=z(lt,1);
                houseinfo(lt,3,nb,t)=z(lt,2);
%                 houseinfo(lt,7,nb,t)=z(lt,3);
                houseinfo(lt,7,nb,t)=mean(coastprox(ibsample(ibsample(:,3)==lt,2)));
                if isempty(find(ils,1))==1
                    houseinfo(lt,[1 4:6],nb,t)=0;
                else
                    houseinfo(lt,1,nb,t)=mean(sampleinfo(ils,1));
                    houseinfo(lt,4,nb,t)=mean(sampleinfo(ils,4));
                    houseinfo(lt,5,nb,t)=mean(sampleinfo(ils,5));
                    houseinfo(lt,6,nb,t)=length(find(ils==1));
                end
                if isempty(find(houseinfo(lt,5,nb,t),1))==0
                    BIDLEVEL(ibsample(ibsample(:,3)==lt,2))=houseinfo(lt,5,nb,t);
                end
                bidlevel(lt,t)=mean(con2lot(lotchoice(:,5)==lt,1)./(ccost(lt)+...
                    discount*Farmstats(FARMERPOS(lotchoice(lotchoice(:,5)==lt,2)),5,1).*z(lt,1)));
%                 bidlevel(lt,t)=mean(con2lot(lotchoice(:,5)==lt,1)./(ccost(lt)+...
%                     (discount*(1+cvht(lt,t-1)))*Farmstats(FARMERPOS(lotchoice(lotchoice(:,5)==lt,2)),5,1).*z(lt,1)));
            end
        end
        warning('on','all');

        

        bb=bbfull(:,:,nb);
        for lt=1:HT
            brokererror(lt,:,nb,t) = (1-DELTA)*brokererror(lt,:,nb,t-1)+...
                DELTA*abs(houseinfo(lt,1,nb,t)-brokerproj(lt,:,nb));
            brokerabserror(lt,:,nb,t)=houseinfo(lt,1,nb,t)-brokerproj(lt,:,nb);
            [brokerbest ibrokerbest] = min(brokererror(lt,:,nb,t),[],2);
            diffbrokererror(lt,:,nb,t)=brokererror(lt,:,nb,t)-brokererror(lt,:,nb,t-1);
            brokerbestabsSAVE(nb,lt,t)=brokerabserror(lt,ibrokerbest,nb,t);
            brokerbestdiffSAVE(nb,lt,t)=diffbrokererror(lt,ibrokerbest,nb,t);
            brokerbestSAVE(nb,lt,t) = brokerbest';
            ibrokerbestSAVE(nb,lt,t) = ibrokerbest';
            brokerprojSAVE(nb,lt,t) = brokerproj(lt,ibrokerbest,nb);
            brokermodelSAVE(nb,lt,t) = brokermodel(lt,ibrokerbest,nb);
            for i = 1:BROKERNUMCLASS
                strb = sprintf('brokerclass%d = find(brokermodel(lt,:,nb) == %d);',i,i);
                eval(strb);
            end
            if houseinfo(lt,1,nb,t) == 0
                bproj(lt,:)=0;
            else
                for i = 1:BROKERNUMCLASS
                    if i == 1
                        % mirror models
                       bproj(lt,brokerclass1) = houseinfo(lt,1,nb,t)+(1-bb...
                            (lt,brokerclass1)).*(0.5*houseinfo(lt,1,nb,t)-...
                            (houseinfo(lt,1,nb,t)-houseinfo(lt,1,nb,t-1)));
                    elseif i == 2
                        % mean model
                        for jl = 1:length(brokerclass2)
                            bproj(lt,brokerclass2(jl)) = mean(houseinfo(lt,1,nb,...
                                t:-1:(t-bb(lt,brokerclass2(jl)))));
                        end
                    elseif i == 3
                        %cycle model
                        bproj(lt,brokerclass3) = houseinfo(lt,1,nb,t-...
                            round(max(1,bb(lt,brokerclass3))));
                    elseif i == 4
                        % projection model
                        for jl = 1:length(brokerclass4)
                            %Nonlinear Forecast
                            indata=houseinfo(lt,1,nb,t-(1+bb(lt,brokerclass4(jl))):t);
                            subindata=reshape(indata,1,length(indata));
                            pcoef=polyfit(1:length(indata),subindata,1);
                            pline=pcoef(1).*(1:length(indata)+1)+pcoef(2);
%                             pcoef=polyfit(1:length(indata),subindata,2);
%                             pline=pcoef(1).*(1:length(indata)+1).^2+pcoef(2).*(1:...
%                                 length(indata)+1)+pcoef(3);
                            bproj(lt,brokerclass4(jl))=pline(length(pline));
                        end
                    elseif i == 5
                        % rescale model
                        bproj(lt,brokerclass5) = bb(lt,brokerclass5)*houseinfo(lt,1,nb,t);
                    elseif i == 6
                        [brows bcols]=ind2sub([nbrokerlong nbrokerwide],nb);
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
            
            isub=(ibsample(:,3) == lt);
            EXPTHOUSE(ibsample(isub,2))=bproj(lt,ibrokerbest);

            numlt(lt,t)=length(find(LOTTYPE==lt))/z(lt,1);   
        end
        brokerproj(:,:,nb)=bproj;
    end
    
    bcheck(:,1:Nbrokers)=houseinfo(:,1,1:Nbrokers,t);
    ibcheck=(bcheck' ~= 0);
    htexist=ismember(1:HT,LOTTYPE);
    hset=htset(htexist);
    for lt=1:HT
        ihtexist=(ibcheck(:,lt)==1);
        if numlt(lt,t) >= 2 && length(find(lotchoice(lotchoice(:,5)==lt,7)==1)) >= 2
            avgbrokervar(lt,t)=0.5*avgbrokervar(lt,t-1)+...
                0.5*mean(var(brokerbestabsSAVE(ihtexist,lt,t)));
            probloss(lt,t)=0.5*probloss(lt,t-1)+...
                0.5*mean((brokerbestabsSAVE(ihtexist,lt,t)>0));
            abserror=brokerbestabsSAVE(ihtexist,lt,t);
%             [muover,sigmaover]=normfit(abserror(abserror >= 0));
%             phatover(lt,:)=[muover sigmaover];
%             [muunder,sigmaunder]=normfit(abserror(abserror <= 0));
%             phatunder(lt,:)=[muunder sigmaunder];
            [mu,sigma]=normfit(abserror);
            phat(lt,:)=[mu sigma];
            probeven(lt,t)=cdf('norm',0,phat(lt,1),phat(lt,2));
            probover(lt,t)=length(find(abserror > 0))/length(abserror);
            probunder(lt,t)=length(find(abserror < 0))/length(abserror);
            overvalue(lt,t)=icdf('norm',max(min(probeven(lt,t)+(1-probeven(lt,t))*...
                probover(lt,t),0.99),0.01),phat(lt,1),phat(lt,2));
            undervalue(lt,t)=icdf('norm',max(min(probeven(lt,t)*...
                (1-probunder(lt,t)),0.99),0.01),phat(lt,1),phat(lt,2));
%             avgover(lt,t)=mean(abserror(abserror >= 0));
%             avgunder(lt,t)=mean(abserror(abserror <= 0));
%             probover(lt,t)=sum(abs(abserror(abserror >= 0)))/sum(abs(abserror));
%             probunder(lt,t)=sum(abs(abserror(abserror <= 0)))/sum(abs(abserror));
%             cvht(lt,t)=0.5*cvht(lt,t-1)+0.5*std(mean(brokerbestabsSAVE(ihtexist,lt,t),3))/...
            mean(con2lot(lotchoice(:,5)==lt,1));
            maxvalue(lt,t)=icdf('norm',max(min(probeven(lt,t)+...
                (1-probeven(lt,t))*0.99,0.99),0.01),phat(lt,1),phat(lt,2));
            minvalue(lt,t)=icdf('norm',max(min(probeven(lt,t)*...
                (1-0.99),0.99),0.01),phat(lt,1),phat(lt,2));
        else
            isimvar=hset(ismember(htset(numlt(:,t)>=2),min(simlotrange(htset(lt),1)):...
                max(simlotrange(htset(lt),2))));
            avgbrokervar(lt,t)=max(avgbrokervar(isimvar,t));
            probloss(lt,t)=alpha_gain/(alpha_gain+alpha_loss);
            overvalue(lt,t)=mean(maxvalue(isimvar,TSTART));
            undervalue(lt,t)=mean(minvalue(isimvar,TSTART));
            probover(lt,t)=alpha_loss/(alpha_gain+alpha_loss);
            probunder(lt,t)=alpha_gain/(alpha_gain+alpha_loss);
%             cvht(lt,t)=max(cvht(isimvar,t));
        end 
        
    end

    numtotbids(:,t)=sum(houseinfo(:,4,:,t),3);
    for lt=1:HT
        iocc=(con2lot(:,2)~=0 & lotchoice(:,5)==lt);
        htincome(lt,t)=mean(Income(con2lot(iocc,2)));
    end
    
    consumerstats(1,t)=Nconsumers;
    consumerstats(4,t)=mean(housemp);
    consumerstats(2,t)=mean(Income(con2lot(ifilled,2)));
    consumerstats(3,t)=mean(AVGUTIL(iurb));
    consumerstats(5,t)=mean(OutIncome(:,1));
    Farmstats(iNfarmers,:,t+1)=Farmstats(iNfarmers,:,t);
    agrland(t)=length(find(BASELAYER==0));
    
    
    if isempty(find(t == testtime,1))==0
        tmark=find(t == testtime);
        BT(:,:,tmark)=BUILDTIME;
        BL(:,:,tmark)=BIDLEVEL;
        Rvacland(:,:,tmark)=VACLAND;
        Rrent(:,:,tmark)=RENT;
        Rlottype(:,:,tmark)=LOTTYPE;  
        Rbaselayer(:,:,tmark)=BASELAYER;
        subprefmap=zeros(NLENGTH,NWIDTH);
        subriskmap=zeros(NLENGTH,NWIDTH);
        for ires=1:length(con2lot(ifilled,1))
%             hopt=((Income(con2lot(ifilled(ires),2),1)-travelcost(lotchoice(ifilled(ires),2))-...
%                 avgrent(:,t)).^CONINFO(con2lot(ifilled(ires),2),1)).*(z(:,2).^...
%                 CONINFO(con2lot(ifilled(ires),2),2)).*(z(:,1).^CONINFO(con2lot(ifilled(ires),2),3)).*...
%                 (z(lt,3).^CONINFO(con2lot(ifilled(ires),2),4));
            hopt=((Income(con2lot(ifilled(ires),2),1)-travelcost(lotchoice(ifilled(ires),2))-...
                avgrent(:,t)).^CONINFO(con2lot(ifilled(ires),2),1)).*(z(:,2).^...
                CONINFO(con2lot(ifilled(ires),2),2)).*(z(:,1).^CONINFO(con2lot(ifilled(ires),2),3)).*...
                (coastprox(lotchoice(ifilled(ires),2)).^CONINFO(con2lot(ifilled(ires),2),4));
            [imaxu,jmaxu]=max(hopt,[],1);
            idealset(jmaxu,tmark)=idealset(jmaxu,tmark)+1;
        end
        
        for i=1:length(ifilled)
            profopt=(avgrent(:,t)-ones(HT,1)*PLAND(lotchoice(ifilled(i),2)).*...
                z(:,1).*discount-ccost)./z(:,1);
%             profopt=(avgrent(:,t)-ones(HT,1)*PLAND(lotchoice(ifilled(i),2)).*...
%                 z(:,1).*(discount*(1+cvht(:,t-1)))-ccost)./z(:,1);
            [imaxp,jmaxp]=max(profopt,[],1);
            profset(jmaxp,tmark)=profset(jmaxp,tmark)+1;
            
            subprefmap(lotchoice(ifilled(i),2))=mean(CONINFO(con2lot(ifilled(i),2),4));
            subriskmap(lotchoice(ifilled(i),2))=mean(CONINFO(con2lot(ifilled(i),2),5));
        end
        PREFMAP(:,:,tmark)=subprefmap;
        SUBRISKMAP(:,:,tmark)=subriskmap;
        for lt=1:HT
            iret=find(LOTTYPE==lt);
            if isempty(iret)==1
                continue
            end
            subRreturn(iret)=(RENT(iret)-((PLAND(iret).*z(lt,1).*...
                discount)+ccost(lt)))./z(lt,1);
%             subRreturn(iret)=(RENT(iret)-((PLAND(iret).*z(lt,1).*...
%                 (discount*(1+cvht(lt,t-1))))+ccost(lt)))./z(lt,1);
        end
        Rreturn(:,:,tmark)=subRreturn;   

        for ii=1:length(con2lot(:,3))
            isamelot=find(Lottype(:,1)==lotchoice(ii,1));
            isamecell=find(lotchoice(:,2)==lotchoice(ii,2));
            if length(isamelot)==1
                if length(isamecell) == 1 && con2lot(isamecell,2)~= 0
                    AVGINCOME(lotchoice(isamecell,2))=Income(con2lot(ii,2));
                elseif length(isamecell) > 1 && isempty(find(con2lot(isamecell,2)~=0,1))==0
                    AVGINCOME(lotchoice(isamecell,2))=mean(Income(con2lot...
                        (isamecell(con2lot(isamecell,2)~=0),2)));
                elseif isempty(find(con2lot(isamecell,2)~=0,1)) == 1
                    AVGINCOME(lotchoice(isamecell,2))=0;
                end
            elseif length(isamelot) > 1
                if con2lot(isamecell,2) ~= 0
                    AVGINCOME(Lottype(isamelot,2))=Income(con2lot(isamecell,2));
                elseif con2lot(isamecell,2) == 0
                    AVGINCOME(lotchoice(isamecell,2))=0;
                end                
            end           
        end
        Rincome(:,:,tmark)=AVGINCOME;
        Rpland(:,:,tmark)=PLAND;
        Rpop(tmark)=length(ifilled);
        Rvacrate(tmark)=vacrate(t);
        Rvaclots(tmark)=vacantlots(t);
        Rnumlot(tmark)=numlt(t);
        Rleftoverpop(tmark)=leftoverpop(t);
    end
    
    ibuildtime=find(BUILDTIME==t);
    iltype=ismember(Lottype(:,2),ibuildtime);
    subdynlt(Lottype(iltype,2))=Lottype(iltype,5);
    subdynrent(Lottype(iltype,2))=con2lot(Lottype(iltype,1),1);
    dynbaselayer=(subdynlt~=0);
    Dynltmap(:,:,t)=subdynlt;
    Dynrentmap(:,:,t)=subdynrent;
    
    for lt=1:HT
        htperyear(lt,t)=length(find(Newlottype(:,:,t)==lt))/z(lt,1);
    end
    alldevdist=10000*ones(NLENGTH,NWIDTH);
    transland=find(Dynltmap(:,:,t)~=0);
    for tl=1:length(transland)
        tdist2hznnei=10000*ones(NLENGTH,NWIDTH);    %distance to horizontal neighbor from icenter
        tdist2vrtnei=10000*ones(NLENGTH,NWIDTH);
        [itlrow,itlcol]=ind2sub([NLENGTH,NWIDTH],transland(tl));
        for col=1:NWIDTH
              tdist2hznnei(1:NLENGTH,col)=abs(col-itlcol).*...
                ones(NLENGTH,1);
        end
        for row=1:NWIDTH
              tdist2vrtnei(row,1:NWIDTH)=abs(row-itlrow).*...
                ones(1,NLENGTH);
        end
        for col=1:NWIDTH
            for row=1:NLENGTH
                alldevdist(row,col)=min(sqrt(tdist2hznnei(row,col)^2+...
                    tdist2vrtnei(row,col)^2),alldevdist(row,col));
            end
        end
    end
    for nf=1:Nfarmers
        Farmdist2dev(nf,t)=mean(alldevdist(FARMERPOS==nf));
    end
    subbmodel=zeros(NLENGTH,NWIDTH);
    subbproj=zeros(NLENGTH,NWIDTH);
    for nb=1:Nbrokers
        for lt=1:HT
            totbrokerrecord(:,lt,t,nb)=[houseinfo(lt,1,nb,t); brokerprojSAVE(nb,lt,t); ...
                brokermodelSAVE(nb,lt,t)];
            ibarea=find(HBROKER==nb);
            idynlt=find(subdynlt==lt);
            itarget=ismember(idynlt,ibarea);
            subbmodel(idynlt(itarget))=brokermodelSAVE(nb,lt,t);
            subbproj(idynlt(itarget))=brokerprojSAVE(nb,lt,t);
            
        end
    end
    Bmodelmap(:,:,t)=subbmodel;
    Bprojmap(:,:,t)=subbproj;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Farm stats maticies
% 1. Farmstats [Farmid acres prod cost return]
% 2. farmsoldinfo [farmid selltime wtaland price epsilon wtpland dist
% pctgain]
%
% 3. totfarmrecord [wtaland wtpland landprojSAVE landmodelSAVE]
% 4. farmacreinfo [Dist ZONED? PROFHT PROFRET PROJHT PROJRET BUILDTIME]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

farmacreinfo=zeros(max(Farmstats(:,2,1)),7,Nfarmers);
farmsoldinfo=zeros([],8);
ifarmsold=find(sellrecord ~= 0);
avgfarmdist=zeros(length(ifarmsold),1);
sellepsilon=zeros(length(ifarmsold),1);
devsellwtp=zeros(length(ifarmsold),1);


for n=1:length(ifarmsold)
    avgfarmdist(n,1)=mean(DISTANCE(FARMERPOS==ifarmsold(n)));
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
WTPlandmap=zeros(80,80,TMAX);
WTAlandmap=zeros(80,80,TMAX);
Landprojmap=zeros(80,80,TMAX);
Landmodmap=zeros(80,80,TMAX);
subwtpmap=zeros(80,80);
subwtamap=zeros(80,80);
sublandprojmap=zeros(80,80);
sublandmodmap=zeros(80,80);

for ts=11:TMAX
    for nf=1:Nfarmers
        Allfarmdist(nf,1)=mean(DISTANCE(FARMERPOS==nf));
        totfarmrecord(:,ts,nf)=[wtaland(nf,ts); wtpland(nf,ts); ...
            landprojSAVE(nf,ts); landmodelSAVE(nf,ts)];
        subwtamap(FARMERPOS==nf)=totfarmrecord(1,ts,nf);
        subwtpmap(FARMERPOS==nf)=totfarmrecord(2,ts,nf);
        sublandprojmap(FARMERPOS==nf)=totfarmrecord(3,ts,nf);
        sublandmodmap(FARMERPOS==nf)=totfarmrecord(4,ts,nf);
    end
    WTAlandmap(:,:,ts)=subwtamap;
    WTPlandmap(:,:,ts)=subwtpmap;
    Landprojmap(:,:,ts)=sublandprojmap;
    Landmodmap(:,:,ts)=sublandmodmap;
    testmaxret=Dynmaxretmap(:,:,ts);
    testmaxret(Dynltmap(:,:,ts)~=0)=0;
    isold=find(totfarmrecord(3,ts,:)==0);
    testmaxret(ismember(setupmap,isold))=0;
    Dynmaxretmap(:,:,ts)=testmaxret;
    
end
avgfarmsize=mean(Farmstats(:,2,1));
stdfarmsize=std(Farmstats(:,2,1));

    
for ttt=1:TMAX
    ibuildtime=find(BUILDTIME==t);
    iltype=ismember(Lottype(:,2),ibuildtime);
    subdynlt(Lottype(iltype,2))=Lottype(iltype,5);
    subdynrent(Lottype(iltype,2))=con2lot(Lottype(iltype,1),1);
    dynbaselayer=(subdynlt~=0);
    Dynltmap(:,:,ttt)=subdynlt;
    Dynrentmap(:,:,ttt)=subdynrent;
    
    for lt=1:HT
        htperyear(lt,ttt)=length(find(Newlottype(:,:,ttt)==lt))/z(lt,1);
    end
    alldevdist=10000*ones(NLENGTH,NWIDTH);
    transland=find(Dynltmap(:,:,ttt)~=0);
    for tl=1:length(transland)
        tdist2hznnei=10000*ones(NLENGTH,NWIDTH);    %distance to horizontal neighbor from icenter
        tdist2vrtnei=10000*ones(NLENGTH,NWIDTH);
        [itlrow,itlcol]=ind2sub([NLENGTH,NWIDTH],transland(tl));
        for col=1:NWIDTH
              tdist2hznnei(1:NLENGTH,col)=abs(col-itlcol).*...
                ones(NLENGTH,1);
        end
        for row=1:NWIDTH
              tdist2vrtnei(row,1:NWIDTH)=abs(row-itlrow).*...
                ones(1,NLENGTH);
        end
        for icol=1:NWIDTH
            for irow=1:NLENGTH
                alldevdist(irow,icol)=min(sqrt(tdist2hznnei(irow,icol)^2+...
                    tdist2vrtnei(irow,icol)^2),alldevdist(irow,icol));
            end
        end
    end
    for nf=1:Nfarmers
        Farmdist2dev(nf,ttt)=mean(alldevdist(FARMERPOS==nf));
    end
    for nb=1:Nbrokers
        for lt=1:HT
            totbrokerrecord(:,lt,ttt,nb)=[houseinfo(lt,1,nb,ttt); brokerprojSAVE(nb,lt,ttt); ...
                brokermodelSAVE(nb,lt,ttt)];
        end
    end
end

testset1=zeros(1,Nbrokers);
testset1pre=zeros(1,Nbrokers);
testset2=zeros(1,Nbrokers);
testset2pre=zeros(1,Nbrokers);
Diffset=zeros(HT,TMAX);
Diffrentset=zeros(HT,TMAX);
Bprojlevel=zeros(HT,TMAX);
Rentlevel=zeros(HT,TMAX);
Bprojleveldiff=zeros(HT,TMAX);
Bprojrentdiff=zeros(HT,TMAX);
Ltperyear=zeros(HT,TMAX);
testrentset=zeros(HT,TMAX);
zpost=[0.25; 0.25; 0.25; 0.5; 0.5; 0.5; 1; 1; 1; 2; 2; 2; 5; 5; 5; 10; 10; 10];
for lt=1:HT
    Ltperyear(lt,10)=length(find(Rlottype(:,:,1)==lt))/zpost(lt);
end

for t=11:TMAX
    for lt=1:HT
        Ltperyear(lt,t)=Ltperyear(lt,10)+length(find(Dynltmap(:,:,t)==lt))/zpost(lt);
        
        testset1(1:Nbrokers)=totbrokerrecord(1,lt,t,:);
        testset1pre(1:Nbrokers)=totbrokerrecord(1,lt,t-1,:);
        testset2(1:Nbrokers)=totbrokerrecord(2,lt,t,:);  %projection for rent(t+1)
        testset2pre(1:Nbrokers)=totbrokerrecord(2,lt,t-1,:); %Related to rent(t), equals Paskhouse(t)
        ibelot=(testset1 ~= 0);
        ibelotpre=(testset1pre ~=0);
        Bprojlevel(lt,t)=mean(testset2(ibelot));
        Rentlevel(lt,t)=mean(testset1(ibelot));
        Bprojleveldiff(lt,t)=mean(testset2(ibelot))-mean(testset2pre(ibelotpre));
        
        Bprojrentdiff(lt,t)=mean(testset2pre(ibelotpre)-testset1pre(ibelotpre));
        testrentset(lt,t)=mean(testset1(ibelotpre)-testset1pre(ibelotpre));
        Diffset(lt,t)=mean(testset2(ibelot))-mean(testset1pre(ibelotpre));
        Diffrentset(lt,t)=mean(testset1(ibelot))-...
            mean(testset1pre(ibelotpre));
    end
end

for nf=1:Nfarmers
    iacres=find(FARMERPOS == nf);
    [nfrow,nfcol]=ind2sub([NLENGTH NWIDTH],iacres);
    bts=unique(BUILDTIME(iacres));
    
    if bts == 0
        continue
    end
    
    bts=bts(bts~=0);
    
    farmacreinfo(1:length(iacres),1,nf)=DISTANCE(iacres);
    farmacreinfo(1:length(iacres),2,nf)=ZONED(iacres);
    
    profltbtmap=zeros(size(BUILDTIME));
    profretbtmap=zeros(size(BUILDTIME));
    projltbtmap=zeros(size(BUILDTIME));
    projretbtmap=zeros(size(BUILDTIME));
    for bt=1:length(bts)
        subprofltmap=Dynltmap(:,:,bts(bt));
        subprofretmap=Dynmaxprofmap(:,:,bts(bt));
        subprojltmap=Dynretltmap(:,:,bts(bt));
        subprojretmap=Dynmaxretmap(:,:,bts(bt));
        
        ibt=(BUILDTIME == bts(bt));
        profltbtmap(ibt)=subprofltmap(ibt);
        profretbtmap(ibt)=subprofretmap(ibt);
        projltbtmap(ibt)=subprojltmap(ibt);
        projretbtmap(ibt)=subprojretmap(ibt);
    end
    farmacreinfo(1:length(iacres),3,nf)=profltbtmap(iacres);
    farmacreinfo(1:length(iacres),4,nf)=profretbtmap(iacres);
    farmacreinfo(1:length(iacres),5,nf)=projltbtmap(iacres);
    farmacreinfo(1:length(iacres),6,nf)=profretbtmap(iacres);
    farmacreinfo(1:length(iacres),7,nf)=BUILDTIME(iacres);
end


%%% PLOTTING COMMANDS %%%

% @@@@@@@@@@@ plot 5x1 htperyear, 4*diff b/w broj(t)-rent(t); @@@@@@
% rent(t)-rent(t-1) <><><><>%
% %@@@@@ by house size @@@@@@@
% 
%%%% Small Houses %%%%
% figure
% subplot(4,2,1)
% plot(TSTART:TMAX-1,htperyear([1 4 7 10],TSTART:TMAX-1))
% title('Number of Small Houses Built')
% legend('HT 1','HT 4','HT 7','HT 10','Best')
% subplot(4,2,3)
% plot(TSTART:TMAX-1,avgrent([1 4 7 10],TSTART:TMAX-1))
% title('Average Rent For Small Houses')
% subplot(4,2,5)
% plot(TSTART:TMAX-1,htincome([1 4 7 10],TSTART:TMAX-1))
% title('Average Income For Small Houses')
% subplot(4,2,7)
% [row,col]=find(bidshare > 0);
% plotrow=unique(row);
% area(bidshare(plotrow,11:30)')
% axis([1 20 0 1])
% colorbar
% title('Share of Bids Per Housing Type')
% xlabel('Time Step')
% subplot(4,2,2)
% %%% Bproj and Rent levels
% % plot(TSTART:TMAX-1,Bprojlevel(1,TSTART:TMAX-1))
% % hold on
% % plot(TSTART:TMAX-1,Rentlevel(1,TSTART:TMAX-1),'r')
% % title('Projected Rent(t)(blue) and Actual Rent(t)(red) for HT 1')    
% 
% % %%% Change in Past Rent Effect on Change in Present Bproj %%%%
% % % plot(TSTART:TMAX-1,Diffset(1,TSTART:TMAX-1))
% % line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% % hold on
% % plot(TSTART+2:TMAX-1,Bprojleveldiff(1,TSTART+3:TMAX),'k')
% % plot(TSTART+2:TMAX-1,Diffrentset(1,TSTART+2:TMAX-1),'r')
% % title('Change between Projected Rent(t+1)-(t)(black) and Rent(t)-(t-1)(red) for HT 1')    
% 
% % %%% Change in Past Bproj Effect on Change in Present Rent %%%%
% % plot(TSTART:TMAX-1,Diffset(1,TSTART:TMAX-1))
% line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% hold on
% plot(TSTART+2:TMAX-1,Bprojleveldiff(1,TSTART+2:TMAX-1),'k')
% plot(TSTART+2:TMAX-1,Diffrentset(1,TSTART+3:TMAX),'r')
% title('Change between Rent(t+1)-(t)(red) and Projected Rent(t)-(t-1)(black) for HT 1')
% 
% subplot(4,2,4)
% %%% Bproj and Rent levels
% % plot(TSTART:TMAX-1,Bprojlevel(4,TSTART:TMAX-1))
% % hold on
% % plot(TSTART:TMAX-1,Rentlevel(4,TSTART:TMAX-1),'r')
% % title('Projected Rent(t)(blue) and Actual Rent(t)(red) for HT 4')    
% 
% % %%% Change in Past Rent Effect on Change in Present Bproj %%%%
% % % plot(TSTART:TMAX-1,Diffset(4,TSTART:TMAX-1))
% % line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% % hold on
% % plot(TSTART+2:TMAX-1,Bprojleveldiff(4,TSTART+3:TMAX),'k')
% % plot(TSTART+2:TMAX-1,Diffrentset(4,TSTART+2:TMAX-1),'r')
% % title('Change between Projected Rent(t+1)-(t)(black) and Rent(t)-(t-1)(red) for HT 4')    
% 
% % %%% Change in Past Bproj Effect on Change in Present Rent %%%%
% % plot(TSTART:TMAX-1,Diffset(4,TSTART:TMAX-1))
% line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% hold on
% plot(TSTART+2:TMAX-1,Bprojleveldiff(4,TSTART+2:TMAX-1),'k')
% plot(TSTART+2:TMAX-1,Diffrentset(4,TSTART+3:TMAX),'r')
% title('Change between Rent(t+1)-(t)(red) and Projected Rent(t)-(t-1)(black) for HT 4')
% 
% subplot(4,2,6)
% %%% Bproj and Rent levels
% % plot(TSTART:TMAX-1,Bprojlevel(7,TSTART:TMAX-1))
% % hold on
% % plot(TSTART:TMAX-1,Rentlevel(7,TSTART:TMAX-1),'r')
% % title('Projected Rent(t)(blue) and Actual Rent(t)(red) for HT 7')    
% 
% % %%% Change in Past Rent Effect on Change in Present Bproj %%%%
% % % plot(TSTART:TMAX-1,Diffset(7,TSTART:TMAX-1))
% % line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% % hold on
% % plot(TSTART+2:TMAX-1,Bprojleveldiff(7,TSTART+3:TMAX),'k')
% % plot(TSTART+2:TMAX-1,Diffrentset(7,TSTART+2:TMAX-1),'r')
% % title('Change between Projected Rent(t+1)-(t)(black) and Rent(t)-(t-1)(red) for HT 7')    
% 
% % %%% Change in Past Bproj Effect on Change in Present Rent %%%%
% % plot(TSTART:TMAX-1,Diffset(7,TSTART:TMAX-1))
% line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% hold on
% plot(TSTART+2:TMAX-1,Bprojleveldiff(7,TSTART+2:TMAX-1),'k')
% plot(TSTART+2:TMAX-1,Diffrentset(7,TSTART+3:TMAX),'r')
% title('Change between Rent(t+1)-(t)(red) and Projected Rent(t)-(t-1)(black) for HT 7')
% 
% subplot(4,2,8)
% %%% Bproj and Rent levels
% % plot(TSTART:TMAX-1,Bprojlevel(10,TSTART:TMAX-1))
% % hold on
% % plot(TSTART:TMAX-1,Rentlevel(10,TSTART:TMAX-1),'r')
% % title('Projected Rent(t)(blue) and Actual Rent(t)(red) for HT 10')    
% 
% % %%% Change in Past Rent Effect on Change in Present Bproj %%%%
% % % plot(TSTART:TMAX-1,Diffset(10,TSTART:TMAX-1))
% % line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% % hold on
% % plot(TSTART+2:TMAX-1,Bprojleveldiff(10,TSTART+3:TMAX),'k')
% % plot(TSTART+2:TMAX-1,Diffrentset(10,TSTART+2:TMAX-1),'r')
% % title('Change between Projected Rent(t+1)-(t)(black) and Rent(t)-(t-1)(red) for HT 10')    
% 
% % %%% Change in Past Bproj Effect on Change in Present Rent %%%%
% % plot(TSTART:TMAX-1,Diffset(10,TSTART:TMAX-1))
% line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% hold on
% plot(TSTART+2:TMAX-1,Bprojleveldiff(10,TSTART+2:TMAX-1),'k')
% plot(TSTART+2:TMAX-1,Diffrentset(10,TSTART+3:TMAX),'r')
% title('Change between Rent(t+1)-(t)(red) and Projected Rent(t)-(t-1)(black) for HT 10')
%  

%%%% Medium Houses %%%%
% figure
% subplot(4,2,1)
% plot(TSTART:TMAX-1,htperyear([2 5 8 11],TSTART:TMAX-1))
% title('Number of Medium Houses Built')
% legend('HT 2','HT 5','HT 8','HT 11','Best')
% subplot(4,2,3)
% plot(TSTART:TMAX-1,avgrent([2 5 8 11],TSTART:TMAX-1))
% title('Average Rent For Medium Houses')
% subplot(4,2,5)
% plot(TSTART:TMAX-1,htincome([2 5 8 11],TSTART:TMAX-1))
% title('Average Income For Medium Houses')
% subplot(4,2,7)
% [row,col]=find(bidshare > 0);
% plotrow=unique(row);
% area(bidshare(plotrow,11:30)')
% axis([1 20 0 1])
% colorbar
% title('Share of Bids Per Housing Type')
% xlabel('Time Step')
% subplot(4,2,2)
% %%% Bproj and Rent levels
% % plot(TSTART:TMAX-1,Bprojlevel(2,TSTART:TMAX-1))
% % hold on
% % plot(TSTART:TMAX-1,Rentlevel(2,TSTART:TMAX-1),'r')
% % title('Projected Rent(t)(blue) and Actual Rent(t)(red) for HT 2')    
% 
% % %%% Change in Past Rent Effect on Change in Present Bproj %%%%
% % % plot(TSTART:TMAX-1,Diffset(2,TSTART:TMAX-1))
% % line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% % hold on
% % plot(TSTART+2:TMAX-1,Bprojleveldiff(2,TSTART+3:TMAX),'k')
% % plot(TSTART+2:TMAX-1,Diffrentset(2,TSTART+2:TMAX-1),'r')
% % title('Change between Projected Rent(t+1)-(t)(black) and Rent(t)-(t-1)(red) for HT 2')    
% 
% % %%% Change in Past Bproj Effect on Change in Present Rent %%%%
% % plot(TSTART:TMAX-1,Diffset(2,TSTART:TMAX-1))
% line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% hold on
% plot(TSTART+2:TMAX-1,Bprojleveldiff(2,TSTART+2:TMAX-1),'k')
% plot(TSTART+2:TMAX-1,Diffrentset(2,TSTART+3:TMAX),'r')
% title('Change between Rent(t+1)-(t)(red) and Projected Rent(t)-(t-1)(black) for HT 2')
% 
% subplot(4,2,4)
% %%% Bproj and Rent levels
% % plot(TSTART:TMAX-1,Bprojlevel(5,TSTART:TMAX-1))
% % hold on
% % plot(TSTART:TMAX-1,Rentlevel(5,TSTART:TMAX-1),'r')
% % title('Projected Rent(t)(blue) and Actual Rent(t)(red) for HT 5')    
% 
% % %%% Change in Past Rent Effect on Change in Present Bproj %%%%
% % % plot(TSTART:TMAX-1,Diffset(5,TSTART:TMAX-1))
% % line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% % hold on
% % plot(TSTART+2:TMAX-1,Bprojleveldiff(5,TSTART+3:TMAX),'k')
% % plot(TSTART+2:TMAX-1,Diffrentset(5,TSTART+2:TMAX-1),'r')
% % title('Change between Projected Rent(t+1)-(t)(black) and Rent(t)-(t-1)(red) for HT 5')    
% 
% % %%% Change in Past Bproj Effect on Change in Present Rent %%%%
% % plot(TSTART:TMAX-1,Diffset(5,TSTART:TMAX-1))
% line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% hold on
% plot(TSTART+2:TMAX-1,Bprojleveldiff(5,TSTART+2:TMAX-1),'k')
% plot(TSTART+2:TMAX-1,Diffrentset(5,TSTART+3:TMAX),'r')
% title('Change between Rent(t+1)-(t)(red) and Projected Rent(t)-(t-1)(black) for HT 5')
% 
% subplot(4,2,6)
% %%% Bproj and Rent levels
% % plot(TSTART:TMAX-1,Bprojlevel(8,TSTART:TMAX-1))
% % hold on
% % plot(TSTART:TMAX-1,Rentlevel(8,TSTART:TMAX-1),'r')
% % title('Projected Rent(t)(blue) and Actual Rent(t)(red) for HT 8')    
% 
% % %%% Change in Past Rent Effect on Change in Present Bproj %%%%
% % % plot(TSTART:TMAX-1,Diffset(8,TSTART:TMAX-1))
% % line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% % hold on
% % plot(TSTART+2:TMAX-1,Bprojleveldiff(8,TSTART+3:TMAX),'k')
% % plot(TSTART+2:TMAX-1,Diffrentset(8,TSTART+2:TMAX-1),'r')
% % title('Change between Projected Rent(t+1)-(t)(black) and Rent(t)-(t-1)(red) for HT 8')    
% 
% % %%% Change in Past Bproj Effect on Change in Present Rent %%%%
% % plot(TSTART:TMAX-1,Diffset(8,TSTART:TMAX-1))
% line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% hold on
% plot(TSTART+2:TMAX-1,Bprojleveldiff(8,TSTART+2:TMAX-1),'k')
% plot(TSTART+2:TMAX-1,Diffrentset(8,TSTART+3:TMAX),'r')
% title('Change between Rent(t+1)-(t)(red) and Projected Rent(t)-(t-1)(black) for HT 8')
% 
% subplot(4,2,8)
% %%% Bproj and Rent levels
% % plot(TSTART:TMAX-1,Bprojlevel(11,TSTART:TMAX-1))
% % hold on
% % plot(TSTART:TMAX-1,Rentlevel(11,TSTART:TMAX-1),'r')
% % title('Projected Rent(t)(blue) and Actual Rent(t)(red) for HT 11')    
% 
% % %%% Change in Past Rent Effect on Change in Present Bproj %%%%
% % % plot(TSTART:TMAX-1,Diffset(11,TSTART:TMAX-1))
% % line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% % hold on
% % plot(TSTART+2:TMAX-1,Bprojleveldiff(11,TSTART+3:TMAX),'k')
% % plot(TSTART+2:TMAX-1,Diffrentset(11,TSTART+2:TMAX-1),'r')
% % title('Change between Projected Rent(t+1)-(t)(black) and Rent(t)-(t-1)(red) for HT 11')    
% 
% % %%% Change in Past Bproj Effect on Change in Present Rent %%%%
% % plot(TSTART:TMAX-1,Diffset(11,TSTART:TMAX-1))
% line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% hold on
% plot(TSTART+2:TMAX-1,Bprojleveldiff(11,TSTART+2:TMAX-1),'k')
% plot(TSTART+2:TMAX-1,Diffrentset(11,TSTART+3:TMAX),'r')
% title('Change between Rent(t+1)-(t)(red) and Projected Rent(t)-(t-1)(black) for HT 11')


%%%% Large Houses %%%%
% figure
% subplot(4,2,1)
% plot(TSTART:TMAX-1,htperyear([3 6 9 12],TSTART:TMAX-1))
% title('Number of Large Houses Built')
% legend('HT 3','HT 6','HT 9','HT 12','Best')
% subplot(4,2,3)
% plot(TSTART:TMAX-1,avgrent([3 6 9 12],TSTART:TMAX-1))
% title('Average Rent For Large Houses')
% subplot(4,2,5)
% plot(TSTART:TMAX-1,htincome([3 6 9 12],TSTART:TMAX-1))
% title('Average Income For Large Houses')
% subplot(4,2,7)
% [row,col]=find(bidshare > 0);
% plotrow=unique(row);
% area(bidshare(plotrow,11:30)')
% axis([1 20 0 1])
% colorbar
% title('Share of Bids Per Housing Type')
% xlabel('Time Step')
% subplot(4,2,2)
% %%% Bproj and Rent levels
% % plot(TSTART:TMAX-1,Bprojlevel(3,TSTART:TMAX-1))
% % hold on
% % plot(TSTART:TMAX-1,Rentlevel(3,TSTART:TMAX-1),'r')
% % title('Projected Rent(t)(blue) and Actual Rent(t)(red) for HT 3')    
% 
% % %%% Change in Past Rent Effect on Change in Present Bproj %%%%
% % % plot(TSTART:TMAX-1,Diffset(3,TSTART:TMAX-1))
% % line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% % hold on
% % plot(TSTART+2:TMAX-1,Bprojleveldiff(3,TSTART+3:TMAX),'k')
% % plot(TSTART+2:TMAX-1,Diffrentset(3,TSTART+2:TMAX-1),'r')
% % title('Change between Projected Rent(t+1)-(t)(black) and Rent(t)-(t-1)(red) for HT 3')    
% 
% % %%% Change in Past Bproj Effect on Change in Present Rent %%%%
% % plot(TSTART:TMAX-1,Diffset(3,TSTART:TMAX-1))
% line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% hold on
% plot(TSTART+2:TMAX-1,Bprojleveldiff(3,TSTART+2:TMAX-1),'k')
% plot(TSTART+2:TMAX-1,Diffrentset(3,TSTART+3:TMAX),'r')
% title('Change between Rent(t+1)-(t)(red) and Projected Rent(t)-(t-1)(black) for HT 3')
% 
% subplot(4,2,4)
% %%% Bproj and Rent levels
% % plot(TSTART:TMAX-1,Bprojlevel(6,TSTART:TMAX-1))
% % hold on
% % plot(TSTART:TMAX-1,Rentlevel(6,TSTART:TMAX-1),'r')
% % title('Projected Rent(t)(blue) and Actual Rent(t)(red) for HT 6')    
% 
% % %%% Change in Past Rent Effect on Change in Present Bproj %%%%
% % % plot(TSTART:TMAX-1,Diffset(6,TSTART:TMAX-1))
% % line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% % hold on
% % plot(TSTART+2:TMAX-1,Bprojleveldiff(6,TSTART+3:TMAX),'k')
% % plot(TSTART+2:TMAX-1,Diffrentset(6,TSTART+2:TMAX-1),'r')
% % title('Change between Projected Rent(t+1)-(t)(black) and Rent(t)-(t-1)(red) for HT 6')    
% 
% % %%% Change in Past Bproj Effect on Change in Present Rent %%%%
% % plot(TSTART:TMAX-1,Diffset(3,TSTART:TMAX-1))
% line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% hold on
% plot(TSTART+2:TMAX-1,Bprojleveldiff(6,TSTART+2:TMAX-1),'k')
% plot(TSTART+2:TMAX-1,Diffrentset(6,TSTART+3:TMAX),'r')
% title('Change between Rent(t+1)-(t)(red) and Projected Rent(t)-(t-1)(black) for HT 6')
% 
% subplot(4,2,6)
% %%% Bproj and Rent levels
% % plot(TSTART:TMAX-1,Bprojlevel(9,TSTART:TMAX-1))
% % hold on
% % plot(TSTART:TMAX-1,Rentlevel(9,TSTART:TMAX-1),'r')
% % title('Projected Rent(t)(blue) and Actual Rent(t)(red) for HT 9')    
% 
% % %%% Change in Past Rent Effect on Change in Present Bproj %%%%
% % % plot(TSTART:TMAX-1,Diffset(9,TSTART:TMAX-1))
% % line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% % hold on
% % plot(TSTART+2:TMAX-1,Bprojleveldiff(9,TSTART+3:TMAX),'k')
% % plot(TSTART+2:TMAX-1,Diffrentset(9,TSTART+2:TMAX-1),'r')
% % title('Change between Projected Rent(t+1)-(t)(black) and Rent(t)-(t-1)(red) for HT 9')    
% 
% % %%% Change in Past Bproj Effect on Change in Present Rent %%%%
% % plot(TSTART:TMAX-1,Diffset(9,TSTART:TMAX-1))
% line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% hold on
% plot(TSTART+2:TMAX-1,Bprojleveldiff(9,TSTART+2:TMAX-1),'k')
% plot(TSTART+2:TMAX-1,Diffrentset(9,TSTART+3:TMAX),'r')
% title('Change between Rent(t+1)-(t)(red) and Projected Rent(t)-(t-1)(black) for HT 9')
% 
% subplot(4,2,8)
% %%% Bproj and Rent levels
% % plot(TSTART:TMAX-1,Bprojlevel(12,TSTART:TMAX-1))
% % hold on
% % plot(TSTART:TMAX-1,Rentlevel(12,TSTART:TMAX-1),'r')
% % title('Projected Rent(t)(blue) and Actual Rent(t)(red) for HT 12')    
% 
% % %%% Change in Past Rent Effect on Change in Present Bproj %%%%
% % % plot(TSTART:TMAX-1,Diffset(12,TSTART:TMAX-1))
% % line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% % hold on
% % plot(TSTART+2:TMAX-1,Bprojleveldiff(12,TSTART+3:TMAX),'k')
% % plot(TSTART+2:TMAX-1,Diffrentset(12,TSTART+2:TMAX-1),'r')
% % title('Change between Projected Rent(t+1)-(t)(black) and Rent(t)-(t-1)(red) for HT 12')    
% 
% % %%% Change in Past Bproj Effect on Change in Present Rent %%%%
% % plot(TSTART:TMAX-1,Diffset(12,TSTART:TMAX-1))
% line(TSTART+2:TMAX-1,zeros(1,length(TSTART+2:TMAX-1)))
% hold on
% plot(TSTART+2:TMAX-1,Bprojleveldiff(12,TSTART+2:TMAX-1),'k')
% plot(TSTART+2:TMAX-1,Diffrentset(12,TSTART+3:TMAX),'r')
% title('Change between Rent(t+1)-(t)(red) and Projected Rent(t)-(t-1)(black) for HT 12')


% 
% plot(TSTART+1:TMAX-1,Avgexptdiff(10,TSTART+1:TMAX-1)./max(Avgexptdiff(10,TSTART+1:TMAX-1)))
% hold on
% plot(TSTART+1:TMAX-1,htperyear(10,TSTART+1:TMAX-1)./max(htperyear(10,TSTART+1:TMAX-1)),'k--')
% plot(TSTART+1:TMAX-1,Bprojlevel(10,TSTART+1:TMAX-1)./max(Bprojlevel(10,TSTART+1:TMAX-1)),'c')
% plot(TSTART+1:TMAX-1,Rentlevel(10,TSTART+1:TMAX-1)./max(Rentlevel(10,TSTART+1:TMAX-1)),'r')
% plot(TSTART+1:TMAX-1,bidshare(10,TSTART+1:TMAX-1),'g')
% legend('Avgexptdiff','htperyear','Bprojlevel','Rentlevel','Bidshare')
% 
% 
% plot(TSTART+1:TMAX-1,htperyear(10,TSTART+1:TMAX-1)./max(htperyear(10,TSTART+1:TMAX-1)),'k--')
% hold on
% plot(TSTART+1:TMAX-1,Avgnewbid(10,TSTART+1:TMAX-1)-1)
% plot(TSTART+1:TMAX-1,Rentlevel(10,TSTART+1:TMAX-1)./max(Rentlevel(10,TSTART+1:TMAX-1)),'r')
% plot(TSTART+1:TMAX-1,bidshare(10,TSTART+1:TMAX-1),'g')
% plot(TSTART+1:TMAX-1,Realavgret(10,TSTART+1:TMAX-1)./max(Realavgret(10,TSTART+1:TMAX-1)),'m')
% legend('htperyear','Avgnewbid','Rentlevel','Bidshare','Return')


% 
% 
% plot(TSTART:TMAX-1,Bprojlevel(10,TSTART:TMAX-1)./max(Bprojlevel(10,TSTART:TMAX-1)),'c')
% plot(TSTART:TMAX-1,Diffrentset(10,TSTART:TMAX-1)./max(Diffrentset(10,TSTART:TMAX-1)),'r')
% legend('Avgexptdiff','htperyear','Bprojdiff','Rentdiff')

% % %Small Houses% figure
% TSTART=10;
% TMAX=30;
% figure
% subplot(5,1,1)
% plot(TSTART:TMAX-1,htperyear([1 4 7 10],TSTART:TMAX-1))
% title('Number of Small Houses Built')
% legend('HT 1','HT 4','HT 7','HT 10','Best')
% subplot(5,1,2)
% plot(TSTART:TMAX-1,Diffset(1,TSTART:TMAX-1))
% hold on
% plot(TSTART:TMAX-1,Diffrentset(1,TSTART:TMAX-1),'r')
% title('Change between Rent(t-1) to Projected Rent(t)(blue) and Rent(t)(red) for HT 1')    
% subplot(5,1,3)
% plot(TSTART:TMAX-1,Diffset(4,TSTART:TMAX-1))
% hold on
% plot(TSTART:TMAX-1,Diffrentset(4,TSTART:TMAX-1),'r')
% title('Change between Rent(t-1) to Projected Rent(t)(blue) and Rent(t)(red) for HT 4')    
% subplot(5,1,4)
% plot(TSTART:TMAX-1,Diffset(7,TSTART:TMAX-1))
% hold on
% plot(TSTART:TMAX-1,Diffrentset(7,TSTART:TMAX-1),'r')
% title('Change between Rent(t-1) to Projected Rent(t)(blue) and Rent(t)(red) for HT 7')    
% subplot(5,1,5)
% plot(TSTART:TMAX-1,Diffset(10,TSTART:TMAX-1))
% hold on
% plot(TSTART:TMAX-1,Diffrentset(10,TSTART:TMAX-1),'r')
% title('Change between Rent(t-1) to Projected Rent(t)(blue) and Rent(t)(red) for HT 10')    

% 
% % %Medium Houses
% figure
% TSTART=10;
% TMAX=30;
% subplot(5,1,1)
% plot(TSTART:TMAX-1,htperyear([2 5 8 11],TSTART:TMAX-1))
% title('Number of Medium Houses Built')
% legend('HT 2','HT 5','HT 8','HT 11','Best')
% subplot(5,1,2)
% plot(TSTART:TMAX-1,Diffset(2,TSTART:TMAX-1))
% hold on
% plot(TSTART:TMAX-1,Diffrentset(2,TSTART:TMAX-1),'r')
% title('Change between Rent(t-1) to Projected Rent(t)(blue) and Rent(t)(red) for HT 2')    
% subplot(5,1,3)
% plot(TSTART:TMAX-1,Diffset(5,TSTART:TMAX-1))
% hold on
% plot(TSTART:TMAX-1,Diffrentset(5,TSTART:TMAX-1),'r')
% title('Change between Rent(t-1) to Projected Rent(t)(blue) and Rent(t)(red) for HT 5')    
% subplot(5,1,4)
% plot(TSTART:TMAX-1,Diffset(8,TSTART:TMAX-1))
% hold on
% plot(TSTART:TMAX-1,Diffrentset(8,TSTART:TMAX-1),'r')
% title('Change between Rent(t-1) to Projected Rent(t)(blue) and Rent(t)(red) for HT 8')    
% subplot(5,1,5)
% plot(TSTART:TMAX-1,Diffset(11,TSTART:TMAX-1))
% hold on
% plot(TSTART:TMAX-1,Diffrentset(11,TSTART:TMAX-1),'r')
% title('Change between Rent(t-1) to Projected Rent(t)(blue) and Rent(t)(red) for HT 11')    
% 
% % %Large Houses
% figure
% TSTART=10;
% TMAX=30;
% subplot(5,1,1)
% plot(TSTART:TMAX-1,htperyear([3 6 9 12],TSTART:TMAX-1))
% title('Number of Large Houses Built')
% legend('HT 3','HT 6','HT 9','HT 12','Best')
% subplot(5,1,2)
% plot(TSTART:TMAX-1,Diffset(3,TSTART:TMAX-1))
% hold on
% plot(TSTART:TMAX-1,Diffrentset(3,TSTART:TMAX-1),'r')
% title('Change between Rent(t-1) to Projected Rent(t)(blue) and Rent(t)(red) for HT 2')    
% subplot(5,1,3)
% plot(TSTART:TMAX-1,Diffset(6,TSTART:TMAX-1))
% hold on
% plot(TSTART:TMAX-1,Diffrentset(6,TSTART:TMAX-1),'r')
% title('Change between Rent(t-1) to Projected Rent(t)(blue) and Rent(t)(red) for HT 6')    
% subplot(5,1,4)
% plot(TSTART:TMAX-1,Diffset(9,TSTART:TMAX-1))
% hold on
% plot(TSTART:TMAX-1,Diffrentset(9,TSTART:TMAX-1),'r')
% title('Change between Rent(t-1) to Projected Rent(t)(blue) and Rent(t)(red) for HT 9')    
% subplot(5,1,5)
% plot(TSTART:TMAX-1,Diffset(12,TSTART:TMAX-1))
% hold on
% plot(TSTART:TMAX-1,Diffrentset(12,TSTART:TMAX-1),'r')
% title('Change between Rent(t-1) to Projected Rent(t)(blue) and Rent(t)(red) for HT 12')    
% %@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

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

avghouserent=zeros(HT,1);
avgrentgrad=zeros(HT,1);
avgrentdynms=zeros(HT,TMAX-TSTART);
rentdynstats=zeros(HT,TMAX-TSTART);
totdev=devcells(1,TMAX);
pctdev=devcells(2,TMAX);
ratediff=zeros(1,TMAX);
maxdist=max(DISTANCE(iurb));
totvacrate=vacrate;
popchange=diff(POP(TSTART+1:TMAX))./POP(TSTART+1:TMAX-1);
lotchange=numlt(:,TMAX)-numlt(:,TSTART);
deltalots=diff(Nlots(TSTART+1:TMAX))./Nlots(TSTART+1:TMAX-1);
ldmean=mean(landdemand);
ldstd=std(landdemand);
deltadev=diff(devcells(2,TSTART:TMAX));
deltaagr=diff(agrland(TSTART+1:TMAX))./agrland(TSTART+1:TMAX-1);
avgrentdynms=diff(avgrent(:,TSTART:TMAX),1,2);
pctbiglots=length(find(LOTTYPE > 9))/length(find(BASELAYER==1));

for lt=1:HT
    rentdynstats(lt,:)=avgrentdynms(lt,:)./avgrent(lt,TSTART:TMAX-1);
end

% Mean dispersion (Irwin and Bockstael, 2007)
altprop=zeros(length(iurblist),1);
for ic=1:length(iurblist)
    [row,col]=ind2sub([NLENGTH NWIDTH],iurblist(ic));
    updir=max(row-2,1);
    dndir=min(row+2,NLENGTH);
    lfdir=max(col-2,1);
    rtdir=min(col+2,NWIDTH);
    altprop(ic)=length(find(BASELAYER(updir:dndir,lfdir:rtdir)==0))/...
        (length(updir:dndir)*length(lfdir:rtdir));
end

meandisp=sum(altprop)/length(iurblist);

maxdevdist=max(DISTANCE(iurb));

distzones=ceil(max(max(DISTANCE))/13);
zonedensity=zeros(distzones,1);
diststart=0;
distlim=indevedge;

for dz=1:distzones
    idistzone=find(DISTANCE > diststart & DISTANCE <= distlim);
    zonedensity(dz)=mean(z(LOTTYPE(idistzone(ismember(idistzone,iurblist))),1));
    diststart=diststart+indevedge;
    distlim=distlim+indevedge;
end

% <><><><><> plot 3x1 avgrent, htperyear,htincome <><><><>%
% % by house size
% 
% %small houses
% figure
% TSTART=10;
% TMAX=30;
% subplot(3,1,1)
% plot(TSTART:TMAX-1,avgrent([1 4 7 10],TSTART:TMAX-1))
% title('Average Rent For Small Houses')
% subplot(3,1,2)
% plot(TSTART:TMAX-1,htperyear([1 4 7 10],TSTART:TMAX-1))
% title('Number of Small Houses Built')
% subplot(3,1,3)
% plot(TSTART:TMAX-1,htincome([1 4 7 10],TSTART:TMAX-1))
% title('Average Income For Small Houses')
% legend('HT 1','HT 4','HT 7','HT 10','SouthEast')
% % 
% % %medium houses
% figure
% TSTART=10;
% TMAX=30;
% subplot(3,1,1)
% plot(TSTART:TMAX-1,avgrent([2 5 8 11],TSTART:TMAX-1))
% title('Average Rent For Medium Houses')
% subplot(3,1,2)
% plot(TSTART:TMAX-1,htperyear([2 5 8 11],TSTART:TMAX-1))
% title('Number of Medium Houses Built')
% subplot(3,1,3)
% plot(TSTART:TMAX-1,htincome([2 5 8 11],TSTART:TMAX-1))
% title('Average Income For Medium Houses')
% legend('HT 2','HT 5','HT 8','HT 11','SouthEast')
% % 
% % %large houses
% figure
% subplot(3,1,1)
% plot(TSTART:TMAX-1,avgrent([3 6 9 12],TSTART:TMAX-1))
% title('Average Rent For Large Houses')
% subplot(3,1,2)
% plot(TSTART:TMAX-1,htperyear([3 6 9 12],TSTART:TMAX-1))
% title('Number of Large Houses Built')
% subplot(3,1,3)
% plot(TSTART:TMAX-1,htincome([3 6 9 12],TSTART:TMAX-1))
% title('Average Income For Large Houses')
% legend('HT 3','HT 6','HT 9','HT 12','SouthEast')

% % by lot size
% figure
% subplot(3,1,1)
% plot(TSTART:TMAX-1,avgrent(10:12,TSTART:TMAX-1))
% title('Average Rent For 2-acre Lots')
% subplot(3,1,2)
% plot(TSTART:TMAX-1,htperyear(10:12,TSTART:TMAX-1))
% title('Number of 2-acre Lots Built')
% subplot(3,1,3)
% plot(TSTART:TMAX-1,htincome(10:12,TSTART:TMAX-1))
% title('Average Income For 2-acre Lots')
% legend('HT 10','HT 11','HT 12','SouthEast')
% % 
% figure
% subplot(3,1,1)
% plot(TSTART:TMAX-1,avgrent(4:6,TSTART:TMAX-1))
% title('Average Rent For 0.5-acre Lots')
% subplot(3,1,2)
% plot(TSTART:TMAX-1,htperyear(4:6,TSTART:TMAX-1))
% title('Number of 0.5-acre Lots Built')
% subplot(3,1,3)
% plot(TSTART:TMAX-1,htincome(4:6,TSTART:TMAX-1))
% title('Average Income For 0.5-acre Lots')
% legend('HT 4','HT 5','HT 6','SouthEast')
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

%Plot model and proj by type over time
ltsz=[0.25 0.25 0.25 0.5 0.5 0.5 1 1 1 2 2 2 5 5 5 10 10 10]';
bmodedist=zeros(HT,LTYPE,TMAX);
projdist=zeros(HT,TMAX);
avgrentdist=zeros(HT,TMAX);

for t=1:TMAX
    for im=1:LTYPE
        for lt=1:HT
            sublotmap=Dynltmap(:,:,t);
            submodelmap=Bmodelmap(:,:,t);
            subprojmap=Bprojmap(:,:,t);
            subrentmap=Dynrentmap(:,:,t);
            ifindlt=(sublotmap==lt);
            bmodedist(lt,im,t)=length(find(submodelmap(ifindlt)==im))/ltsz(lt);
            projdist(lt,t)=mean(subprojmap(ifindlt));
            avgrentdist(lt,t)=mean(subrentmap(ifindlt));
        end
    end
end
projdist(isnan(projdist))=0;
bprojdiff=diff(projdist(1:HT,TSTART+1:TMAX),1,2);
avgrentdiff=diff(avgrent(1:HT,TSTART+1:TMAX),1,2);


bidshare=zeros(HT,30);
buildshare=zeros(HT,30);
for t=11:30
    subshare=numtotbids(:,t)./Ltperyear(:,t);
    ibidshare=(isinf(subshare) == 0 & isnan(subshare)==0);
    bidshare(ibidshare,t)=subshare(ibidshare)/sum(subshare(ibidshare),1);
    buildshare(:,t)=htperyear(:,t)/sum(htperyear(:,t),1);
end

%<><><><><><><><><><><><><><><><>
% figure
% area(bidshare(1:12,11:30)')
% axis([1 20 0 1])
% title('Share of Bids Per Housing Type')
% xlabel('Time Step')
% 
% figure
% area(buildshare(1:12,11:30)')
% axis([1 20 0 1])
% title('Share of Housing Type Built')
% xlabel('Time Step')
% <><><><><><><><><><><><><><><><>
Buildtimemap=zeros(80,80);
Buildtimemap(Dynltmap(:,:,11)~=0)=11;
% Bprojdiffmap=zeros(80,80,TMAX);
% Dynrentdiffmap=zeros(80,80,TMAX);
% Diffmap=zeros(80,80,TMAX);

Bprojdiffmap=zeros(80,80);
Dynrentdiffmap=zeros(80,80);
Diffmap=zeros(80,80);


for ti=TSTART+2:TMAX
    subprojdiffmap=zeros(80,80);
    subrentdiffmap=zeros(80,80);
    subdiffmap=zeros(80,80);
    Buildtimemap((Dynltmap(:,:,ti)-Dynltmap(:,:,ti-1))~=0)=ti;
    ibtpre=(Buildtimemap==ti-1);
    ibtpost=(Buildtimemap==ti);
    for lt=1:HT
        submap5=Dynltmap(:,:,ti)-Dynltmap(:,:,ti-1);
        submap4=Dynltmap(:,:,ti-1);
        submappre=Bprojmap(:,:,ti-1);
        submappost=Bprojmap(:,:,ti);
        subrentpre=Dynrentmap(:,:,ti-1);
        subrentpost=Dynrentmap(:,:,ti);
        ifltpre=find(submap4 == lt & Buildtimemap == ti-1);
        ifltpost=find(submap5 == lt & Buildtimemap == ti);
        
        if isempty(find(ifltpost,1)) == 0
        % Absolute change
        Bprojdiffmap(ifltpre)=mean(submappost(ifltpost))-...
                mean(subrentpost(ifltpost));
        Dynrentdiffmap(ifltpost)=mean(subrentpost(ifltpost))-mean(...
            subrentpre(ifltpre));
        Diffmap(ifltpre)=mean(submappost(ifltpost))-...
                mean(subrentpost(ifltpost));
        Diffmap(ifltpost)=mean(subrentpost(ifltpost))-mean(...
                subrentpre(ifltpre));
            
        Diffset(lt,1,ti)=mean(submappost(ifltpost))-...
                mean(subrentpost(ifltpost));
        Diffset(lt,2,ti)=mean(subrentpost(ifltpost))-mean(...
            subrentpre(ifltpre));
         
         end

        
        % percent change
%         Bprojdiffmap(ifltpost)=(mean(submappost(ifltpost))-mean(...
%             subrentpost(ifltpost)))/mean(submappre(ifltpost));
%         Dynrentdiffmap(ifltpost)=(mean(subrentpost(ifltpost))-mean(...
%             subrentpre(ifltpre)))/mean(subrentpre(ifltpre));
%         Bprojdiffmap(ifltpost)=(mean(submappost(ifltpost))-mean(...
%             submappre(ifltpre)))/mean(submappre(ifltpre));
        
%         subprojdiffmap(ifltpost)=(mean(submappost(ifltpost))-mean(...
%             submappre(ifltpre)))/mean(submappre(ifltpre));
%         subrentdiffmap(ifltpost)=(mean(subrentpost(ifltpost))-mean(...
%             subrentpre(ifltpre)))/mean(subrentpre(ifltpre));
    end
%     Diffmap(:,:,ti)=subdiffmap;
%     Bprojdiffmap(:,:,ti)=subprojdiffmap;
%     Dynrentdiffmap(:,:,ti)=subrentdiffmap;
end

%%% Ideal housing offerings
% idealset=zeros(HT,1);
% 
% for ires=1:length(con2lot(ifilled,1))
%     hopt=((Income(con2lot(ifilled(ires),2),1)-travelcost(lotchoice(ifilled(ires),2))-...
%         avgrent(:,TMAX)).^CONINFO(con2lot(ifilled(ires),2),1)).*(z(:,2).^...
%         CONINFO(con2lot(ifilled(ires),2),2)).*(z(:,1).^CONINFO(con2lot(ifilled(ires),2),3)).*...
%         (z(lt,3).^CONINFO(con2lot(ifilled(ires),2),4));
%     [imaxu,jmaxu]=max(hopt,[],1);
%     idealset(jmaxu)=idealset(jmaxu)+1; 
% end
subrmap=SUBRISKMAP(:,:,5);
colmean=zeros(1,NWIDTH);
colmedian=zeros(1,NWIDTH);
colmax=zeros(1,NWIDTH);
colmin=zeros(1,NWIDTH);
for col=1:NWIDTH
    inonzero=(subrmap(:,col)>0);
    if isempty(find(inonzero,1))==1
        continue
    end
    colmean(1,col)=mean(subrmap(inonzero,col));
    colmedian(1,col)=median(subrmap(inonzero,col));
    colmin(1,col)=min(subrmap(inonzero,col));
    colmax(1,col)=max(subrmap(inonzero,col));
end
damcoef_map=zeros(NLENGTH,NWIDTH);
for k=1:length(ifilled)
    isamecell=find(Lottype(:,1)==con2lot(ifilled(k),3),1,'first');
    subdamcoef=damcoef(:,:,con2lot(ifilled(k),2),BUILDTIME(Lottype(isamecell,2)));
    damcoef_map(Lottype(isamecell,2))=subdamcoef(Lottype(isamecell,2));
end