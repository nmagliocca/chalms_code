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
COAST=zeros(NLENGTH,NWIDTH);
COAST(:,1:5)=1;
SCAPE(COAST~=1)=1;
iscape=(SCAPE==1);
iscapelist=find(iscape==1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Distance matrices   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

icenterrow=1;
% icentercol=round(NWIDTH/2);
icentercol=6;

for col=1:NWIDTH
    dist2hznnei(1:NLENGTH,col)=abs(col-icentercol).*ones(NLENGTH,1);
end

for row=1:NWIDTH
    dist2vrtnei(row,1:NWIDTH)=abs(row-icenterrow).*ones(1,NLENGTH);
end

for col=1:NWIDTH
    for row=1:NLENGTH
        DISTANCE(row,col)=sqrt(dist2hznnei(row,col)^2+dist2vrtnei(row,col)^2);   
    end
end

travelcost(:,:)=margtc.*DISTANCE;
coastdist=cumsum(SCAPE,2);
coastprox=10*(max(max(coastdist))+1-coastdist);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%     Impact Surface    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxPflood=0.7;
highrisk=30;
risksurf=1./(1+exp((coastdist-highrisk)/10));
Pflood=maxPflood*(risksurf+(1-max(max(risksurf))));
TSI=ones(size(SCAPE));
% siterisk=0.1*randn(size(SCAPE));
siterisk=zeros(size(SCAPE));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ZONES Layer    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nzoneslong=sqrt(NZONES);
nzoneswide=sqrt(NZONES);
extralengthz=rem(NLENGTH,nzoneslong);
extrawidthz=rem(NWIDTH,nzoneswide);
izonelength=(length(BASELAYER(:,1))-extralengthz)/nzoneslong;
izonewidth=(length(BASELAYER(1,:))-extrawidthz)/nzoneswide;

zonemarklong=1;
zonemarkwide=1;
for ii=1:nzoneswide
    for jj=1:nzoneslong
        ZONES(zonemarklong:zonemarklong+izonelength-1,...
            zonemarkwide:zonemarkwide+izonewidth-1)=...
            ii*nzoneslong-(nzoneslong-jj);        
        
        if jj==nzoneslong && extralengthz > 0
            ZONES(izonelength*jj+1:izonelength*jj+extralengthz,...
                zonemarkwide:zonemarkwide+izonewidth-1)=...
                ii*nzoneslong-(nzoneslong-jj);        
        end
        if ii==nzoneswide && extrawidthz > 0
            ZONES(zonemarklong:zonemarklong+izonelength-1,...
                izonewidth*ii+1:izonewidth*ii+extrawidthz)=...
                ii*nzoneslong-(nzoneslong-jj);
        end
        zonemarklong=mod(izonelength*jj+1,izonelength*nzoneslong);
    end     
    if jj==nzoneslong && ii==nzoneswide && extralengthz > 0
        ZONES(izonelength*jj+1:izonelength*jj+extralengthz,...
            izonewidth*ii+1:izonewidth*ii+extrawidthz)=...
            ii*nzoneslong-(nzoneslong-jj);
    end
    zonemarkwide=mod(izonewidth*ii+1,izonewidth*nzoneswide);
end

for n=1:NZONES
    str= sprintf('izone%d = find(ZONES == %d);',n,n);
    eval(str);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ZONING   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

house2acre=[2 0.5 0.2]';
zs=min(max(1./house2acre,0),5);

for zt=1:NZONES
%     zoning(zt,:)=[1 1];
    if isempty(find(zt==[1:2 6:8 11:13 16:18 21:22],1))==0
        zoning(zt,:)=[min(zs) max(zs)]; %no zoning
        ZONEMAP(ZONES==zt)=zoning(zt,1);
        ZONED(ZONES==zt)=0;
    elseif isempty(find(zt==[3:5 9:10 14:15 19:20 23:25],1))==0
%         zoning(zt,:)=[5 max(zs)];
        zoning(zt,:)=[min(zs) max(zs)];
        ZONEMAP(ZONES==zt)=zoning(zt,1);
        ZONED(ZONES==zt)=1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%    Broker Layer    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbrokerlong=(NLENGTH/5);
nbrokerwide=(NWIDTH/5);
extralengthb=rem(NLENGTH,nbrokerlong);
extrawidthb=rem(NWIDTH,nbrokerwide);
ibrokerlength=(length(BASELAYER(:,1))-extralengthb)/nbrokerlong;
ibrokerwidth=(length(BASELAYER(1,:))-extrawidthb)/nbrokerwide;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%% Landscape Template %%%%%%%%%%%%%%%%%%%%%%%%%%%%

stream.Substream=86;
stream.State=repeatstate1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%    Housing Layer    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nlots=zeros(1,TMAX);
inurbarea=288;
globeind=[(1:NLENGTH*NWIDTH)' reshape(DISTANCE,NLENGTH*NWIDTH,1)];
globeind=globeind(ismember(globeind(:,1),iscapelist),:);
iurbprox=sortrows(globeind,2);
BASELAYER(iurbprox(1:inurbarea,1))=1;
% inurbpct=0.166*max(DISTANCE);    %initial urban extent
% for ii=1:NWIDTH
%     BASELAYER(:,ii)=(DISTANCE(:,ii) <= inurbpct(ii) & SCAPE(:,ii) == 1);
% end

iurb=(BASELAYER==1 & SCAPE==1);
iurblist=find(iurb==1);
iagr=(BASELAYER==0 & SCAPE==1);
iagrlist=find(iagr==1);


indevedge=find(BASELAYER(:,icentercol)==1,1,'last');


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
LOTS(iagr)=Nlots(1)+1000;
LOTS(SCAPE==0)=Nlots(1)+1000;

%%% Assign 2 acre lots
adjcells=diff(find(LOTS==0));
if isempty(find(adjcells==1,1))==0
    for f=1:Nlots(1)
        mylots=find(LOTS==f);
        for ml=1:length(mylots)
            zoneid=ZONES(mylots(ml));
            growthlim=2-length(mylots);
            if growthlim == 0
                continue
            else
                [frows,fcols]=ind2sub([NLENGTH,NWIDTH],mylots(ml));
                rnei=(fcols+1)*NLENGTH-(NLENGTH-frows);
                lnei=(fcols-1)*NLENGTH-(NLENGTH-frows);
                upnei=fcols*NLENGTH-(NLENGTH-(frows-1));
                dnnei=fcols*NLENGTH-(NLENGTH-(frows+1));
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
ilotfill=find(LOTS==0);
LOTS(ilotfill)=Nlots(1)+(1:length(ilotfill));
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

housesize=[1500 2000 2500]';
house2cells=unique(LOTS((LOTS < 250)));

%%%% Amenity Level
% amlevel=[100 300 200 100 300 200 400]';
% amlevel=[100 200 150 100 200 150 400]';

for ih=1:length(house2cells)
    isamecell=find(LOTS == house2cells(ih));
    HOUSESIZE(isamecell)=housesize(ceil(length(housesize)*rand(1)));
%     AMLEVEL(isamecell)=amlevel(ceil(length(housesize)*rand(1)));
    %%% Starting with all conventional subdivision types %%%
%     HOUSESIZE(isamecell)=housesize(1);
%     AMLEVEL(isamecell)=amlevel(1);
%     AMLEVEL(isamecell)=max(max(coastdist))-mean(coastdist(isamecell));
    AMLEVEL(isamecell)=mean(coastprox(isamecell));
end

z=[reshape(repmat(zs,1,3)',HT,1) repmat(housesize,3,1)];

% devpref=[6; 4.67; 5.75; 4.67; 4.58; 4.33; 4.11]./7;
% devpref=[0.851; 0.646; 0.823; 0.620; 0.604; 0.601; 0.641];
devpref=ones(HT,1);

%%% Define and locate similar lots
for lt=1:HT
    if z(lt,1) < 1
        simlotrange(lt,1)=1;
        simlotrange(lt,2)=find(z(:,1) < 1,1,'last');
    elseif z(lt,1) >=1 && z(lt,1) <=2
        simlotrange(lt,1)=find(z(:,1) >= 1,1,'first');
        simlotrange(lt,2)=find(z(:,1) <= 2,1,'last');
    elseif z(lt,1) > 1
        simlotrange(lt,1)=find(z(:,1) > 1,1,'first');
        simlotrange(lt,2)=HT;
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

Lottype=zeros(length(Lotinfo(:,1)),8);     %[id index lotsize housesize ltype ccost occ/vac amlevel]
lotchoice=zeros([],8);
cellinfo=zeros(length(iurblist),8);

Lottype(:,[1:2 4 8])=Lotinfo;
for lt=1:length(Lotinfo)
    isamecell=(Lottype(:,1)==Lotinfo(lt,1));
    Lottype(lt,3)=length(find(Lotinfo(:,1)==Lotinfo(lt,1)))/...
        length(find(Lotinfo(:,2)==Lotinfo(lt,2)));
    Lottype(lt,5)=find(z(:,1)==Lottype(lt,3) & z(:,2)==Lottype(lt,4));
end

[lotids,ilots,jnum]=unique(Lottype(:,1));
lotchoice(lotids,:)=Lottype(ilots,:);

[cellids,icells,jnumcell]=unique(Lottype(:,2));
cellinfo(1:length(iurblist),:)=Lottype(icells,:);
LOTTYPE(cellinfo(:,2))=cellinfo(:,5);


Nlots(1:TSTART+1)=length(lotchoice(:,1));
indlots=lotchoice(:,3);

m=zeros(Nlots(1),1);
m(1:Nlots(1))=lotchoice(1:Nlots(1),4)./(lotchoice(1:Nlots(1),3)*acre2sqft);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  Agricultural Layer   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stream.Substream=35;
stream.State=repeatstate2;

load FARMMAP

% %<><><><><><><><><><> Randomized Land Allocation <><><><><><><><><><><><><>
% startpos=zeros(Nfarmers,1);
% highpop=round(Nfarmers*0.65);
% posid=1:NLENGTH*NWIDTH;
% posid=reshape(posid,NLENGTH,NWIDTH);
% subtop=posid(1:NLENGTH/2,:);
% subbot=posid(NLENGTH/2+1:NLENGTH,:);
% 
% 
% startpos(1:highpop,1)=subtop(ceil((NLENGTH*NWIDTH)/2*rand(highpop,1)));
% startpos(highpop+1:Nfarmers,1)=subbot(ceil((NLENGTH*NWIDTH)/2*rand(Nfarmers-highpop,1)));
% 
% for np=1:length(startpos)
%     check=find(startpos==startpos(np));
%     while length(check)>1
%         ipick=ceil(length(check)*rand(1));
%         if ismember(startpos(check(ipick)),subtop)==1
%             startpos(check(ipick))=subtop(ceil((NLENGTH*NWIDTH)/2*rand(1)));
%         else
%             startpos(check(ipick))=subbot(ceil((NLENGTH*NWIDTH)/2*rand(1)));
%         end
%         check=find(startpos==startpos(np));
%     end
% end
% 
% while length(unique(startpos)) ~= length(startpos)
%     startpos(BASELAYER(startpos)==1)=min(ceil(NLENGTH*NWIDTH*rand(1)),NLENGTH*NWIDTH);
%     if (length(unique(startpos)) ~= length(startpos))==1
%         for k=1:length(startpos)
%             check=find(startpos(k)==startpos);
%             if length(check) >1
%                 startpos(check)=min(ceil(NLENGTH*NWIDTH*rand(length(check),1)),NLENGTH*NWIDTH);
%             end
%         end
%     end
% end
% AGLAYER(startpos)=(1:Nfarmers);
% 
% while isempty(find(AGLAYER==0,1))==0
%     for f=1:Nfarmers
%         myland=find(AGLAYER==f);
%         for ml=1:length(myland)
%             [frows,fcols]=ind2sub([NLENGTH,NWIDTH],myland(ml));
%             rnei=(fcols+1)*NLENGTH-(NLENGTH-frows);
%             lnei=(fcols-1)*NLENGTH-(NLENGTH-frows);
%             upnei=fcols*NLENGTH-(NLENGTH-(frows-1));
%             dnnei=fcols*NLENGTH-(NLENGTH-(frows+1));
%             distcells=find(abs(DISTANCE-DISTANCE(myland(ml)))<=1 & AGLAYER==0);
%             edgecells=find(distcells==rnei | distcells==lnei | distcells==upnei | ...
%                 distcells==dnnei);
%             AGLAYER(distcells(edgecells))=f;
%         end
%     end
% end

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@    AGENTS    @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% <><><><><><><><><><><><><><>    FARMERS    <><><><><><><><><><><><><><><>

idNfarmers=unique(AGLAYER(iscape));

Farmstats(iNfarmers,1,1)=idNfarmers(iNfarmers);
Farmstats(iNfarmers,3,1)=FARMPROD+PRODSTD*randn(Nfarmers,1);
Farmstats(iNfarmers,4,1)=FARMCOST+COSTSTD*randn(Nfarmers,1);
Farmstats(iNfarmers,5,1)=normrnd(AVGFARMRETURN,STDFARMRETURN,Nfarmers,1);

%@@@@@@@ Farmer Projection Models @@@@@@@@@@@@@@@@@@

%Distance-Discounting models for farmers
distcoeff=mincoeff+(maxcoeff-mincoeff)*rand(Nfarmers,NUMMODELDIST);
%distance coefficient is in $1000/acre_distance

landmodel = ceil(FARMNUMCLASS*rand(Nfarmers,NUMMODEL));
for i = 1:FARMNUMCLASS
    strl = sprintf('landclass%d = find(landmodel == %d);',i,i);
    eval(strl);
end

aa = zeros(Nfarmers,NUMMODEL);    %land models
for i = 1:FARMNUMCLASS
    if i == 1
        % mirror model
        aa(landclass1) = rand(1); % fraction that pred is away from 1/2 from mirror image
    elseif i == 2
        % mean model
        aa(landclass2) = ceil(MAXMEANMODEL*rand(length(landclass2),1));
    elseif i == 3
        %cycle model
        aa(landclass3) = ceil(MAXCYCLEMODEL*rand(length(landclass3),1));
    elseif i == 4
        % projection model
        aa(landclass4) = ceil(2+((MAXPROJECT-1)-2)*rand(length(landclass4),1));
    elseif i == 5
        % rescale model
        aa(landclass5) = 2*rand(length(landclass5),1);
    elseif i == 6
        %local or regional trends
        aa(landclass6) = round(rand(length(landclass6),1));
    end
end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%               Learning Period for Financial Prediction Models
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


FarmerModule_Coast_1v1

%%%%%%%%%%%%%%%   ADD LEARNED INFO   %%%%%%%%%%%%%%%

% %Implicitly added: 
% landerror
% distcoeff

fitness(:,:,TSTART)=fitnessSAVE;
landerror=learnlanderror;
wtaland(iNfarmers,1:TSTART+1)=learnwtaland(iNfarmers,tlearn-TSTART:tlearn);
Paskland(iNfarmers,1:TSTART+1)=Paskland(iNfarmers,tlearn-TSTART:tlearn);

landbestSAVE(iNfarmers,1:TSTART)=learnlandbestSAVE(iNfarmers,tlearn-TSTART:tlearn-1);
ilandbestSAVE(iNfarmers,1:TSTART)=ilearnlandbestSAVE(iNfarmers,tlearn-TSTART:tlearn-1);
landprojbestSAVE(iNfarmers,1:TSTART)=learnlandprojbestSAVE(iNfarmers,tlearn-TSTART:...
    tlearn-1);
landmodelbestSAVE(iNfarmers,1:TSTART)=learnlandmodelbestSAVE(iNfarmers,...
    tlearn-TSTART:tlearn-1);
Plandproj(iNfarmers,1:TSTART)=learnlandprojbestSAVE(iNfarmers,tlearn-TSTART:tlearn-1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%    Developers    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stream.Substream=4;
stream.State=repeatstate3;
%%%% Developer's Population Prediction Models %%%%%%%

classagentmodel = ceil(POPNUMCLASS*rand(Ndevelopers,NUMMODEL));

for i = 1:POPNUMCLASS
    str = sprintf('indclass%d = find(classagentmodel == %d);',i,i);
    eval(str);
end

dd = zeros(Ndevelopers,NUMMODEL);
for i = 1:POPNUMCLASS
    if i == 1
        % mirror model
        dd(indclass1) = 0.40+(0.60-0.40)*rand(1,length(indclass1)); % fraction that pred is away from 1/2 from mirror image
    elseif i == 2
        % mean model
        dd(indclass2) = ceil(MAXMEANMODEL*rand(length(indclass2),1));
    elseif i == 3
        %cycle model
        dd(indclass3) = ceil(MAXCYCLEMODEL*rand(length(indclass3),1));
    elseif i == 4
        % projection model
        dd(indclass4) = 1+ceil((MAXPROJECT-1)*rand(length(indclass4),1));
    elseif i == 5
        % rescale model
        dd(indclass5) = 2*rand(length(indclass5),1);
    end
end

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
    for i = 1:BROKERNUMCLASS
        strb = sprintf('brokerclass%d = find(brokermodel(:,:,nb) == %d);',i,i);
        eval(strb);
    end
    bb=zeros(HT,NUMMODEL);
    for i = 1:BROKERNUMCLASS
%         bb=bbfull(:,:,nb);
        if i == 1
            % mirror model
            bb(brokerclass1) = rand(1); % fraction that pred is away from 1/2 from mirror image
        elseif i == 2
            % mean model
            bb(brokerclass2) = ceil(MAXMEANMODEL*rand(length(brokerclass2),1));
        elseif i == 3
            %cycle model
            bb(brokerclass3) = ceil(MAXCYCLEMODEL*rand(length(brokerclass3),1));
        elseif i == 4
            % projection model
            bb(brokerclass4) = ceil(2+((MAXPROJECT-1)-2)*rand(length(brokerclass4),1));
        elseif i == 5
            % rescale model
            bb(brokerclass5) = 2*rand(length(brokerclass5),1);
        elseif i == 6
            %local or regional trends
            bb(brokerclass6) = ceil(rand(length(brokerclass6),1));
        end
    end
    bbfull(:,:,nb)=bb;
end
clear brokerclass1 brokerclass2 brokerclass3 brokerclass4 brokerclass5 brokerclass6 bb

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%               Learning Period for Broker Prediction Models
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


BrokersModule_Coast_0v3

%%%%%%%%%%%%%%%   ADD LEARNED INFO   %%%%%%%%%%%%%%%
brokererror(:,:,:,1:TSTART)=learnbrokererror(:,:,:,tlearn-9:tlearn);
brokerbestdiffSAVE(:,:,1:TSTART)=learnbestdiffSAVE(:,:,max(12,tlearn-9):tlearn);
brokerabserror(:,:,:,1:TSTART)=learnabserror(:,:,:,tlearn-9:tlearn);
brokerbestabsSAVE(:,:,1:TSTART)=learnbestabsSAVE(:,:,max(12,tlearn-9):tlearn);
% <><><><><><><><><><><><><>    Consumers    <><><><><><><><><><><><><><><>

% Nconsumers=round(Nlots(TSTART)/0.67);
Nconsumers=Nlots(TSTART);
housepref=zeros(Nconsumers,1);
HU=zeros(Nconsumers,1);
subincome=round(exp(lognrnd(parmhat(1),parmhat(2),Nconsumers,1)));
Income(1:Nconsumers,1)=min(max(subincome,minwage),maxwage);
Income(1:Nconsumers,2)=TSTART+ceil(searchtimemin+(searchtimemax-searchtimemin)*rand(Nconsumers,1));

income3=(Income(:,1) >= lowminwage & Income(:,1) <= lowmaxwage);     %income 1 2 3 = hi mid low wages
income2=(Income(:,1) >= midminwage & Income(:,1) <= midmaxwage);
income1=(Income(:,1) >= himinwage & Income(:,1) <= himaxwage);
housepref(income1)=HIBETA(1)+(HIBETA(2)-HIBETA(1))*rand(length(find(income1==1)),1);
housepref(income2)=MIDBETA(1)+(MIDBETA(2)-MIDBETA(1))*rand(length(find(income2==1)),1);
housepref(income3)=LOWBETA(1)+(LOWBETA(2)-LOWBETA(1))*rand(length(find(income3==1)),1);

% %%%%%%%%%%%%%%%% Representative agent set %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Income(income1)=mean([himaxwage himinwage]);
% Income(income2)=mean([midmaxwage midminwage]);
% Income(income3)=mean([lowmaxwage lowminwage]);
% housepref(income1)=0.15+(0.15-0.15)*rand(length(find(income1==1)),1);
% housepref(income2)=0.25+(0.25-0.25)*rand(length(find(income2==1)),1);
% housepref(income3)=0.35+(0.35-0.35)*rand(length(find(income3==1)),1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subjective risk perception
damcoef(:,:,1:Nconsumers,1:TSTART)=1;
% subrisk=ones(Nconsumers,1);
subrisk=rand(Nconsumers,1);

Hbeta=housepref;
beta=(0.1+(0.9-0.1)*rand(length(Hbeta),1)).*Hbeta;     %cite Carliner
gamma=rand(length(Hbeta),1).*(Hbeta-beta)/2;   %cite Carliner
ampref=Hbeta-(beta+gamma);
alpha=(1-housepref);