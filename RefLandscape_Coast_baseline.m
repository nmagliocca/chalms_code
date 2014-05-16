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
COAST=zeros(NLENGTH,NWIDTH);
SCAPE=zeros(NLENGTH,NWIDTH);
COAST(:,1:5)=1;
icoast=find(COAST==1);
SCAPE(COAST~=1)=1;
iscape=(SCAPE==1);
iscapelist=find(iscape==1);

coastdist=cumsum(SCAPE,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Distance matrices   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load DIST2CBD
% icenterrow=1;
% % icentercol=round(NWIDTH/2);
% icentercol=6;
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
am_int=500000;
travelcost(iscapelist)=num2cell(margtc*dist2cbd(iscapelist));
travelcost(icoast)=num2cell(10000*ones(length(icoast),1));
% coastprox=num2cell(reshape(10*(max(max(coastdist))+1-coastdist),NCELLS,1));
coastprox=num2cell(reshape(am_int*exp(-0.1*coastdist),NCELLS,1));
% coastprox=num2cell(reshape(0.1*(max(max(coastdist))+1-coastdist),NCELLS,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%     Impact Surface    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxPflood=0.7;
highrisk=30;
risksurf=1./(1+exp((coastdist-highrisk)/10));
Pflood=num2cell(reshape(maxPflood*(risksurf+(1-max(max(risksurf)))),NCELLS,1));
% Pflood=mat2cell([(1:NCELLS)' reshape(maxPflood*(risksurf+(1-...
%     max(max(risksurf)))),NCELLS,1)],ones(NCELLS,1),2);
TSI=num2cell(ones(NCELLS,1));
% siterisk=0.1*randn(size(SCAPE));
% siterisk=zeros(size(SCAPE));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ZONES Layer    %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nzoneslong=sqrt(NZONES);
nzoneswide=sqrt(NZONES);
extralengthz=rem(NLENGTH,nzoneslong);
extrawidthz=rem(NWIDTH,nzoneswide);
izonelength=(NLENGTH-extralengthz)/nzoneslong;
izonewidth=(NWIDTH-extrawidthz)/nzoneswide;
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

house2acre=[2 0.5 0.2]';
zs=min(max(1./house2acre,0),5);

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

stream.Substream=86;
stream.State=repeatstate1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%    Housing Layer    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nlots=zeros(1,TMAX);
inurbarea=288;
globeind=[(1:NLENGTH*NWIDTH)' reshape(dist2cbd,NLENGTH*NWIDTH,1)];
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
ilotfill=find(LOTS==0 & SCAPE==1);
LOTS(ilotfill)=Nlots(1)+(1:length(ilotfill));
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

housesize=[1500 2000 2500]';
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
% Lottype=zeros(length(Lotinfo(:,1)),8);     %[id index lotsize housesize ltype ccost occ/vac amlevel]
% lotchoice=zeros([],8);
% cellinfo=zeros(length(iurblist),8);
% 
% Lottype(:,[1:2 4 8])=Lotinfo;
% for lt=1:length(Lotinfo)
%     isamecell=(Lottype(:,1)==Lotinfo(lt,1));
%     Lottype(lt,3)=length(find(Lotinfo(:,1)==Lotinfo(lt,1)))/...
%         length(find(Lotinfo(:,2)==Lotinfo(lt,2)));
%     Lottype(lt,5)=find(z(:,1)==Lottype(lt,3) & z(:,2)==Lottype(lt,4));
% end
% 
% [lotids,ilots,jnum]=unique(Lottype(:,1));
% lotchoice(lotids,:)=Lottype(ilots,:);
% 
% [cellids,icells,jnumcell]=unique(Lottype(:,2));
% cellinfo(1:length(iurblist),:)=Lottype(icells,:);
% LOTTYPE(cellinfo(:,2))=cellinfo(:,5);


% Nlots(1:TSTART+1)=length(lotchoice(:,1));
% indlots=lotchoice(:,3);
% 
% m=zeros(Nlots(1),1);
% m(1:Nlots(1))=lotchoice(1:Nlots(1),4)./(lotchoice(1:Nlots(1),3)*acre2sqft);


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
LANDINFO(1,1:TMAX)=mat2cell(repmat(reshape(AGLAYER,NCELLS,1),...
    1,TMAX),NCELLS,ones(1,TMAX));
sublandvalue=zeros(NCELLS,1);
subpland=zeros(NCELLS,1);
for nf=1:Nfarmers
    farmacres=find(LANDINFO{1,1}==nf);
%     farmprod=ones(length(farmacres),1)*FARMPROD+PRODSTD*randn(1);
%     farmcost=ones(length(farmacres),1)*FARMCOST+COSTSTD*randn(1);
%     farmret=ones(length(farmacres),1)*normrnd(AVGFARMRETURN,STDFARMRETURN,1,1);
    farmprod=ones(length(farmacres),1)*FARMPROD;
    farmcost=ones(length(farmacres),1)*FARMCOST;
    farmret=ones(length(farmacres),1)*AVGFARMRETURN;
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


FarmerModule_Coast_base

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

load master_dist