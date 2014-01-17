%%%%%%%%%%%%%% House Market %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wtpcon=zeros(length(Income),Nlots(t));
wtpconstar=zeros(length(Income),Nlots(t));
wtbcon=zeros(length(Income),Nlots(t));
Rn=zeros(length(Income),Nlots(t));
Phousebid=zeros(length(Income),Nlots(t));
subnhouselook=zeros(length(Income),Nlots(t));
notherbuyers=zeros(length(Income),Nlots(t));
iotherbuyers=zeros(length(Income),Nlots(t));
housemp=zeros(length(Income),1);
EUmit=zeros(length(Income),Nlots(t));
EUnomit=zeros(length(Income),Nlots(t));
U=zeros(length(Income),Nlots(t));
Unorm=zeros(length(Income),Nlots(t));
isamecell=ismember(Lottype(:,1),inewlots);
maxUset=zeros(length(Income),LTYPE*HTYPE);
exptCmit=zeros(length(Income),Nlots(t));
exptCnomit=zeros(length(Income),Nlots(t));
mitchoice=zeros(length(Income),Nlots(t));
exptcost=zeros(2,Nlots(t),length(Income));
for c=1:length(inewcon)
    subdamcoef=damcoef(:,:,inewcon(c),t);
    exptCmit(inewcon(c),inewlots)=subdamcoef(lotchoice(inewlots,2)).*Pflood(lotchoice(inewlots,2))*...
        (Cmit+miteff*Cdam)+(1-subdamcoef(lotchoice(inewlots,2))).*(1-...
        Pflood(lotchoice(inewlots,2)))*Cmit;
    exptCnomit(inewcon(c),inewlots)=subdamcoef(lotchoice(inewlots,2)).*...
        Pflood(lotchoice(inewlots,2))*Cdam;
    exptcost(:,:,inewcon(c))=[exptCmit(inewcon(c),:); exptCnomit(inewcon(c),:)];
    wtpcon(inewcon(c),inewlots)=(Income(inewcon(c))-travelcost(lotchoice(inewlots,2)))*...
        (CONINFO(inewcon(c),2)+CONINFO(inewcon(c),3)+CONINFO(inewcon(c),4))';
%     ihouseout=(Paskhouse(inewlots)' > wtpcon(inewcon(c),inewlots));
%     U(inewcon(c),inewlots)=((Income(inewcon(c))-Paskhouse(inewlots)-...
%         travelcost(lotchoice(inewlots,2))).^CONINFO(inewcon(c),1)).*...
%         (lotchoice(inewlots,4).^CONINFO(inewcon(c),2)).*(lotchoice(inewlots,3).^...
%         CONINFO(inewcon(c),3)).*(lotchoice(inewlots,8).^CONINFO(inewcon(c),4));
    EUmit(inewcon(c),inewlots)=(max(Income(inewcon(c))-travelcost(lotchoice(inewlots,2))-...
        Paskhouse(inewlots)-exptCmit(inewcon(c),inewlots)',0).^CONINFO(inewcon(c),1)).*...
        (lotchoice(inewlots,4).^CONINFO(inewcon(c),2)).*(lotchoice(inewlots,3).^CONINFO(inewcon(c),3)).*...
        (lotchoice(inewlots,8).^CONINFO(inewcon(c),4));
    EUnomit(inewcon(c),inewlots)=(max(Income(inewcon(c))-travelcost(lotchoice(inewlots,2))-...
        Paskhouse(inewlots)-exptCnomit(inewcon(c),inewlots)',0).^CONINFO(inewcon(c),1)).*...
        (lotchoice(inewlots,4).^CONINFO(inewcon(c),2)).*(lotchoice(inewlots,3).^CONINFO(inewcon(c),3)).*...
        (lotchoice(inewlots,8).^CONINFO(inewcon(c),4));
    [Umax,mitcheck]=max([EUmit(inewcon(c),inewlots); EUnomit(inewcon(c),inewlots)]);
%     U(inewcon(c),inewlots)=Umax;
%     mitchoice(inewcon(c),inewlots)=mitcheck;
    U(inewcon(c),inewlots)=EUnomit(inewcon(c),inewlots);
    mitchoice(inewcon(c),inewlots)=2*ones(1,length(inewlots));
    ihousein=(BIDLEVEL(lotchoice(inewlots,2)).*Paskhouse(inewlots) < ...
        wtpcon(inewcon(c),inewlots)' & U(inewcon(c),inewlots)' > ...
        CONINFO(inewcon(c),1)*max(U(inewcon(c),inewlots))');
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
for nl=1:length(inewlots) 
    if BUILDTIME(lotchoice(inewlots(nl),2)) == t
        iunderbid(:,inewlots(nl))=(subPhousebid(:,inewlots(nl)) < Paskhouse...
            (inewlots(nl))/(1+discount));
    elseif BUILDTIME(lotchoice(inewlots(nl),2)) < t
        iunderbid(:,inewlots(nl))=(subPhousebid(:,inewlots(nl)) < Paskhouse...
            (inewlots(nl))*BIDLEVEL(lotchoice(inewlots(nl),2)));
    end
    subPhousebid((iunderbid(:,inewlots(nl))==1),inewlots(nl))=0;  
end

while openhouse > 0
    if isempty(find(subPhousebid > 0,1))==1
        break
    end
    wincon=zeros(1,length(inewlots));
    istillopen=find(con2lot(inewlots,2)==0);
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
        subexptcost=exptcost(:,:,conset(cs));
        iexptcost=sub2ind(size(subexptcost),mitchoice(conset(cs),...
            inewlots(ilotmatch)),inewlots(ilotmatch)');
        uset=(((Income(conset(cs))-travelcost(lotchoice(inewlots(ilotmatch),2))-...
           subPhousebid(conset(cs),inewlots(ilotmatch))'-subexptcost(iexptcost)').^CONINFO(conset(cs),1)).*...
           (lotchoice(inewlots(ilotmatch),4).^CONINFO(conset(cs),2)).*...
           (lotchoice(inewlots(ilotmatch),3).^CONINFO(conset(cs),3)).*...
           (lotchoice(inewlots(ilotmatch),8).^CONINFO(conset(cs),4)))';
%         uset=(((Income(conset(cs))-travelcost(lotchoice(inewlots(ilotmatch),2))-...
%            subPhousebid(conset(cs),inewlots(ilotmatch))').^CONINFO(conset(cs),1)).*...
%            (lotchoice(inewlots(ilotmatch),4).^CONINFO(conset(cs),2)).*...
%            (lotchoice(inewlots(ilotmatch),3).^CONINFO(conset(cs),3)).*...
%            (lotchoice(inewlots(ilotmatch),8).^CONINFO(conset(cs),4)))';

        ilotid=find(uset==max(uset));
        
        if length(ilotid) > 1
            ipick=find(subPhousebid(conset(cs),inewlots(ilotmatch(ilotid)))==...
                max(subPhousebid(conset(cs),inewlots(ilotmatch(ilotid)))));
            if length(ipick) > 1
                ipick=ipick(ceil(length(ipick)*rand(1)));
            end
            lotid=inewlots(ilotmatch(ilotid(ipick)));
        else
            lotid=inewlots(ilotmatch(ilotid));
        end
        
        conid=conset(cs);

        con2lot(lotid,1)=subPhousebid(conid,lotid);
        con2lot(lotid,2)=conid;
        con2lot(lotid,3)=lotid;
        lotchoice(lotid,7)=1;
        Lottype(Lottype(:,1)==lotid,7)=1;
        cellinfo(:,7)=Lottype(icells,7);
        con2lot(lotid,4)=max(ceil(t+avgrestime/2+normrnd(avgrestime/2,stdrestime/2,1,1)),t+1);
        RESTIME(Lottype(Lottype(:,1)==lotid,2))=con2lot(lotid,4);
        subPhousebid(conid,:)=0;
        subPhousebid(:,lotid)=0;
        subU(conid,:)=0;
        subU(:,lotid)=0;
        openhouse=openhouse-1;
        
        isamelot=(Lottype(:,1)==lotid);
        MITIGATE(Lottype(isamelot,2))=(mitchoice(conid,lotid)==1);
    end
end

% VACANT HOUSES
ifilled=find(con2lot(:,2)~=0);
istillvac=find(con2lot(:,2)==0);
popinhouse=ismember(1:length(Income),con2lot(:,2))';
ioldpop=find(popinhouse==0 & Income(:,1)~=0);

Income(con2lot(ifilled,2),2)=con2lot(ifilled,4);

for sv=1:length(istillvac)
   if BUILDTIME(lotchoice(istillvac(sv),2)) == t 
       con2lot(istillvac(sv),1)=Paskhouse(istillvac(sv))./(1+discount);
   elseif BUILDTIME(lotchoice(istillvac(sv),2)) < t
%        con2lot(istillvac(sv),1)=min((Paskhouse(istillvac(sv))*...
%            BIDLEVEL(lotchoice(istillvac(sv),2)))/((1+discount)),...
%            Paskhouse(istillvac(sv))/((1+discount)));
       con2lot(istillvac(sv),1)=max(Paskhouse(istillvac(sv))-...
           discount*ccost(lotchoice(con2lot(istillvac(sv),3),5))*...
           (t-BUILDTIME(lotchoice(istillvac(sv),2))),0);
   end
end
con2lot(istillvac,3)=istillvac;
con2lot(istillvac,4)=t+1;

for kk=1:length(iurblist)
   isamecell=ismember(Lottype(:,2),cellinfo(kk,2));
   if length(find(isamecell==1)) > 1
       RENT(cellinfo(kk,2))=mean(con2lot(Lottype(isamecell,1),1));
   elseif length(find(isamecell==1)) == 1
       RENT(cellinfo(kk,2))=con2lot(cellinfo(kk,1),1);
   end
end


%%%%%%%% Utility check %%%%%%%%%%%
realulot=zeros(length(Income),3);
conlist=(1:length(Income))';
inhouselist=conlist(popinhouse);
maxulot=zeros(1,[]);
imaxulot=zeros(1,[]);
for c=1:length(Income)
    [maxulot(c),imaxulot(c)]=max(U(c,:),[],2);   
end
for c=1:length(inhouselist)
    ireallot=find(con2lot(:,2)==inhouselist(c));
    realulot(inhouselist(c),1:3)=[ireallot lotchoice(ireallot,5) U(inhouselist(c),ireallot)];
end
maxuset=[imaxulot' lotchoice(imaxulot,5) maxulot'];
fullset=[conlist realulot maxuset];
utildiff(t)=mean(fullset(fullset(:,4)~=0,7)-fullset(fullset(:,4)~=0,4));
pctutildiff(t)=mean(fullset(fullset(:,4)~=0,4)./fullset(fullset(:,4)~=0,7));

Ufinset=zeros(length(ifilled),1);
for i=1:length(ifilled)
    conid=con2lot(ifilled(i),2);
    lotid=con2lot(ifilled(i),3);
    Ufinset(i)=(((Income(conid)-travelcost(lotchoice(ifilled(i),2))-...
        con2lot(ifilled(i),1)).^CONINFO(conid,1)).*...
        (lotchoice(ifilled(i),4).^CONINFO(conid,2)).*(lotchoice(ifilled(i),3).^...
        CONINFO(conid,3)).*(lotchoice(ifilled(i),8).^CONINFO(conid,4)))';
end
Ufinset=sort(Ufinset,'descend');
Incomeset=sort(Income(con2lot(ifilled,2)),'descend');
utilgini(t)=(length(Ufinset)+1)/(length(Ufinset)-1)-2/(length(Ufinset)*...
    (length(Ufinset)-1)*mean(Ufinset))*sum((1:length(Ufinset))'.*Ufinset);
incgini(t)=(length(Incomeset)+1)/(length(Incomeset)-1)-2/(length(Incomeset)*...
    (length(Incomeset)-1)*mean(Incomeset))*sum((1:length(Incomeset))'.*Incomeset);
