%%%%%%%%%%%%%% House Market Initial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paskhouse(1:Nlots(1),1)=INITIALPASK(iurb);
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
for c=1:length(CONINFO(:,1))
    wtpcon(c,:)=(CONINFO{c,1}-cat(1,travelcost{cat(1,lotchoice{:,2})}))*(CONINFO{c,4}+...
        CONINFO{c,5}+CONINFO{c,6});
    U(c,:)=((CONINFO{c,1}-cat(1,travelcost{cat(1,lotchoice{:,2})})-...
        Paskhouse).^CONINFO{c,3}).*(cat(1,Lottype{:,4}).^CONINFO{c,4}).*...
        (cat(1,Lottype{:,3}).^CONINFO{c,5}).*(cat(1,Lottype{:,7}).^CONINFO{c,6});
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
%     U(c,:)=EUnomit(c,:);
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
% istillvac=find(con2lot(:,2)==0);
% popout=find(ismember(1:Nconsumers,con2lot(:,2))'==0);
% popin=find((ismember(1:Nconsumers,popout)==0)==1)';
% maxbids(ifilled)=max(Phousebid(:,ifilled),[],1);
% maxbids(istillvac)=max(wtpconstar(popout,istillvac),[],1)';

% Income(con2lot(ifilled,2),2)=con2lot(ifilled,4);
% 
% con2lot(istillvac,1)=Paskhouse(istillvac)./(1+discount);
% con2lot(istillvac,3)=istillvac;
% con2lot(istillvac,4)=TSTART+1;

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
% 
% Ufinset=zeros(length(ifilled),1);
% for i=1:length(ifilled)
%     conid=con2lot(ifilled(i),2);
%     lotid=con2lot(ifilled(i),3);
%     Ufinset(i)=(((Income(conid)-travelcost(lotchoice(ifilled(i),2))-...
%         con2lot(ifilled(i),1)).^CONINFO(conid,1)).*...
%         (lotchoice(ifilled(i),4).^CONINFO(conid,2)).*(lotchoice(ifilled(i),3).^...
%         CONINFO(conid,3)).*(lotchoice(ifilled(i),8).^CONINFO(conid,4)))';
% end
% Ufinset=sort(Ufinset,'descend');
% Incomeset=sort(Income(con2lot(ifilled,2)),'descend');
% utilgini(TSTART)=(length(Ufinset)+1)/(length(Ufinset)-1)-2/(length(Ufinset)*...
%     (length(Ufinset)-1)*mean(Ufinset))*sum((1:length(Ufinset))'.*Ufinset);
% incgini(TSTART)=(length(Incomeset)+1)/(length(Incomeset)-1)-2/(length(Incomeset)*...
%     (length(Incomeset)-1)*mean(Incomeset))*sum((1:length(Incomeset))'.*Incomeset);
