%%%%%%%%%%%%%% House Market %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
for c=1:length(inewcon)
    if isempty(find(inewlots,1))==1
        break
    end
    wtpcon(inewcon(c),inewlots)=(CONINFO{inewcon(c),1}-cat(1,travelcost{cat(1,lotchoice{inewlots,2})}))*...
        (CONINFO{inewcon(c),4}+CONINFO{inewcon(c),5}+CONINFO{inewcon(c),6});
    U(inewcon(c),inewlots)=((CONINFO{inewcon(c),1}-cat(1,travelcost{cat(1,lotchoice{inewlots,2})})-...
        Paskhouse(inewlots)).^CONINFO{inewcon(c),3}).*(cat(1,Lottype{inewlots,4}).^...
        CONINFO{inewcon(c),4}).*(cat(1,Lottype{inewlots,3}).^CONINFO{inewcon(c),5}).*...
        (cat(1,Lottype{inewlots,7}).^CONINFO{inewcon(c),6});
    %%%%%%%%%% Expected Costs %%%%%%%%%%%%
%     subdamcoef=damcoef{inewcon(c),t};
%     exptCmit(inewcon(c),inewlots)=subdamcoef(inewlots).*cat(1,Pflood{cat(1,lotchoice{inewlots,2})}).*...
%         (Cmit+miteff*Cdam*Paskhouse(inewlots))+(1-subdamcoef(inewlots)).*(1-cat(1,Pflood{cat(1,lotchoice{inewlots,2})}))*Cmit;
%     exptCnomit(inewcon(c),inewlots)=subdamcoef(inewlots).*cat(1,Pflood{cat(1,lotchoice{inewlots,2})}).*...
%         (Cdam*Paskhouse(inewlots));
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
%     U(inewcon(c),inewlots)=EUnomit(inewcon(c),inewlots);
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
% utilgini(t)=(length(Ufinset)+1)/(length(Ufinset)-1)-2/(length(Ufinset)*...
%     (length(Ufinset)-1)*mean(Ufinset))*sum((1:length(Ufinset))'.*Ufinset);
% incgini(t)=(length(Incomeset)+1)/(length(Incomeset)-1)-2/(length(Incomeset)*...
%     (length(Incomeset)-1)*mean(Incomeset))*sum((1:length(Incomeset))'.*Incomeset);
