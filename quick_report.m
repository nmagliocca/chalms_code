
figure(1)
imagesc(reshape(LOTTYPE(:,15),80,80))
% surf(reshape(LOTTYPE(:,25),80,80))
axis ij
view(0,90)
set(gca,'clim',[1 HT])
colorbar
title('t=15')

figure(2)
plot(10:30,avgrent(:,10:30),'-')
legend('HT 1','HT 2','HT 3','HT 4','HT 5','HT 6','HT 7','HT 8','HT 9')
ylabel('Rent')
xlabel('Time Step')
title('Housing Prices')
xlabel('Time Step')
ylabel('Rent')

incometest=zeros(80,80);
preftest=zeros(80,80);
houseprice=zeros(80,80);
inotvac=find(cat(1,lotchoice{:,4})==1);
for i=1:length(cat(1,Lottype{inotvac,1}))
    incometest(Lottype{inotvac(i),2})=CONINFO{lotchoice{inotvac(i),5},1}*...
        ones(length(Lottype{inotvac(i),2}),1);
    preftest(Lottype{inotvac(i),2})=CONINFO{lotchoice{inotvac(i),5},6}*...
        ones(length(Lottype{inotvac(i),2}),1);
    houseprice(Lottype{inotvac(i),2})=lotchoice{inotvac(i),7};
end

figure(3)
imagesc(incometest)
% surf(incometest)
axis ij
view(0,90)
set(gca,'clim',[0 200000])
colorbar
title('Consumer Income')

figure(4)
imagesc(preftest)
% surf(preftest)
axis ij
view(0,90)
set(gca,'clim',[0 0.5])
colorbar
title('Consumer Coastal Pref')

figure(5)
imagesc(reshape(BUILDTIME,80,80))
% surf(reshape(BUILDTIME,80,80))
axis ij
view(0,90)
set(gca,'clim',[10 30])
colorbar
title('Build Time')

figure(6)
% surf(houseprice)
imagesc(houseprice)
axis ij
view(0,90)
colorbar
title('House Prices at t=30')

landsales=zeros(80,80);
for nf=1:length(farmsoldinfo(:,1))
    ifarm=find(cat(1,LANDINFO{1,10})==farmsoldinfo(nf,1));
    landsales(ifarm)=farmsoldinfo(nf,4);
end
figure(7)
% surf(landsales)
imagesc(landsales)
axis ij
view(0,90)
colorbar
title('Land Sales')

%%%% Viz of farmer expectations
fexptmap=zeros(NLENGTH,NWIDTH,length(TSTART+1:TMAX));
fsoldmap=zeros(NLENGTH,NWIDTH,length(TSTART+1:TMAX));
for it=TSTART+1:TMAX
   ifarmer=find(landprojSAVE(:,it)~=0 & isnan(landprojSAVE(:,it))==0);
   subfexptmap=zeros(NLENGTH,NWIDTH);
   for ii=1:length(ifarmer)
       nfinfo=cat(2,Farminfo{ifarmer(ii),2});
       subfexptmap(nfinfo(:,1))=landprojSAVE(ifarmer(ii),it);
   end
   fexptmap(:,:,it-TSTART)=subfexptmap;
   ifsold=find(farmsoldinfo(:,2)==it);
   if it>1
           subfsoldmap=fsoldmap(:,:,it-TSTART);
       else
           subfsoldmap=fsoldmap(:,:,1);
   end
   if isempty(find(ifsold,1))==1
       fsoldmap(:,:,it-TSTART+1)=subfsoldmap;
       continue
   end
   for j=1:length(ifsold)
       sdinfo=cat(2,Farminfo{farmsoldinfo(ifsold(j),1),2});
       
%        subfsoldmap=fsoldmap(:,:,farmsoldinfo(ifsold(j),2)-TSTART);
       subfsoldmap(sdinfo(:,1))=farmsoldinfo(ifsold(j),4); 
   end
   fsoldmap(:,:,it-TSTART+1)=subfsoldmap;
end

lpsort=sortrows(farmsoldinfo,2);

writerObj=VideoWriter('fexpt_viz2.mp4','MPEG-4');
writerObj.FrameRate=2;
open(writerObj);

figure(8)
maxpoint=max(max(farmsoldinfo(:,4)),max(max(landprojSAVE)));
clim=[0 maxpoint];
nframes=length(TSTART+1:TMAX);
F(1:nframes)=struct('cdata',[],'colormap',[]);
for tt=1:length(TSTART+1:TMAX)
    subplot(3,1,1)
    axis([0 80 0 80 0 maxpoint]);
    imagesc(fexptmap(:,:,tt),clim)
    h=get(gcf,'CurrentAxes');
    set(h,'clim',clim)
    colorbar
    title('Farmer Price Expectations')
    drawnow
    hold on
    subplot(3,1,2)
    axis([0 80 0 80 0 maxpoint]);
    imagesc(fsoldmap(:,:,tt),clim)
    h=get(gcf,'CurrentAxes');
    set(h,'clim',clim)
    colorbar
    title('Land Sale Price')
    drawnow
    hold on
    subplot(3,1,3)
    axis([11 30 0 maxpoint]);
    set(gca,'nextplot','replacechildren')
    iplot=find(lpsort(:,2)<tt+TSTART,1,'last');
    plot(lpsort(1:iplot,2),lpsort(1:iplot,4),'*')
    drawnow
    hold on
%     ylabel('Land Price')
%     xlabel('Time Step')
    title('Land Price vs. Time')
%     frm=figure(8);
%     F(tt)=getframe(frm);
%     writeVideo(writerObj,F(tt));
end
close(writerObj);

movie2avi(F,'farmer_expt_viz.avi','fps',2);

% M=movie(F,4,2);
writerObj=VideoWriter('farmer_viz.mp4');
open(writerObj);
writeVideo(writerObj,F)


%%% Check Hedonics
% mdl=fitlm([cat(1,Lottype{ifilled,3}) cat(1,Lottype{ifilled,4})...
%         cat(1,Lottype{ifilled,8}) cat(1,Lottype{ifilled,7})],...
%         cat(1,lotchoice{ifilled,7}))
% try adding consumer attributes, which can be replaced by broker averages
% for prediction
mdl=fitlm([cat(1,Lottype{ifilled,3}) cat(1,Lottype{ifilled,4})...
        cat(1,Lottype{ifilled,8}) cat(1,Lottype{ifilled,7}) ...
        cat(1,CONINFO{cat(1,lotchoice{ifilled,5}),1})],...
        cat(1,lotchoice{ifilled,7}));
rentmdl=fitlm([cat(1,Lottype{ifilled,3}) ...
        cat(1,Lottype{ifilled,8}) cat(1,Lottype{ifilled,7}) ...
        cat(1,CONINFO{cat(1,lotchoice{ifilled,5}),1})],...
        cat(1,lotchoice{ifilled,7}));   

figure
plot(rentmdl)