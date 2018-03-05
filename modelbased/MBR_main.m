function MBR_main(buffer_cinemri3,rangey,rangex,referenceFrame,numSkip,numPreContrast_bld,numPreContrast_tiss,noi,slice,upsample_flag,upsample_factor)

temp=strcat('Output1/',int2str(slice));
mkdir(temp)
cinemri1=buffer_cinemri3(:,:,:);   
[RV LV]=FindLVRV(cinemri1);
x=floor(mean([RV(1),LV(1)]));
y=floor(mean([RV(2),LV(2)]));

[shifts_out, cinemri_least_squares] = NeighboringTemporalSmoothingRegistration(cinemri1, x-40:x+40, y-40:x+40, 15);

temp=strcat('Output1/',int2str(slice),'/cinemri1.mat');
save(temp,'cinemri1');
temp=strcat('Output1/',int2str(slice),'/cinemri_least_squares.mat');
save(temp,'cinemri_least_squares');  
    imagesc(cinemri1(:,:,referenceFrame)),colormap gray,brighten(0.3)
    temp=strcat('Output1/',int2str(slice),'/cinemri_least_squares.mat');
    load(temp)
       close all
       [bw x]=auto_roi_mbr(cinemri_least_squares);
    temp=strcat('Output1/',int2str(slice),'/bw_rv.mat');
save(temp,'bw','x','y');
g_sh=zeros(size(cinemri1,3),2,noi);

for MBR_iter_no=1:noi    
    tic
    
      temp=strcat('Output1/',int2str(slice),'/cinemri_least_squares.mat');
      load(temp)
    cinemri1=cinemri_least_squares;
    % stage 3.2 - extracting curves
    disp(strcat('Doing: extracting curves - iteration #',int2str(MBR_iter_no)))    
    ii=find(bw);    
    bldcurve=zeros(size(cinemri1,3),1);    
    parfor i=1:size(cinemri1,3)
        bldcurve(i)=sum(sum(bw.*cinemri1(:,:,i)))/length(ii);
    end
    tmpImg=ones(size(cinemri1(:,:,1)));    
    [Y,X]=find(tmpImg);    
    nX = length(X);    
    tisscurves=zeros(nX,size(cinemri1,3));    
    parfor j = 1 : nX
        tisscurves(j,:) = cinemri1(Y(j), X(j), :);
    end    
    curves=[bldcurve';tisscurves];    
    temp=strcat('Output1/',int2str(slice),'/curves.mat');
save(temp,'curves');
    % done
    
    % stage 3.4 - deltaSIcurves    
    disp(strcat('Doing: generating deltaSIcurves - iteration #',int2str(MBR_iter_no)))    
    si_curves = curves;
    [nRegs, nTimes]=size(curves);
    bldcurve=si_curves(1,:)';
    nRegs=nRegs-1;
    tisscurves = si_curves(2:nRegs+1,:);
    tisscurve=tisscurves';
    init_bld=sum(bldcurve(numSkip+1:numPreContrast_bld+numSkip))/numPreContrast_bld;
    tiss_avgSpatial=0; counter=0;
    
    init_tiss=zeros(nRegs,1);
    
    parfor ireg=1:nRegs
        init_tiss(ireg)=sum(tisscurve(numSkip+1:numPreContrast_tiss+numSkip,ireg))/numPreContrast_tiss;
        if(init_tiss(ireg)>0)
            tiss_avgSpatial=tiss_avgSpatial+init_tiss(ireg);
            counter=counter+1;
        end
    end
    tiss_avgSpatial=tiss_avgSpatial/counter;    
    
    deltaSI_bldcurve= (bldcurve - init_bld);
    deltaSI_bldcurve=tiss_avgSpatial*deltaSI_bldcurve./(tiss_avgSpatial*ones(length(bldcurve),1)) ;
    temp=strcat('Output1/',int2str(slice),'/init_bld.mat');
    save(temp,'init_bld');
    
    deltaSI_tisscurve=zeros(nTimes,nRegs);
    parfor ireg=1:nRegs
        deltaSI_tisscurve(:,ireg)=tisscurve(:,ireg)-init_tiss(ireg);
        init_tiss(ireg);
        if(init_tiss(ireg)~=0)
            deltaSI_tisscurve(:,ireg)=tiss_avgSpatial*deltaSI_tisscurve(:,ireg)./(init_tiss(ireg)*ones(length(tisscurve(:,ireg)),1)) ;
        end
    end
    
    temp=strcat('Output1/',int2str(slice),'/init_tiss.mat');
save(temp,'init_tiss');
    temp=strcat('Output1/',int2str(slice),'/tiss_avgSpatial.mat');
save(temp,'tiss_avgSpatial');
    deltaSI_curves =[ deltaSI_bldcurve'; deltaSI_tisscurve'];  
    temp=strcat('Output1/',int2str(slice),'/deltaSI_curves.mat');
save(temp,'deltaSI_curves');
    % done
    
    %  Murase fitting(compartment model fitting)
    disp(strcat('Doing: model fitting - iteration #',int2str(MBR_iter_no)))    
    murase_fitg(deltaSI_curves,slice);   %measured = fv*bld + ( conv( bld , k1*exp(-k2*(t-t0))) )    
    
    % Replacing with Murase fits
    disp(strcat('Doing: model image generation - iteration #',int2str(MBR_iter_no)))    
    cine_curv_fit(0,MBR_iter_no,Y,X,slice);
    % done   
    toc    
end
 temp=strcat('Output1/',int2str(slice),'/g_sh.mat');
save(temp,'g_sh');



