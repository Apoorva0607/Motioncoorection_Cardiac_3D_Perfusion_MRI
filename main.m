clc;
close all;
clear ;
step=0.05;
% Load the file you wish to register , although it might happen the file
% you load might have a different name , so be careful with the filenames.
load('RegisteredData_P031417_MID114.mat');
%% Rigid Registration
UnregisteredData=Data.Systole;
% To select a mask around the myocardium
figure(20);imagesc(UnregisteredData(:,:,4,32));
Mask=roipoly;
%
[RegisteredRigidData,Offsets]=jRegister3D(UnregisteredData,Mask);

%% Model images
mkdir Output1
buffer_cinemri=abs((RegisteredRigidData(:,:,:,:)));
for i=1:8
    referenceFrame=20;
    numPreContrast_bld=5;
    numPreContrast_tiss=7;
    numSkip=1;
    pd_end=10;
    seriesNum=1000*i;
    %%%%%%
    buffer_cinemri1(:,:,:)=buffer_cinemri(:,:,i,:);
    [ sx sy st]=size(buffer_cinemri1);
    rangex=1:sx;
    rangey=1:sy;
    
    %%%Model based registration ( cinemri_curv_fit1.mat is the model image for
    %%%every slice)
    run_registration(buffer_cinemri1(:,:,pd_end+1:end),rangey,rangex,referenceFrame,  numSkip,numPreContrast_bld,numPreContrast_tiss,seriesNum);
end
%%
cd('Output')
% Creating reference volume of model images from all slices.This is done because the model based registration is done after the
% proton density frames.
model3d_images(:,:,:,1:pd_end)=RegisteredRigidData(:,:,:,1:pd_end);

for i=1:8
    seriesNum=num2str(1000*i);
    cd(seriesNum);
    load('cinemri_curv_fit1.mat');
    
    model3d_images(:,:,i,pd_end+1:size(RegisteredRigidData,4))=cinemri_curv_fit;
    cd ..
end
cd ..
%% Deformable registration
I0=RegisteredRigidData;
Final_registered=zeros(size(I0,1),size(I0,2),size(I0,3),size(I0,4));
X=zeros(size(I0));
Y=zeros(size(I0));
Z=zeros(size(I0));

tstart=tic;
parfor i=pd_end+1:size(I0,4)
    
    [I2,x,y,z,E] = Reg_GS_3D(I0(:,:,:,i),model3d_images(:,:,:,i),step);
    
    Final_registered(:,:,:,i)=I2;
    X(:,:,:,i)=x;
    Y(:,:,:,i)=y;
    Z(:,:,:,i)=z;
end
telapsed=toc(tstart);
%%
% Results:
for i=1:size(I0,4)
    figure(30);imagesc(squeeze(Final_registered(:,:,4,i)));colormap gray;colorbar;title(['time frame ',num2str(i)]);
    figure(40);imagesc(squeeze(Data.Systole(:,:,4,i)));colormap gray;colorbar;title(['time frame ',num2str(i)]);
    pause(0.2);
end


