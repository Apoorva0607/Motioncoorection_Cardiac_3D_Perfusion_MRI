% --------------------------------------------------------------------------------
% [RegisteredData,Offsets]=jRegister2D(UnregisteredData,Mask,ReferenceData)
% Jason Mendes - UCAIR - University of Utah - 2017
% --------------------------------------------------------------------------------
% Use cross correlation to find pixel offsets.
% ----------------------------------------------------------------------
function [RegisteredData,Offsets]=jRegister3D(UnregisteredData,Mask,ReferenceData)

% Get data from the selected region
NumberOfIterations=4;
MaxShiftPerIteration=3;
SearchWindowX=8;
SearchWindowZ=1;
[Nx,Ny,Nz,Nm]=size(UnregisteredData);
KspaceFilter=repmat(reshape(gausswin(Nx,10),Nx,1),[1 Ny Nz]).*repmat(reshape(gausswin(Ny,10),1,Ny),[Nx 1 Nz]).*repmat(reshape(gausswin(Nz,10),1,1,Nz),[Nx Ny 1]);
SearchX=(-SearchWindowX:SearchWindowX)+floor(Nx/2)+1;
SearchZ=(-SearchWindowZ:SearchWindowZ)+floor(Nz/2)+1;
Offsets=[];
Offsets.X=zeros(NumberOfIterations,Nm);
Offsets.Y=zeros(NumberOfIterations,Nm);
Offsets.Z=zeros(NumberOfIterations,Nm);
if (exist('RefernceData','var')~=1)
    BaseData=UnregisteredData;
else
    BaseData=ReferenceData;
end
PhaseX=repmat(reshape(-2i*pi*((1:Nx)-floor(Nx/2)-1)/Nx,Nx,1),[1 Ny Nz]);
PhaseY=repmat(reshape(-2i*pi*((1:Ny)-floor(Ny/2)-1)/Ny,1,Ny),[Nx 1 Nz]);
PhaseZ=repmat(reshape(-2i*pi*((1:Nz)-floor(Nz/2)-1)/Nz,1,1,Nz),[Nx Ny 1]);
RegisteredData=UnregisteredData;
ImageMaskX=1.0*(conv2(single(Mask),ones(SearchWindowX,SearchWindowX)/(SearchWindowX^2),'same')>0);
ImageMaskX=conv2(ImageMaskX,ones(SearchWindowX,SearchWindowX)/(SearchWindowX^2),'same');
ImageMaskZ=triang(Nz);
ImageMaskZ=min(min(ImageMaskZ(SearchZ)),ImageMaskZ);
ImageMaskZ=reshape(ImageMaskZ/max(ImageMaskZ),1,1,Nz);
ImageMask=repmat(ImageMaskZ,[Nx Ny 1]).*repmat(ImageMaskX,[1 1 Nz]);
for MeasIndex=Nm:-1:1
    for IterIndex=1:NumberOfIterations
        CurrentImage=jFFT(ImageMask.*BaseData(:,:,:,MeasIndex));
        ReferenceImage=jFFT(ImageMask.*BaseData(:,:,:,min(MeasIndex+1,Nm)));
        FilteredCrossCorrelation=abs(jIFFT(exp(1i*angle(ReferenceImage.*conj(CurrentImage))).*KspaceFilter));
        ReducedCrossCorrelation=FilteredCrossCorrelation(SearchX,SearchX,SearchZ);
        PeakLocation=find(ReducedCrossCorrelation(:)==max(ReducedCrossCorrelation(:)),1);
        [IndexX,IndexY,IndexZ]=ind2sub(size(ReducedCrossCorrelation),PeakLocation);
        IndexX=min(max(mean(IndexX)-floor(length(SearchX)/2)+floor(Nx/2),2),Nx-1);
        IndexY=min(max(mean(IndexY)-floor(length(SearchX)/2)+floor(Ny/2),2),Ny-1);
        IndexZ=min(max(mean(IndexZ)-floor(length(SearchZ)/2)+floor(Nz/2),2),Nz-1);
        LocalPeak=sum(sum(FilteredCrossCorrelation((-1:1)+IndexX,(-1:1)+IndexY,(-1:1)+IndexZ),2),3);
        Offsets.X(IterIndex,MeasIndex)=max(min(IndexX-0.5*(LocalPeak(3)-LocalPeak(1))/(LocalPeak(1)+LocalPeak(3)-2*LocalPeak(2))-floor(Nx/2)-1,MaxShiftPerIteration),-MaxShiftPerIteration);
        LocalPeak=sum(sum(FilteredCrossCorrelation((-1:1)+IndexX,(-1:1)+IndexY,(-1:1)+IndexZ),1),3);
        Offsets.Y(IterIndex,MeasIndex)=max(min(IndexY-0.5*(LocalPeak(3)-LocalPeak(1))/(LocalPeak(1)+LocalPeak(3)-2*LocalPeak(2))-floor(Ny/2)-1,MaxShiftPerIteration),-MaxShiftPerIteration);
        LocalPeak=sum(sum(FilteredCrossCorrelation((-1:1)+IndexX,(-1:1)+IndexY,(-1:1)+IndexZ),1),2);
        Offsets.Z(IterIndex,MeasIndex)=max(min(IndexZ-0.5*(LocalPeak(3)-LocalPeak(1))/(LocalPeak(1)+LocalPeak(3)-2*LocalPeak(2))-floor(Nz/2)-1,MaxShiftPerIteration),-MaxShiftPerIteration);       
        BaseData(:,:,:,MeasIndex)=abs(jIFFT(jFFT(BaseData(:,:,:,MeasIndex)).*exp(PhaseX.*Offsets.X(IterIndex,MeasIndex)+PhaseY.*Offsets.Y(IterIndex,MeasIndex)+PhaseZ.*Offsets.Z(IterIndex,MeasIndex))));
        RegisteredData(:,:,:,MeasIndex)=abs(jIFFT(jFFT(RegisteredData(:,:,:,MeasIndex)).*exp(PhaseX.*Offsets.X(IterIndex,MeasIndex)+PhaseY.*Offsets.Y(IterIndex,MeasIndex)+PhaseZ.*Offsets.Z(IterIndex,MeasIndex))));   
    end
end