clear all; close all; clc

%% read data from MRtrix files
[y,vox,v2w,grad] = MRtrix.readToMatlab('dwi.mif');
y = double(y); % dwi data, make sure it is double
grad = double(grad); % gradient directions and b-value
grad(:,4) = grad(:,4)/1000; % express b-value in ms/um^2 for better conditioning
mask = MRtrix.readToMatlab('mask.mif'); % logical processing mask

%% smoothing image prior to fit can improve DTI fit
for i = 1:size(y,4)
    y(:,:,:,i) = smooth3(y(:,:,:,i),'gaussian',3,0.65);
end

%% set all voxels outside processing mask to NaN
y = Volumes.mask(y,mask);

%% fit model
model = DTI(grad);
x = model.solve(y);

%% calculate scalar metrics
m = model.metrics(x);

%% visual inspection of scalar metrics
figure; colormap gray; sgtitle('DTI metrics')
h(1)=subplot(3,3,1); imagesc(m.b0(:,:,11),[0 inf]); title('b0')
h(2)=subplot(3,3,2); imagesc(m.ad(:,:,11),[0 3.5]); title('ad')
h(3)=subplot(3,3,3); imagesc(m.rd(:,:,11),[0 3.5]); title('rd')
h(4)=subplot(3,3,4); imagesc(m.md(:,:,11),[0 3.5]); title('md')
h(5)=subplot(3,3,5); imagesc(m.fa(:,:,11),[0 1.0]); title('fa')
axis(h,'equal','off')

%% write metrics to MRtrix files
MRtrix.writeFromMatlab(m.b0 ,vox,v2w,[],'dki_b0.mif')
MRtrix.writeFromMatlab(m.ad ,vox,v2w,[],'dki_ad.mif')
MRtrix.writeFromMatlab(m.rd ,vox,v2w,[],'dki_rd.mif')
MRtrix.writeFromMatlab(m.md ,vox,v2w,[],'dki_md.mif')
MRtrix.writeFromMatlab(m.fa ,vox,v2w,[],'dki_fa.mif')
