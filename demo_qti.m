clear all; close all; clc

%% read data from MRtrix files
[y,vox,v2w,grad] = MRtrix.readToMatlab('dwi.mif');
y = double(y); % dwi data, make sure it is double
grad = double(grad); % gradient directions and b-value
grad(:,4) = grad(:,4)%/1000; % express b-value in ms/um^2 for better conditioning
mask = MRtrix.readToMatlab('mask.mif'); % logical processing mask

%% smoothing image prior to fit can improve QTI fit
for i = 1:size(y,4)
    y(:,:,:,i) = smooth3(y(:,:,:,i),'gaussian',3,0.65);
end

%% set all voxels outside processing mask to NaN
y = Volumes.mask(y,mask);

%% fit model
model = QTI(grad);
x = model.solve(y);

%% calculate scalar metrics
m = model.metrics(x);

%% visual inspection of scalar metrics
figure; colormap gray; sgtitle('QTI metrics')
h(1)=subplot(3,3,1); imagesc(m.b0(:,:,49),[0 inf]); title('b0')
h(2)=subplot(3,3,2); imagesc(m.ad(:,:,49),[0 3.5]); title('ad')
h(3)=subplot(3,3,3); imagesc(m.rd(:,:,49),[0 3.5]); title('rd')
h(4)=subplot(3,3,4); imagesc(m.md(:,:,49),[0 3.5]); title('md')
h(5)=subplot(3,3,5); imagesc(m.fa(:,:,49),[0 1.0]); title('fa')
h(6)=subplot(3,3,6); imagesc(m.mki(:,:,49),[0 3.5]); title('mki')
h(7)=subplot(3,3,7); imagesc(m.mka(:,:,49),[0 3.5]); title('mka')
h(8)=subplot(3,3,8); imagesc(m.mkt(:,:,49),[0 3.5]); title('mkt')
h(9)=subplot(3,3,9); imagesc(m.ufa(:,:,49),[0 1.0]); title('ufa')
axis(h,'equal','off')

%% write metrics to MRtrix files
MRtrix.writeFromMatlab(m.b0 ,vox,v2w,[],'qti_b0.mif')
MRtrix.writeFromMatlab(m.ad ,vox,v2w,[],'qti_ad.mif')
MRtrix.writeFromMatlab(m.rd ,vox,v2w,[],'qti_rd.mif')
MRtrix.writeFromMatlab(m.md ,vox,v2w,[],'qti_md.mif')
MRtrix.writeFromMatlab(m.fa ,vox,v2w,[],'qti_fa.mif')
MRtrix.writeFromMatlab(m.mki ,vox,v2w,[],'qti_mki.mif')
MRtrix.writeFromMatlab(m.mka ,vox,v2w,[],'qti_mka.mif')
MRtrix.writeFromMatlab(m.mkt ,vox,v2w,[],'qti_mkt.mif')
MRtrix.writeFromMatlab(m.ufa ,vox,v2w,[],'qti_ufa.mif')