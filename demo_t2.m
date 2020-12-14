clear all; close all; clc

%% read data from MRtrix files
[y,vox,v2w,te] = MRtrix.readToMatlab('t2.mif');
y = double(y); % t2 data, make sure it is double
te = double(te); % TE values
mask = MRtrix.readToMatlab('mask.mif'); % logical processing mask

%% smoothing image prior to fit can improve T2 fit
for i = 1:size(y,4)
    y(:,:,:,i) = smooth3(y(:,:,:,i),'gaussian',3,0.65);
end

%% set all voxels outside processing mask to NaN
y = Volumes.mask(y,mask);

%% fit model
model = T2(te);
x = model.solve(y);

%% calculate scalar metrics
m = model.metrics(x);

%% visual inspection of scalar metrics
figure; colormap gray; sgtitle('T2 metrics')
h(1)=subplot(3,3,1); imagesc(m.rho(:,:,11),[0 inf]); title('rho')
h(2)=subplot(3,3,2); imagesc(m.t2(:,:,11),[0 3.5]); title('t2')
axis(h,'equal','off')

%% write metrics to MRtrix files
MRtrix.writeFromMatlab(m.rho ,vox,v2w,[],'t2_rho.mif')
MRtrix.writeFromMatlab(m.t2 ,vox,v2w,[],'t2_t2.mif')
