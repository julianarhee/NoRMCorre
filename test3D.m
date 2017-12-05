clear
gcp;

source = '/nas/volume1/2photon/RESDATA/test_motion_correction_3D/DATA/Raw';
output_path = fullfile(source, 'NoRMCorre', 'results');
if ~exist(output_path)
    mkdir(output_path)
end
filename = 'fov1_bar037Hz_run4_File001.tif';
name = fullfile(source, filename);

memmap = false %true;
nonrigid = false;

tic; Y = read_file(name); toc; % read the file (optional, you can also pass the path in the function instead of Y)

nslices = 12;
Y = Y(:,:,1:2:end);
Yr = reshape(Y, [size(Y,1), size(Y,2), nslices, size(Y,3)/nslices]);
size(Yr)
d1 = size(Yr,1);
d2 = size(Yr,2);
d3 = size(Yr,3);

if memmap
    Yr = single(Yr);                 % convert to single precision 
    T = size(Yr,ndims(Yr));
    Y = matfile(fullfile(output_path, 'Y_memmap.mat'), 'Writable', true);
    Y.Y = Yr;
    clear Yr;
else
    Y = single(Yr);
    T = size(Y, ndims(Y))
end
%Y = Y - min(Y(:));
%% set parameters (first try out rigid motion correction)

paramname = 'shift30dev5';

%options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',50,'max_shift',15,'us_fac',50);
%options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'d3',size(Y,3),'bin_width',50,'max_shift',30,'us_fac',50,...
                                    %'correct_bidir', false, 'max_dev', 5);
if memmap
    options_rigid = NoRMCorreSetParms('d1',d1,'d2',d2,'d3',d3,'bin_width',50,'max_shift',30,'us_fac',50,...
                                    'correct_bidir', false, 'max_dev', 5);
else
    options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'d3',size(Y,3),'bin_width',50,'max_shift',30,'us_fac',50,'correct_bidir', false, 'max_dev', 5);

end

%% perform motion correction
tic; [M1,shifts1,template1] = normcorre(Y,options_rigid); toc

%% now try non-rigid motion correction (also in parallel)
%options_nonrigid =
%NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[32,32],'mot_uf',4,'bin_width',50,'max_shift',15,'max_dev',3,'us_fac',50);

if nonrigid
if memmap
    options_nonrigid = NoRMCorreSetParms('d1',d1,'d2',d2,'d3',d3,...
                                    'grid_size',[32,32],'mot_uf',4,'bin_width',50,...
                                    'max_shift',30,'max_dev',5,'us_fac',50,'correct_bidir',false);
else
    options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'d3',size(Y,3),...
                                    'grid_size',[32,32],'mot_uf',4,'bin_width',50,...
                                    'max_shift',30,'max_dev',5,'us_fac',50,'correct_bidir',false);
end

tic; [M2,shifts2,template2] = normcorre_batch(Y,options_nonrigid); toc
end

%% compute metrics
get_metrics = false

if get_metrics
nnY = quantile(Y(:),0.005);
mmY = quantile(Y(:),0.995);

[cY,mY,vY] = motion_metrics(Y,10);
[cM1,mM1,vM1] = motion_metrics(M1,10);
if nonrigid
    [cM2,mM2,vM2] = motion_metrics(M2,10);
end
T = length(cY);
%% plot metrics
figure;
    ax1 = subplot(2,3,1); imagesc(mY,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean raw data','fontsize',14,'fontweight','bold')
    ax2 = subplot(2,3,2); imagesc(mM1,[nnY,mmY]);  axis equal; axis tight; axis off; title('mean rigid corrected','fontsize',14,'fontweight','bold')
    if nonrigid
        ax3 = subplot(2,3,3); imagesc(mM2,[nnY,mmY]); axis equal; axis tight; axis off; title('mean non-rigid corrected','fontsize',14,'fontweight','bold')
    end
    if nonrigid
        subplot(2,3,4); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
    else
        subplot(2,3,4); plot(1:T,cY,1:T,cM1); legend('raw data','rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
    end

    subplot(2,3,5); scatter(cY,cM1); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('raw data','fontsize',14,'fontweight','bold'); ylabel('rigid corrected','fontsize',14,'fontweight','bold');
    if nonrigid
        subplot(2,3,6); scatter(cM1,cM2); hold on; plot([0.9*min(cY),1.05*max(cM1)],[0.9*min(cY),1.05*max(cM1)],'--r'); axis square;
        xlabel('rigid corrected','fontsize',14,'fontweight','bold'); ylabel('non-rigid corrected','fontsize',14,'fontweight','bold');
        linkaxes([ax1,ax2,ax3],'xy')
     
    end
        savefig(fullfile(output_path, sprintf('metrics_%s_%s', filename(1:end-4), paramname)));
close;
end

%% plot shifts        

shifts_r = squeeze(cat(3,shifts1(:).shifts));
if nonrigid
shifts_nr = cat(ndims(shifts2(1).shifts)+1,shifts2(:).shifts);
shifts_nr = reshape(shifts_nr,[],ndims(Y)-1,T);

shifts_x = squeeze(shifts_nr(:,1,:))';
shifts_y = squeeze(shifts_nr(:,2,:))';

patch_id = 1:size(shifts_x,2);
str = strtrim(cellstr(int2str(patch_id.')));
str = cellfun(@(x) ['patch # ',x],str,'un',0);

figure;
    ax1 = subplot(311); plot(1:T,cY,1:T,cM1,1:T,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax2 = subplot(312); plot(shifts_x); hold on; plot(shifts_r(:,1),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
            set(gca,'Xtick',[])
    ax3 = subplot(313); plot(shifts_y); hold on; plot(shifts_r(:,2),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
            xlabel('timestep','fontsize',14,'fontweight','bold')
    linkaxes([ax1,ax2,ax3],'x')

savefig(fullfile(output_path, sprintf('shifts_%s_%s', filename(1:end-4), paramname)));
close;
end


%% plot a movie with the results

figure;
for t = 1:1:T
    subplot(121);imagesc(Y(:,:,t),[nnY,mmY]); xlabel('raw data','fontsize',14,'fontweight','bold'); axis equal; axis tight;
    title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    if nonrigid
        subplot(122);imagesc(M2(:,:,t),[nnY,mmY]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
        title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    else
        subplot(122);imagesc(M2(:,:,t),[nnY,mmY]); xlabel('non-rigid corrected','fontsize',14,'fontweight','bold'); axis equal; axis tight;
        title(sprintf('Frame %i out of %i',t,T),'fontweight','bold','fontsize',14); colormap('bone')
    end
    set(gca,'XTick',[],'YTick',[]);
    drawnow;
    pause(0.02);
end
if nonrigid
    tiffWrite(M2, sprintf('nonrigid_%s_%s', filename, paramname), output_path);
else
    tiffWrite(M1, sprintf('rigid_%s_%s', filename, paramname), output_path);
end
close;
