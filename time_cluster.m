clear

% detrend_n=1; FWHM = 2; cluster_size_limit = 5; winsz = 7;
detrend_n=1; FWHM =1; cluster_size_limit =4; winsz =3.5;

zscore_limit =100;

close all
% data = load('Mxy_100inhom_800TR_phantom.mat').Echo_Mxy;
data = load('phantom_nonhom.mat').Echo_Mxy;

data = data(:,:,:,761:end);
data(isnan(data)) = 0;
data = squeeze(data);

mask = ones(size(data,1),size(data,2));
[row,col,tp,repet] = size(data);
data = reshape(data,[row col tp repet]);
%%

% TR =0.005;
% Ts = TR;
% fs = 1/Ts;
% f = (-40/2:40/2-1)*fs/40;
% 
% stop_data=[];
% tic;
% for i =1:row
%     for j=1:col
%         for k =1:repet
%         stop_data(i,j,:,k) = filtering_signal(squeeze(data(i,j,:,k)));
%         end
%     end
%     a = toc;
%     b= sprintf('%.2d 초 경과, i=%d',a,i);
%     disp(b)
%     tic;
% end
% data = stop_data;
%%
gaussian_data = [];
for a=1:row
    for b=1:col
            for i=1:repet
                gaussian_data(a,b,:,i)=tgaussmooth(squeeze(data(a,b,:,i)),1,'omitnan',winsz,1);
                % gaussian_data(a,b,:,i)=smoothdata(squeeze(data(a,b,:,i)),1,'gaussian',winsz);
            end
    end
end
blur_data = [];
for i=1:tp
    for j=1:repet
            blur_data(:,:,i,j)=imgaussfilt(squeeze(gaussian_data(:,:,i,j)),FWHM);
    end
end

blur_data = mean(blur_data,4);

% gaussian_data = mean(gaussian_data,4);
psc_data = [];
for i=1:row
    for j=1:col
        a = squeeze(blur_data(i,j,:))';
        m=mean(a(:,1:10),2);
        % m=mean(a(:,1:20),2);
        psc=100*(a-repmat(m,1,tp))./repmat(m,1,tp);
        psc_data(i,j,:) = psc;
    end
end


detrend_data=[];
for i=1:row
    for j=1:col
            dt=detrend(squeeze(psc_data(i,j,:)),detrend_n)';
            detrend_data(i,j,:)=dt-repmat(mean(dt(:,1:10),2),1,tp);
    end
end

%%
detrend_data = repmat(mask,[1 1 tp]).*detrend_data;
detrend_data((detrend_data)==0) = NaN;

mean_data = detrend_data;
mean_value = mean(mean_data( ~isnan(mean_data) ));
std_value = std(mean_data( ~isnan(mean_data) ));
zscore_val = (mean_data-mean_value)/std_value;

% outlier_voxels = sum((squeeze(sum((zscore_val < -zscore_limit) | (zscore_val >  zscore_limit), 3)) ~= 0),3);

outlier_voxels = logical(sum((squeeze(sum((zscore_val < -zscore_limit) | (zscore_val >  zscore_limit), 3)) ~= 0),3));

template = mean(data,[3 4]);
red = cat(3,ones(size(template)),zeros(size(template)),zeros(size(template)));

figure;
imshow(template,[]);
% hold on
h1 = imshow(red);
set(h1,'Alphadata',(outlier_voxels))
%%
detrend_data(repmat(outlier_voxels,[1 1 tp])) =NaN ;
[M,I] = max(detrend_data(:,:,11:end),[],3);
% [M,I] = min(detrend_data(:,:,11:40),[],3);



index_value = I.*mask;
index_value = index_value.*(~outlier_voxels);

figure;
ax1 = subplot(1,1,1);
imagesc(template); colormap(ax1,'gray'); axis(ax1,'image'); 
set(gca,'XTickLabel',[],'YTickLabel',[],'YTick',[],'XTick',[]);
ax2 = axes; imagesc(ax2,index_value*5,'Alphadata',index_value>0); colormap(ax2,'jet'); 
axis(ax2,'image'); ax2.Visible = 'off'; linkprop([ax1 ax2],'Position'); linkaxes([ax1 ax2],'xy')
set(gcf,'color',[1 1 1])

% figure;
for i = 1:tp-10
    iv(:,:,i) = (index_value == i);
    % subplot(5,6,i);imagesc(flip(rot90(iv(:,:,i)),1));
    % title(strcat(num2str(i*5-5),'ms'));
end
for i = 1:tp-10
L(:,:,i) = bwlabel(iv(:,:,i),4);
end

%%
for t = 1:tp-10
    label = L(:,:,t);
for i = 1:max(label(:))
    if (size(find(label==i),1)<=cluster_size_limit)
    label(label==i)=0;
    end
end
L(:,:,t) = label;
end
figure;imagesc(sum(L,3));


LL = (logical(L));
for i=1:tp-10
    ll = LL(:,:,i);
    L_time(:,:,i) = ll*5*i;
end
figure;
ax1 = subplot(1,1,1);
imagesc(template); colormap(ax1,'gray'); axis(ax1,'image'); 
set(gca,'XTickLabel',[],'YTickLabel',[],'YTick',[],'XTick',[]);
ax2 = axes; imagesc(ax2,sum(L_time,3),'Alphadata',sum(L_time,3)>0); colormap(ax2,'jet'); 
axis(ax2,'image'); ax2.Visible = 'off'; linkprop([ax1 ax2],'Position'); linkaxes([ax1 ax2],'xy')
set(gcf,'color',[1 1 1])




figure;
for i = 1:tp-10
ax1 = subplot(5,6,i);
imagesc(template); colormap(ax1,'gray'); axis(ax1,'image'); 
set(gca,'XTickLabel',[],'YTickLabel',[],'YTick',[],'XTick',[]);
title(strcat(num2str(i*5),' ms'));

ax2 = axes; imagesc(ax2,L(:,:,i),'Alphadata',L(:,:,i)>0); colormap(ax2); 
axis(ax2,'image'); ax2.Visible = 'off'; linkprop([ax1 ax2],'Position'); linkaxes([ax1 ax2],'xy')
set(gcf,'color',[1 1 1])
end

%% plot
color ="#05a61f"; %green
figure;
for i = 1:25
subplot(5,5,i);

roi = zeros(size(data,1),size(data,2));
roi(:,i) = 1;
% roi(i,:) = 1;

ss =gaussian_data(:,:,:,:).*roi;
sss = squeeze(sum(ss,[1 2]))/size(find(roi),1);
psc=detrend(100*(sss-mean(sss(1:10,:)))./mean(sss(1:10,:)),detrend_n);
psc=psc-repmat(mean(psc(1:10,:),1),40,1);
plot(-45:5:150,mean(psc,2),'LineWidth',2.5,'color',color);
title(strcat('PE direction:',num2str(i),'st line'));

end

figure;
for i = 26:50
subplot(5,5,i-25);

roi = zeros(size(data,1),size(data,2));
roi(:,i) = 1;
% roi(i,:) = 1;

ss =gaussian_data(:,:,:,:).*roi;
sss = squeeze(sum(ss,[1 2]))/size(find(roi),1);
psc=detrend(100*(sss-mean(sss(1:10,:)))./mean(sss(1:10,:)),detrend_n);
psc=psc-repmat(mean(psc(1:10,:),1),40,1);
plot(-45:5:150,mean(psc,2),'LineWidth',2.5,'color',color);
title(strcat('PE direction:',num2str(i),'st line'));

end




%% plot data
dd = data;
roi = (L(:,:,2)==5); % roi = load('s1a_roi.mat').roi;
color ="#05a61f"; %green
ss =gaussian_data(:,:,:,:).*roi;
sss = squeeze(sum(ss,[1 2]))/size(find(roi),1);
psc=detrend(100*(sss-mean(sss(1:10,:)))./mean(sss(1:10,:)),detrend_n);
psc=psc-repmat(mean(psc(1:10,:),1),40,1);
 % figure;plot(-45:5:150,psc,'LineWidth',2.5,'color',color);


figure;plot(-45:5:150,mean(psc,2),'LineWidth',2.5,'color',color);
hold on
   error = std(psc(:,:),0,2)./sqrt(size(psc(:,:),2));
    errorbar(-45:5:150,mean(psc,2),error,'.b', 'LineWidth',0.7, 'MarkerSize',0.7,'color',color);
xlabel('Time [ms]');ylabel('Percent signal change [%]');xline(0,'r');yline(0,'--k');box off; fontsize(15,"points");


