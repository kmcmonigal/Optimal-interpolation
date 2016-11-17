% optimally interpolate our fake, subsampled grid
% keep in mind that vertical is NOT isotropic.
% temperature

fake_grid.temp(250,8)=NaN;
fake_grid.temp(300,8)=NaN;
fake_grid.sal(250,8)=NaN;
fake_grid.sal(300,8)=NaN;
fake_grid.theta(250,8)=NaN;
fake_grid.theta(300,8)=NaN;

small_scale_ratio=.20; %.01, .15 (and twice)
clear x
clear xgrid
clear p
clear pgrid

% basically do OI twice: once with a Gaussian with large decay scales, once
% with "real statistics" of smaller scale
fake_grid.pres=nan(420,9);
for i=1:420
    for j=1:9
        fake_grid.pres(i,j)=10*i;
    end
end

c2_pos = [27.5167,-33.4232];
c3_pos = [27.5698,-33.5112];
a_pos = [27.595,-33.5583];
b_pos = [27.6428,-33.6674];
c_pos = [27.7152,-33.7996];
d_pos = [27.8603,-34.0435];
e_pos = [28.0316,-34.29];
f_pos = [28.17,-34.5380];
g_pos = [28.3453,-34.8213];

instr_pos = [c2_pos;c3_pos;a_pos;b_pos;c_pos;d_pos;e_pos;f_pos;g_pos];
x=0:500:1000*sw_dist([instr_pos(1,2) instr_pos(9,2)],[instr_pos(1,1) instr_pos(9,1)],'km');
p=10:10:4200;

clear dist
for i=1:9
    dist(i)=1000*sw_dist([instr_pos(1,2) instr_pos(i,2)],[instr_pos(1,1) instr_pos(i,1)],'km');
end

for i=1:length(p)
    xgrid(i,:)=x(:);
end
for i=1:length(x)
    pgrid(:,i)=p(:);
end

xc_l=130*1000; %horizontal correlation length from cross correlations (m)
zc_l=790; %vertical correlation length from cross correlations (m)

x_corr_func=@(x) exp(-(x(:)/xc_l).^2).*cos(pi.*x(:)./(2.*xc_l)); %horiz correlation function
z_corr_func=@(z) exp(-(z(:)/zc_l).^2); %vertical correlation function


% set up the array of dx of observations
total=sum(~isnan(fake_grid.temp(:,:)));
clear dx_obs
for i=1:sum(total)
    for j=1:length(total)
        if sum(total(1:j))==i
            dx_obs(i)=dist(j);
        end
    end
end

for i=61:-1:1
    %if sum(dx_obs==0)>0
        if dx_obs(i)==0
            dx_obs(i)=dx_obs(i+1);
        end
    %end
end

rng=s;
a = 0.01;
b = .1;
r = (b-a).*rand(1000,1) + a;

clear ratio
ratio=zeros(length(dx_obs),length(dx_obs)); % just do the ratio since 
for i=1:length(dx_obs)
    for j=1:length(dx_obs)
        if i==j
            ratio(i,j)=10*r(i);
        end
    end
end

% need a dz of obs
clear dp_obs
dp_obs=[fake_grid.pres(~isnan(fake_grid.temp(:,1)),1);fake_grid.pres(~isnan(fake_grid.temp(:,2)),2);fake_grid.pres(~isnan(fake_grid.temp(:,3)),3);fake_grid.pres(~isnan(fake_grid.temp(:,4)),4);fake_grid.pres(~isnan(fake_grid.temp(:,5)),5);fake_grid.pres(~isnan(fake_grid.temp(:,6)),6);fake_grid.pres(~isnan(fake_grid.temp(:,7)),7);fake_grid.pres(~isnan(fake_grid.temp(:,8)),8);fake_grid.pres(~isnan(fake_grid.temp(:,9)),9)];

% need array of the observations without any NaN
temp_obs=[fake_grid.theta(~isnan(fake_grid.temp(:,1)),1);fake_grid.theta(~isnan(fake_grid.temp(:,2)),2);fake_grid.theta(~isnan(fake_grid.temp(:,3)),3);fake_grid.theta(~isnan(fake_grid.temp(:,4)),4);fake_grid.theta(~isnan(fake_grid.temp(:,5)),5);fake_grid.theta(~isnan(fake_grid.temp(:,6)),6);fake_grid.theta(~isnan(fake_grid.temp(:,7)),7);fake_grid.theta(~isnan(fake_grid.temp(:,8)),8);fake_grid.theta(~isnan(fake_grid.temp(:,9)),9)];

% would loop over time for microcat data but we only have 1 snapshot
clear weight_corr
for i=1:length(ratio)
    for j=1:size(pgrid,1)
        for k=1:size(pgrid,2)
            weight_corr(i,j,k)=x_corr_func(abs(xgrid(j,k)-dx_obs(1,i)))*z_corr_func(abs(pgrid(j,k)-dp_obs(i,1)));
        end
    end
end

clear cross_corr
for i=1:length(ratio)
    for j=1:length(ratio)
        cross_corr(i,j)=x_corr_func(abs(dx_obs(i)-dx_obs(j)))*z_corr_func(abs(dp_obs(i)-dp_obs(j)));
    end
end

clear weights
for j=1:size(pgrid,1)
    for k=1:size(pgrid,2)
        weights(:,j,k)=(ratio+cross_corr)\weight_corr(:,j,k);
    end
end

for j=1:size(pgrid,1)
    for k=1:size(pgrid,2)
        int_t(j,k)=weights(:,j,k).'*temp_obs(:); %check to see if dimensions work?
    end
end

%now subtract and do smaller scale
for i=1:9
    dist_coast(i)=1000*sw_dist([c2_pos(2) fake_grid.lat(i)],[c2_pos(1) fake_grid.lon(i)],'km');
end

for i=1:9
    for j=1:size(xgrid,2)
        grid_dist(i,j)=abs(xgrid(1,j)-dist_coast(i));
    end
end

for i=1:9
    [M(i),I(i)] = min(grid_dist(i,:));
end

for i=1:9
    for j=1:size(xgrid,1)
        fake_grid.t_anom(j,i)=fake_grid.theta(j,i)-int_t(j,I(i));
    end
end

xc=130*1000*small_scale_ratio; %horizontal correlation length from cross correlations (m)
zc=790*small_scale_ratio; %vertical correlation length from cross correlations (m)
ratio=ratio*small_scale_ratio;

t_anom_obs=[fake_grid.t_anom(~isnan(fake_grid.temp(:,1)),1);fake_grid.t_anom(~isnan(fake_grid.temp(:,2)),2);fake_grid.t_anom(~isnan(fake_grid.temp(:,3)),3);fake_grid.t_anom(~isnan(fake_grid.temp(:,4)),4);fake_grid.t_anom(~isnan(fake_grid.temp(:,5)),5);fake_grid.t_anom(~isnan(fake_grid.temp(:,6)),6);fake_grid.t_anom(~isnan(fake_grid.temp(:,7)),7);fake_grid.t_anom(~isnan(fake_grid.temp(:,8)),8);fake_grid.t_anom(~isnan(fake_grid.temp(:,9)),9)];

% would loop over time for microcat data but we only have 1 snapshot
for i=1:length(ratio)
    for j=1:size(pgrid,1)
        for k=1:size(pgrid,2)
            weight_corr(i,j,k)=x_corr_func(abs(xgrid(j,k)-dx_obs(1,i)))*z_corr_func(abs(pgrid(j,k)-dp_obs(i,1)));
        end
    end
end

for i=1:length(ratio)
    for j=1:length(ratio)
        cross_corr(i,j)=x_corr_func(abs(dx_obs(i)-dx_obs(j)))*z_corr_func(abs(dp_obs(i)-dp_obs(j)));
    end
end

for j=1:size(pgrid,1)
    for k=1:size(pgrid,2)
        weights(:,j,k)=(ratio+cross_corr)\weight_corr(:,j,k);
    end
end

for j=1:size(pgrid,1)
    for k=1:size(pgrid,2)
        int_t2(j,k)=weights(:,j,k).'*t_anom_obs(:); %check to see if dimensions work?
    end
end

background_error=nan(size(xgrid,1),size(xgrid,2));
for i=1:size(pgrid,1)
    for j=1:size(pgrid,2)
        background_error(i,j)=1; %let's say it's .01 C everywhere, see what happens
    end
end

clear analysis_error_t
for i=1:size(pgrid,1)
    for j=1:size(pgrid,2)
        analysis_error_t(i,j)=background_error(i,j)-(weight_corr(:,i,j).'*inv(ratio+cross_corr)*weight_corr(:,i,j));
    end
end

% compare large, small scale fields, reality
clear n
for i=300:313
    n(i-299)=1000*sw_dist([ctdg.lat(i) instr_pos(1,2)],[ctdg.lon(i) instr_pos(1,1)],'km');
end

for i=1:501
    for j=1:2415
        ctdg.theta(i,j)=gsw_pt_from_t(ctdg.sal(i,j),ctdg.temp(i,j),ctdg.pres(i,j));
    end
end

figure
subplot(2,2,1)
hold on
ax = gca;
ax.FontSize = 24;
[C,h]=contourf(n/1000,ctdg.pgrid(:,1),ctdg.theta(:,300:313),0:1:30)
clabel(C,h,0:5:30,'FontSize',20)
axis 'ij'
%xlabel('Distance (km)','FontSize',24)
ylabel('Pressure','FontSize',24)
%colorbar
cmocean('thermal')
title('2013 hydrography subsampled to ASCA instrument locations (temperature)','FontSize',24)
brighten(.3) %how to shift whole colormap to lower values?
for i=1:9
    plot([dist(i)/1000 dist(i)/1000],[0 4200],'k')
end
for i=1:420
    for j=1:9
        if ~isnan(fake_grid.temp(i,j))
            scatter(dist(j)/1000,fake_grid.depth(i,j),'k','filled')
        end
    end
end
caxis([0 25])
ylim([0 500])
xlim([0 180])
hold off

subplot(2,2,3)
hold on
ax = gca;
ax.FontSize = 24;
[C,h]=contourf(n/1000,ctdg.pgrid(:,1),ctdg.theta(:,300:313),0:1:30)
clabel(C,h,0:5:30,'FontSize',20)
axis 'ij'
xlabel('Distance (km)','FontSize',24)
ylabel('Pressure','FontSize',24)
%colorbar
cmocean('thermal')
%title('2013 hydrography subsampled to ASCA instrument locations (temperature)','FontSize',24)
brighten(.3) %how to shift whole colormap to lower values?
for i=1:9
    plot([dist(i)/1000 dist(i)/1000],[0 4200],'k')
end
for i=1:420
    for j=1:9
        if ~isnan(fake_grid.temp(i,j))
            scatter(dist(j)/1000,fake_grid.depth(i,j),'k','filled')
        end
    end
end
caxis([0 25])
ylim([500 4200])
xlim([0 180])
hold off

% subplot(2,2,2)
% hold on
% [C,h]=contourf(xgrid,pgrid,int_t,0:2.5:30)
% axis 'ij'
% cmocean('thermal')
% xlabel('Distance (m)','FontSize',24)
% ylabel('Pressure','FontSize',24)
% title('Large scale','FontSize',24)
% colorbar
% clabel(C,h,0:5:30)
% brighten(.3)
% caxis([0 25])
% xlim([0 180*1000])
% hold off
% 
% 
% subplot(2,2,3)
% hold on
% [C,h]=contourf(xgrid,pgrid,int_t2,-10:.25:10)
% axis 'ij'
% xlabel('Distance (m)','FontSize',24)
% ylabel('Pressure','FontSize',24)
% title('Small scale','FontSize',24)
% colorbar
% clabel(C,h,-10:.5:10)
% cmocean('balance','zero')
% %caxis([0 25])
% brighten(.3)
% xlim([0 180*1000])
% hold off

subplot(2,2,2)
hold on
ax = gca;
ax.FontSize = 24;
[C,h]=contourf(xgrid/1000,pgrid,int_t+int_t2,0:1:30)
clabel(C,h,0:5:30,'FontSize',20)
axis 'ij'
cmocean('thermal')
%xlabel('Distance (km)','FontSize',24)
%ylabel('Pressure','FontSize',24)
title('Optimally interpolated temperature','FontSize',24)
colorbar
caxis([0 25])
xlim([0 180])
ylim([0 500])
brighten(.3)

subplot(2,2,4)
hold on
ax = gca;
ax.FontSize = 24;
[C,h]=contourf(xgrid/1000,pgrid,int_t+int_t2,0:1:30)
clabel(C,h,0:5:30,'FontSize',20)
axis 'ij'
cmocean('thermal')
xlabel('Distance (km)','FontSize',24)
%ylabel('Pressure','FontSize',24)
%title('Optimally interpolated temperature','FontSize',24)
colorbar
caxis([0 25])
xlim([0 180])
ylim([500 4200])
brighten(.3)

% map of error given from small scale OI
% figure
% hold on
% [C,h]=contourf(xgrid,pgrid,analysis_error,-10:.01:10)
% axis 'ij'
% cmocean('curl','zero')
% xlabel('Distance (m)','FontSize',24)
% ylabel('Pressure','FontSize',24)
% title('Expected error ratio: temperature','FontSize',24)
% colorbar
% clabel(C,h,-10:.05:10)
% caxis([0 .3])
% xlim([0 280*1000])
% brighten(.3)


%% salinity
xc_l=77*1000*(1/small_scale_ratio); %horizontal correlation length from cross correlations (m)
zc_l=171*(1/small_scale_ratio); %vertical correlation length from cross correlations (m)

x_corr_func=@(x) exp(-(x(:)/xc_l).^2).*cos(pi.*x(:)./(2.*xc_l)); %horiz correlation function
z_corr_func=@(z) exp(-(z(:)/zc_l).^2); %vertical correlation function

% set up the array of dx of observations
clear dx_obs
total=sum(~isnan(fake_grid.temp(:,:)));
for i=1:sum(total)
    for j=1:length(total)
        if sum(total(1:j))==i
            dx_obs(i)=dist(j);
        end
    end
end

for i=61:-1:1
    %if sum(dx_obs==0)>0
        if dx_obs(i)==0
            dx_obs(i)=dx_obs(i+1);
        end
    %end
end

rng=t;
a = 0.01;
b = .1;
r = (b-a).*rand(1000,1) + a;

ratio=zeros(length(dx_obs),length(dx_obs)); % just do the ratio since 
for i=1:length(dx_obs)
    for j=1:length(dx_obs)
        if i==j
            ratio(i,j)=10*r(i);
        end
    end
end

% need a dz of obs
dp_obs=[fake_grid.pres(~isnan(fake_grid.temp(:,1)),1);fake_grid.pres(~isnan(fake_grid.temp(:,2)),2);fake_grid.pres(~isnan(fake_grid.temp(:,3)),3);fake_grid.pres(~isnan(fake_grid.temp(:,4)),4);fake_grid.pres(~isnan(fake_grid.temp(:,5)),5);fake_grid.pres(~isnan(fake_grid.temp(:,6)),6);fake_grid.pres(~isnan(fake_grid.temp(:,7)),7);fake_grid.pres(~isnan(fake_grid.temp(:,8)),8);fake_grid.pres(~isnan(fake_grid.temp(:,9)),9)];

% need array of the observations without any NaN
sal_obs=[fake_grid.sal(~isnan(fake_grid.sal(:,1)),1);fake_grid.sal(~isnan(fake_grid.sal(:,2)),2);fake_grid.sal(~isnan(fake_grid.sal(:,3)),3);fake_grid.sal(~isnan(fake_grid.sal(:,4)),4);fake_grid.sal(~isnan(fake_grid.sal(:,5)),5);fake_grid.sal(~isnan(fake_grid.sal(:,6)),6);fake_grid.sal(~isnan(fake_grid.sal(:,7)),7);fake_grid.sal(~isnan(fake_grid.sal(:,8)),8);fake_grid.sal(~isnan(fake_grid.sal(:,9)),9)];

sal_obs=sal_obs-35;

% would loop over time for microcat data but we only have 1 snapshot
clear weight_corr
for i=1:length(ratio)
    for j=1:size(pgrid,1)
        for k=1:size(pgrid,2)
            weight_corr(i,j,k)=x_corr_func(abs(xgrid(j,k)-dx_obs(1,i)))*z_corr_func(abs(pgrid(j,k)-dp_obs(i,1)));
        end
    end
end

clear cross_corr
for i=1:length(ratio)
    for j=1:length(ratio)
        cross_corr(i,j)=x_corr_func(abs(dx_obs(i)-dx_obs(j)))*z_corr_func(abs(dp_obs(i)-dp_obs(j)));
    end
end

clear weights
for j=1:size(pgrid,1)
    for k=1:size(pgrid,2)
        weights(:,j,k)=(ratio+cross_corr)\weight_corr(:,j,k);
    end
end

for j=1:size(pgrid,1)
    for k=1:size(pgrid,2)
        int_s(j,k)=weights(:,j,k).'*sal_obs(:); %check to see if dimensions work?
    end
end

% now subtract and do smaller scale
for i=1:9
    dist_coast(i)=1000*sw_dist([c2_pos(2) fake_grid.lat(i)],[c2_pos(1) fake_grid.lon(i)],'km');
end

for i=1:9
    for j=1:347
        grid_dist(i,j)=abs(xgrid(1,j)-dist_coast(i));
    end
end

for i=1:9
    [M(i),I(i)] = min(grid_dist(i,:));
end

for i=1:9
    for j=1:420
        fake_grid.s_anom(j,i)=fake_grid.sal(j,i)-35-int_s(j,I(i));
    end
end

xc=77*1000;%*small_scale_ratio; %horizontal correlation length from cross correlations (m) * .1
zc=171; %vertical correlation length from cross correlations (m) * .1
ratio=ratio*small_scale_ratio;

s_anom_obs=[fake_grid.s_anom(~isnan(fake_grid.sal(:,1)),1);fake_grid.s_anom(~isnan(fake_grid.sal(:,2)),2);fake_grid.s_anom(~isnan(fake_grid.sal(:,3)),3);fake_grid.s_anom(~isnan(fake_grid.sal(:,4)),4);fake_grid.s_anom(~isnan(fake_grid.sal(:,5)),5);fake_grid.s_anom(~isnan(fake_grid.sal(:,6)),6);fake_grid.s_anom(~isnan(fake_grid.sal(:,7)),7);fake_grid.s_anom(~isnan(fake_grid.sal(:,8)),8);fake_grid.s_anom(~isnan(fake_grid.sal(:,9)),9)];

% would loop over time for microcat data but we only have 1 snapshot
for i=1:length(ratio)
    for j=1:size(pgrid,1)
        for k=1:size(pgrid,2)
            weight_corr(i,j,k)=x_corr_func(abs(xgrid(j,k)-dx_obs(1,i)))*z_corr_func(abs(pgrid(j,k)-dp_obs(i,1)));
        end
    end
end

for i=1:length(ratio)
    for j=1:length(ratio)
        cross_corr(i,j)=x_corr_func(abs(dx_obs(i)-dx_obs(j)))*z_corr_func(abs(dp_obs(i)-dp_obs(j)));
    end
end

for j=1:size(pgrid,1)
    for k=1:size(pgrid,2)
        weights(:,j,k)=(ratio+cross_corr)\weight_corr(:,j,k);
    end
end

for j=1:size(pgrid,1)
    for k=1:size(pgrid,2)
        int_s2(j,k)=weights(:,j,k).'*s_anom_obs(:); %check to see if dimensions work?
    end
end

% % let's do error from this step
background_error=nan(420,347);
for i=1:size(pgrid,1)
    for j=1:size(pgrid,2)
        background_error(i,j)=1; 
    end
end

clear analysis_error_s
for i=1:size(pgrid,1)
    for j=1:size(pgrid,2)
        analysis_error_s(i,j)=background_error(i,j)-(weight_corr(:,i,j).'*inv(ratio+cross_corr)*weight_corr(:,i,j));
    end
end

%compare large, small scale fields, reality
clear n
for i=300:313
    n(i-299)=1000*sw_dist([ctdg.lat(i) instr_pos(1,2)],[ctdg.lon(i) instr_pos(1,1)],'km');
end

figure
subplot(2,2,1)
hold on
ax = gca;
ax.FontSize = 24;
[C,h]=contourf(n/1000,ctdg.pgrid(:,1),ctdg.sal(:,300:313),30:.1:40)
clabel(C,h,30:.2:40,'FontSize',20)
axis 'ij'
%xlabel('Distance (km)','FontSize',24)
ylabel('Pressure','FontSize',24)
%colorbar
title('2013 hydrography subsampled to ASCA instrument locations (salinity)','FontSize',24)
for i=1:9
    plot([dist(i)/1000 dist(i)/1000],[0 4200],'k')
end
for i=1:420
    for j=1:9
        if ~isnan(fake_grid.temp(i,j))
            scatter(dist(j)/1000,fake_grid.depth(i,j),'k','filled')
        end
    end
end
caxis([34.3 35.7])
ylim([0 500])
xlim([0 180])
cmocean('haline')
brighten(.3)
hold off

subplot(2,2,3)
hold on
ax = gca;
ax.FontSize = 24;
[C,h]=contourf(n/1000,ctdg.pgrid(:,1),ctdg.sal(:,300:313),30:.1:40)
clabel(C,h,30:.2:40,'FontSize',20)
axis 'ij'
xlabel('Distance (km)','FontSize',24)
ylabel('Pressure','FontSize',24)
%colorbar
%title('2013 hydrography subsampled to ASCA instrument locations (salinity)','FontSize',24)
for i=1:9
    plot([dist(i)/1000 dist(i)/1000],[0 4200],'k')
end
for i=1:420
    for j=1:9
        if ~isnan(fake_grid.temp(i,j))
            scatter(dist(j)/1000,fake_grid.depth(i,j),'k','filled')
        end
    end
end
caxis([34.3 35.7])
ylim([500 4200])
xlim([0 180])
cmocean('haline')
brighten(.3) %how to shift whole colormap to lower values?
hold off

% figure
% hold on
% [C,h]=contourf(xgrid,pgrid,int_s,0:.25:40)
% axis 'ij'
% cmocean('haline')
% xlabel('Distance (m)','FontSize',24)
% ylabel('Pressure','FontSize',24)
% title('Large scale','FontSize',24)
% colorbar
% clabel(C,h,0:.5:40)
% brighten(.3)
% caxis([34 36])
% xlim([0 180*1000])
% hold off


% figure
% hold on
% [C,h]=contourf(xgrid,pgrid,int_s2,-10:.25:10)
% axis 'ij'
% cmocean('balance','zero')
% xlabel('Distance (m)','FontSize',24)
% ylabel('Pressure','FontSize',24)
% title('Small scale','FontSize',24)
% colorbar
% clabel(C,h,-10:.5:10)
% caxis([0 25])
% brighten(.3)
% xlim([0 180*1000])
% hold off

subplot(2,2,2)
hold on
ax = gca;
ax.FontSize = 24;
[C,h]=contourf(xgrid/1000,pgrid,int_s2+int_s+35,30:.1:40)
axis 'ij'
%xlabel('Distance (km)','FontSize',24)
ylabel('Pressure','FontSize',24)
title('Optimally interpolated salinity','FontSize',24)
colorbar
clabel(C,h,30:.2:40,'FontSize',20)
caxis([34.3 35.7])
xlim([0 180])
ylim([0 500])
cmocean('haline')
brighten(.3)
hold off

subplot(2,2,4)
hold on
ax = gca;
ax.FontSize = 24;
[C,h]=contourf(xgrid/1000,pgrid,int_s2+int_s+35,30:.1:40)
axis 'ij'
xlabel('Distance (km)','FontSize',24)
ylabel('Pressure','FontSize',24)
%title('Optimally interpolated salinity','FontSize',24)
colorbar
clabel(C,h,30:.2:40,'FontSize',20)
caxis([34.3 35.7])
xlim([0 180])
ylim([500 4200])
cmocean('haline')
brighten(.3)

% 
% % figure
% % hold on
% % [C,h]=contourf(xgrid,zgrid,int_s3+int_s2+int_s,30:.25:40)
% % axis 'ij'
% % cmocean('haline')
% % xlabel('Distance (m)','FontSize',24)
% % ylabel('Pressure','FontSize',24)
% % title('Large scale + small scale x2','FontSize',24)
% % colorbar
% % clabel(C,h,30:.5:40)
% % caxis([34 36])
% % xlim([0 180*1000])
% % brighten(.3)
% 
% 
% % map of error given from small scale OI
% figure
% hold on
% [C,h]=contourf(xgrid,zgrid,analysis_error,-10:.005:10)
% axis 'ij'
% cmocean('curl','zero')
% xlabel('Distance (m)','FontSize',24)
% ylabel('Pressure','FontSize',24)
% title('Expected error ratio: salinity','FontSize',24)
% colorbar
% clabel(C,h,-10:.05:10)
% caxis([0 .15])
% xlim([0 180*1000])
% brighten(.3)
%% u and v

clear dist
for i=1:9
    dist(i)=1000*sw_dist([instr_pos(1,2) instr_pos(i,2)],[instr_pos(1,1) instr_pos(i,1)],'km');
end

fake_grid.vel=nan(420,9,2);
fake_grid.vel(:,:,1)=fake_grid.u(:,:);
fake_grid.vel(:,:,2)=fake_grid.v(:,:);

% to do super obs
for i=6:10:46
    for l=1:2
        for j=4:9
            fake_grid.vel(i,j,l)=.1*fake_grid.vel(i-4,j,l)+.2*fake_grid.vel(i-2,j,l)+.4*fake_grid.vel(i,j,l)+.2*fake_grid.vel(i+2,j,l)+.1*fake_grid.vel(i+4,j,l);
        end
    end
end

for i=[4,8,10,12,14,18,20,22,24,28,30,32,34,38,40,42,44,48]
    for l=1:2
        for j=4:9
            fake_grid.vel(i,j,l)=NaN;
        end
    end
end

for l=1:2
    for j=1:2
        fake_grid.vel(6,j,l)=.2*fake_grid.vel(4,j,l)+.6*fake_grid.vel(6,j,l)+.2*fake_grid.vel(8,j,l);
    end
    fake_grid.vel(6,3,l)=.1*fake_grid.vel(2,j,l)+.2*fake_grid.vel(4,j,l)+.4*fake_grid.vel(6,j,l)+.2*fake_grid.vel(8,j,l)+.1*fake_grid.vel(10,j,l);
    fake_grid.vel(16,3,l)=.1*fake_grid.vel(12,j,l)+.2*fake_grid.vel(14,j,l)+.4*fake_grid.vel(16,j,l)+.2*fake_grid.vel(18,j,l)+.1*fake_grid.vel(20,j,l);
    fake_grid.vel(25,3,l)=(1/6)*fake_grid.vel(22,3,l)+(1/3)*fake_grid.vel(24,3,l)+(1/3)*fake_grid.vel(26,3,l)+(1/6)*fake_grid.vel(28,3,l);
end

for l=1:2
    for j=1:2
        fake_grid.vel(4,j,l)=NaN;
        fake_grid.vel(8,j,l)=NaN;
    end
    for j=[2,4,8,10,12,14,18,20,22,24,26,28]
        fake_grid.vel(j,3,l)=NaN;
    end
end

% do simultaneously as talked about by mariano

xc_l=55*1000; %horizontal correlation length from cross correlations (m)
zc_l=2200; %vertical correlation length from cross correlations (m)

x_corr_func=@(x) exp(-(x(:)/xc_l).^2).*cos(pi.*x(:)./(2.*xc_l)); %horiz correlation function
z_corr_func=@(z) exp(-(z(:)/zc_l).^2); %vertical correlation function

% set up the array of dx of observations
clear dx_obs
total=sum(~isnan(fake_grid.vel(:,:,1)));
for i=1:sum(total)
    for j=1:length(total)
        if sum(total(1:j))==i
            dx_obs(i)=dist(j);
        end
    end
end

for i=length(dx_obs):-1:1
    %if sum(dx_obs==0)>0
        if dx_obs(i)==0
            dx_obs(i)=dx_obs(i+1);
        end
    %end
end

rng=s;
a = 0.01;
b = .1;
r = (b-a).*rand(1000,1) + a;

clear ratio
ratio=zeros(length(dx_obs),length(dx_obs)); % just do the ratio since 
for i=1:length(dx_obs)
    for j=1:length(dx_obs)
        if i==j
            ratio(i,j)=10*r(i);
        end
    end
end

clear x
clear xgrid
clear z
clear zgrid
x=0:500:1000*sw_dist([instr_pos(1,2) instr_pos(9,2)],[instr_pos(1,1) instr_pos(9,1)],'km');
z=10:10:4200; % then later will make anything under topography into NaN
for i=1:length(z)
    xgrid(i,:)=x(:);
end
for i=1:length(x)
    zgrid(:,i)=z(:);
end

% need array of the observations without any NaN
clear vel_obs
for i=1:2
    vel_obs(:,i)=[fake_grid.vel(~isnan(fake_grid.vel(:,1)),1,i);fake_grid.vel(~isnan(fake_grid.vel(:,2)),2,i);fake_grid.vel(~isnan(fake_grid.vel(:,3)),3,i);fake_grid.vel(~isnan(fake_grid.vel(:,4)),4,i);fake_grid.vel(~isnan(fake_grid.vel(:,5)),5,i);fake_grid.vel(~isnan(fake_grid.vel(:,6)),6,i);fake_grid.vel(~isnan(fake_grid.vel(:,7)),7,i);fake_grid.vel(~isnan(fake_grid.vel(:,8)),8,i);fake_grid.vel(~isnan(fake_grid.vel(:,9)),9,i)];
end

% dz of this
clear dz_obs
for i=1
    dz_obs=[fake_grid.depth(~isnan(fake_grid.vel(:,1)),1,i);fake_grid.depth(~isnan(fake_grid.vel(:,2)),2,i);fake_grid.depth(~isnan(fake_grid.vel(:,3)),3,i);fake_grid.depth(~isnan(fake_grid.vel(:,4)),4,i);fake_grid.depth(~isnan(fake_grid.vel(:,5)),5,i);fake_grid.depth(~isnan(fake_grid.vel(:,6)),6,i);fake_grid.depth(~isnan(fake_grid.vel(:,7)),7,i);fake_grid.depth(~isnan(fake_grid.vel(:,8)),8,i);fake_grid.depth(~isnan(fake_grid.vel(:,9)),9,i)];
end

% would loop over time for microcat data but we only have 1 snapshot
clear weight_corr
for i=1:length(ratio)
    for j=1:size(zgrid,1)
        for k=1:size(zgrid,2)
            for l=1:2
                weight_corr(i,j,k,l)=x_corr_func(abs(xgrid(j,k)-dx_obs(1,i)))*z_corr_func(abs(zgrid(j,k)-dz_obs(i,1)));
            end %weight_corr would depend on l (U or V measurement) if we make the decorrelation length or function anisotropic
        end
    end
end

clear cross_corr
for i=1:length(ratio)
    for j=1:length(ratio)
        for l=1:2
            cross_corr(i,j,l)=x_corr_func(abs(dx_obs(i)-dx_obs(j)))*z_corr_func(abs(dz_obs(i)-dz_obs(j)));
        end %similar - could be dependent on l if we want
    end
end

clear weights
for j=1:size(zgrid,1)
    for k=1:size(zgrid,2)
        for l=1:2
            weights(:,j,k,l)=(ratio(:,:)+cross_corr(:,:,l))\weight_corr(:,j,k,l);
        end % again - could depend on l
    end
end

clear int_vel
for j=1:size(zgrid,1)
    for k=1:size(zgrid,2)
        for l=1:2
            int_vel(j,k,l)=weights(:,j,k,l).'*vel_obs(:,l); %check to see if dimensions work?
        end
    end
end

% now subtract and do smaller scale
for i=1:9
    dist_coast(i)=1000*sw_dist([c2_pos(2) fake_grid.lat(i)],[c2_pos(1) fake_grid.lon(i)],'km');
end

for i=1:9
    for j=1:347
        grid_dist(i,j)=abs(xgrid(1,j)-dist_coast(i));
    end
end

clear M
clear I
for i=1:9
    [M(i),I(i)] = min(grid_dist(i,:));
end

for i=1:9
    for j=1:420
        for l=1:2
            fake_grid.vel_anom(j,i,l)=fake_grid.vel(j,i,l)-int_vel(j,I(i),l);
        end
    end
end

clear vel_anom_obs
for i=1:2
    vel_anom_obs(:,i)=[fake_grid.vel_anom(~isnan(fake_grid.vel(:,1)),1,i);fake_grid.vel_anom(~isnan(fake_grid.vel(:,2)),2,i);fake_grid.vel_anom(~isnan(fake_grid.vel(:,3)),3,i);fake_grid.vel_anom(~isnan(fake_grid.vel(:,4)),4,i);fake_grid.vel_anom(~isnan(fake_grid.vel(:,5)),5,i);fake_grid.vel_anom(~isnan(fake_grid.vel(:,6)),6,i);fake_grid.vel_anom(~isnan(fake_grid.vel(:,7)),7,i);fake_grid.vel_anom(~isnan(fake_grid.vel(:,8)),8,i);fake_grid.vel_anom(~isnan(fake_grid.vel(:,9)),9,i)];
end

xc=55*1000*small_scale_ratio; %horizontal correlation length from cross correlations (m)
zc=2200*small_scale_ratio; %vertical correlation length from cross correlations (m)

x_corr_func=@(x) exp(-(x(:)/xc).^2).*cos(pi.*x(:)./(2.*xc)); %horiz correlation function
z_corr_func=@(z) exp(-(z(:)/zc).^2); %vertical correlation function
ratio=ratio*small_scale_ratio;

% would loop over time for microcat data but we only have 1 snapshot
for i=1:length(ratio)
    for j=1:size(zgrid,1)
        for k=1:size(zgrid,2)
            for l=1:2
                weight_corr(i,j,k,l)=x_corr_func(abs(xgrid(j,k)-dx_obs(1,i)))*z_corr_func(abs(zgrid(j,k)-dz_obs(i,1)));
            end %weight_corr would depend on l (U or V measurement) if we make the decorrelation length or function anisotropic
        end
    end
end

for i=1:length(ratio)
    for j=1:length(ratio)
        for l=1:2
            cross_corr(i,j,l)=x_corr_func(abs(dx_obs(i)-dx_obs(j)))*z_corr_func(abs(dz_obs(i)-dz_obs(j)));
        end %similar - could be dependent on l if we want
    end
end

clear weights
for j=1:size(zgrid,1)
    for k=1:size(zgrid,2)
        for l=1:2
            weights(:,j,k,l)=(ratio(:,:)+cross_corr(:,:,l))\weight_corr(:,j,k,l);
        end % again - could depend on l
    end
end

for j=1:size(zgrid,1)
    for k=1:size(zgrid,2)
        for l=1:2
            int_vel2(j,k,l)=weights(:,j,k,l).'*vel_anom_obs(:,l); %check to see if dimensions work?
        end
    end
end

% compare large, small scale fields, reality
clear n
for i=1:19
    n(i)=1000*sw_dist([ldp.lat(i) instr_pos(1,2)],[ldp.lon(i) instr_pos(1,1)],'km');
end

figure
subplot(2,2,1)
hold on
ax = gca;
ax.FontSize = 24;
[C,h]=contourf(n/1000,20:20:4620,ldp.v(:,1:19),-2:.1:3)
clabel(C,h,-2:.2:3,'FontSize',20)
axis 'ij'
%xlabel('Distance (km)','FontSize',24)
ylabel('Depth','FontSize',24)
%colorbar
caxis([-2 .2])
cmocean('balance','zero')
title('2013 LADCP subsampled to ASCA instrument locations (v)','FontSize',24)
brighten(.3) %how to shift whole colormap to lower values?
for i=1:9
    plot([dist(i)/1000 dist(i)/1000],[0 4200],'k')
end
for i=1:420
    for j=1:9
        if ~isnan(fake_grid.u(i,j))
            scatter(dist(j)/1000,fake_grid.depth(i,j),'k','filled')
        end
    end
end
ylim([0 500])
xlim([0 180])
hold off

subplot(2,2,3)
hold on
ax = gca;
ax.FontSize = 24;
[C,h]=contourf(n/1000,20:20:4620,ldp.v(:,1:19),-2:.1:3)
clabel(C,h,-2:.2:3,'FontSize',20)
axis 'ij'
xlabel('Distance (km)','FontSize',24)
ylabel('Depth','FontSize',24)
%colorbar
caxis([-2 .2])
cmocean('balance','zero')
%title('2013 LADCP subsampled to ASCA instrument locations (u)','FontSize',24)
brighten(.3) %how to shift whole colormap to lower values?
for i=1:9
    plot([dist(i)/1000 dist(i)/1000],[0 4200],'k')
end
for i=1:420
    for j=1:9
        if ~isnan(fake_grid.u(i,j))
            scatter(dist(j)/1000,fake_grid.depth(i,j),'k','filled')
        end
    end
end
ylim([500 4200])
xlim([0 180])
hold off
% 
% figure
% hold on
% [C,h]=contourf(xgrid,zgrid,int_vel(:,:,2),-2:.1:3)
% axis 'ij'
% caxis([-2 .2])
% cmocean('balance','zero')
% xlabel('Distance (m)','FontSize',24)
% ylabel('Depth','FontSize',24)
% title('Large scale v','FontSize',24)
% colorbar
% clabel(C,h,-2:.5:3)
% brighten(.3)
% xlim([0 180*1000])
% hold off

% 
% figure
% hold on
% [C,h]=contourf(xgrid,zgrid,int_vel2(:,:,2),-3:.1:3)
% axis 'ij'
% cmocean('balance','zero')
% xlabel('Distance (m)','FontSize',24)
% ylabel('Pressure','FontSize',24)
% title('Small scale','FontSize',24)
% colorbar
% clabel(C,h,-2:.2:3)
% %caxis([0 25])
% brighten(.3)
% xlim([0 180*1000])
% hold off
subplot(2,2,2)
hold on
ax = gca;
ax.FontSize = 24;
[C,h]=contourf(xgrid/1000,zgrid,int_vel2(:,:,2)+int_vel(:,:,2),-3:.1:2)
axis 'ij'
caxis([-2 .2])
cmocean('balance','zero')
%xlabel('Distance (km)','FontSize',24)
ylabel('Depth','FontSize',24)
title('Optimally interpolated v velocity','FontSize',24)
colorbar
clabel(C,h,-3:.2:2,'FontSize',20)
xlim([0 180])
ylim([0 500])
brighten(.3)

subplot(2,2,4)
hold on
ax = gca;
ax.FontSize = 24;
[C,h]=contourf(xgrid/1000,zgrid,int_vel2(:,:,2)+int_vel(:,:,2),-3:.1:2)
axis 'ij'
caxis([-2 .2])
cmocean('balance','zero')
xlabel('Distance (km)','FontSize',24)
ylabel('Depth','FontSize',24)
%title('Optimally interpolated v velocity','FontSize',24)
colorbar
clabel(C,h,-3:.2:2,'FontSize',20)
xlim([0 180])
ylim([500 4200])
brighten(.3)
% 
% % figure
% % hold on
% % [C,h]=contourf(xgrid,zgrid,int_vel3(:,:,1)+int_vel2(:,:,1)+int_vel(:,:,1),-3:.1:2)
% % axis 'ij'
% % caxis([-2 .2])
% % cmocean('balance','zero')
% % xlabel('Distance (m)','FontSize',24)
% % ylabel('Pressure','FontSize',24)
% % title('Large scale + small scale x2','FontSize',24)
% % colorbar
% % clabel(C,h,-3:.5:2)
% % xlim([0 180*1000])
% % brighten(.3)
% 
% % map of error given from small scale OI
% figure
% hold on  
% [C,h]=contourf(xgrid,zgrid,analysis_error(:,:,1),-10:.1:10)
% axis 'ij'
% cmocean('curl','zero')
% xlabel('Distance (m)','FontSize',24)
% ylabel('Pressure','FontSize',24)
% title('Expected error ratio: large scale u','FontSize',24)
% colorbar
% clabel(C,h,-10:.1:10)
% caxis([0 1])
% xlim([0 180*1000])
% brighten(.3)
% 
% figure
% hold on  
% [C,h]=contourf(xgrid,zgrid,analysis_error2(:,:,1),-10:.1:10)
% axis 'ij'
% cmocean('curl','zero')
% xlabel('Distance (m)','FontSize',24)
% ylabel('Pressure','FontSize',24)
% title('Expected error ratio: small scale u','FontSize',24)
% colorbar
% clabel(C,h,-10:.1:10)
% caxis([0 1])
% xlim([0 180*1000])
% brighten(.3)

%% get rid of stuff under topography

for i=1:347
    for j=1:420
        if p(j)>abs(topo(i))
            int_s(j,i)=NaN;
            int_s2(j,i)=NaN;
            int_t(j,i)=NaN;
            int_t2(j,i)=NaN;
        end
        for k=1:2
            if z(j)>abs(topo(i))
                int_vel(j,i,k)=NaN;
                int_vel2(j,i,k)=NaN;
            end
        end
    end
end

% plot up analysis errors to look at instrument placement

figure
hold on
[C,h]=contourf(xgrid,pgrid,analysis_error_t)
axis 'ij'
cmocean('curl','zero')
xlabel('Distance (m)','FontSize',24)
ylabel('Pressure','FontSize',24)
title('Expected error ratio: temperature','FontSize',24)
colorbar
clabel(C,h,-10:.1:10)
caxis([0 1])
xlim([0 180*1000])
%brighten(.3)
for i=1:9
    plot([dist(i) dist(i)],[0 4200],'k')
end
for i=1:420
    for j=1:9
        if ~isnan(fake_grid.temp(i,j))
            scatter(dist(j),fake_grid.depth(i,j),'k','filled')
        end
    end
end

figure
hold on
[C,h]=contourf(xgrid,pgrid,analysis_error_s)
axis 'ij'
cmocean('curl','zero')
xlabel('Distance (m)','FontSize',24)
ylabel('Pressure','FontSize',24)
title('Expected error ratio: salinity','FontSize',24)
colorbar
clabel(C,h,-10:.2:10)
caxis([0 1])
xlim([0 180*1000])
brighten(.3)
for i=1:9
    plot([dist(i) dist(i)],[0 4200],'k')
end
for i=1:420
    for j=1:9
        if ~isnan(fake_grid.sal(i,j))
            scatter(dist(j),fake_grid.depth(i,j),'k','filled')
        end
    end
end


figure
hold on
[C,h]=contourf(xgrid,zgrid,.5*(analysis_error_vel(:,:,1)+analysis_error_vel(:,:,2)))
axis 'ij'
cmocean('curl','zero')
xlabel('Distance (m)','FontSize',24)
ylabel('Depth','FontSize',24)
title('Expected error ratio: velocity','FontSize',24)
colorbar
clabel(C,h,-10:.1:10)
caxis([0 1])
xlim([0 180*1000])
brighten(.3)
for i=1:9
    plot([dist(i) dist(i)],[0 4200],'k')
end
for i=1:420
    for j=1:9
        if ~isnan(fake_grid.u(i,j))
            scatter(dist(j),fake_grid.depth(i,j),'k','filled')
        end
    end
end
