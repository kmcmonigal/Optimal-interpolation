% figures for plotting the optimally interpolated data
% also some analysis
% check to make sure none of this is needed in the OI itself 
%% graph ave T, S
load('ACTGEMData_skye.mat')

instr_distance = 1000.*[sw_dist([coast_lat -33.5583],[coast_lon 27.595],'km'),sw_dist([coast_lat -33.6674],[coast_lon 27.6428],'km'),sw_dist([coast_lat -33.7996],[coast_lon 27.7152],'km'),sw_dist([coast_lat -34.0435],[coast_lon 27.8603],'km')];
    
instr_distance_rep = [instr_distance(1,1),instr_distance(1,2),instr_distance(1,2),instr_distance(1,3),instr_distance(1,3),instr_distance(1,3),instr_distance(1,3),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4)];
instr_depth = [250,650,1000,700,1000,1500,2000,700,1000,1500,2000,2500,3000];

micro_distance = [instr_distance(1,2),instr_distance(1,2),instr_distance(1,2),instr_distance(1,2),instr_distance(1,2),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4)];
micro_depth = [750,750,750,1000,1300,750,900,1100,1500];

%dif=nanmean(int_t_OI_clim(1:246,:,:),3)-nanmean(int_t_johns,3);

figure
hold on
contourf(xgrid,pgrid,nanmean(int_t,3))
axis 'ij'
cmocean('thermal')
colorbar
ylim([0 4500])
title('OI temperature')
xlabel('Distance offshore (m)')
ylabel('Pressure (dbar)')
scatter(micro_distance,micro_depth,50,'k','filled','Visible','on')

% dif=nanmean(int_s_OI(1:246,:,:),3)-nanmean(int_s,3);
figure
hold on
contourf(xgrid,pgrid,nanmean(int_s,3))
axis 'ij'
cmocean('haline')
colorbar
ylim([0 4500])
title('OI salinity')
xlabel('Distance offshore (m)')
ylabel('Pressure (dbar)')
scatter(micro_distance,micro_depth,50,'k','filled','Visible','on')

%% do same interpolation with ctd profiles, see how well they match up
% subsample some real profiles
% will have to consider them as totally seperate - basically only thinking
% about vertical interpolation
ctdg.months=month(datetime(ctdg.datenum,'ConvertFrom','datenum'));
for i=235:342 % these are the real full depth profiles from ACT
    micro_fake_d.temp(i-234,1) = ctdg.temp(66,i);
    micro_fake_d.sal(i-234,1) = ctdg.sal(66,i);
    micro_fake_d.pres(i-234,1) = ctdg.pres(66,i);
    micro_fake_d.months(i-234,1) = ctdg.months(i); %650
    
    micro_fake_d.temp(i-234,2) = ctdg.temp(81,i);
    micro_fake_d.sal(i-234,2) = ctdg.sal(81,i);
    micro_fake_d.pres(i-234,2) = ctdg.pres(81,i);
    micro_fake_d.months(i-234,2) = ctdg.months(i); %800
    
    micro_fake_d.temp(i-234,3) = ctdg.temp(101,i);
    micro_fake_d.sal(i-234,3) = ctdg.sal(101,i);
    micro_fake_d.pres(i-234,3) = ctdg.pres(101,i);
    micro_fake_d.months(i-234,3) = ctdg.months(i); %1000
    
    micro_fake_d.temp(i-234,4) = ctdg.temp(151,i);
    micro_fake_d.sal(i-234,4) = ctdg.sal(151,i);
    micro_fake_d.pres(i-234,4) = ctdg.pres(151,i);
    micro_fake_d.months(i-234,4) = ctdg.months(i); %1500
    
    micro_fake_b.temp(i-234,1) = ctdg.temp(71,i);
    micro_fake_b.sal(i-234,1) = ctdg.sal(71,i);
    micro_fake_b.pres(i-234,1) = ctdg.pres(71,i);
    micro_fake_b.months(i-234,1) = ctdg.months(i); %700
    
    micro_fake_b.temp(i-234,2) = ctdg.temp(74,i);
    micro_fake_b.sal(i-234,2) = ctdg.sal(74,i);
    micro_fake_b.pres(i-234,2) = ctdg.pres(74,i);
    micro_fake_b.months(i-234,2) = ctdg.months(i); %730
    
    micro_fake_b.temp(i-234,3) = ctdg.temp(77,i);
    micro_fake_b.sal(i-234,3) = ctdg.sal(77,i);
    micro_fake_b.pres(i-234,3) = ctdg.pres(77,i);
    micro_fake_b.months(i-234,3) = ctdg.months(i); %760
    
    micro_fake_b.temp(i-234,4) = ctdg.temp(91,i);
    micro_fake_b.sal(i-234,4) = ctdg.sal(91,i);
    micro_fake_b.pres(i-234,4) = ctdg.pres(91,i);
    micro_fake_b.months(i-234,4) = ctdg.months(i); %900
    
    micro_fake_b.temp(i-234,5) = ctdg.temp(121,i);
    micro_fake_b.sal(i-234,5) = ctdg.sal(121,i);
    micro_fake_b.pres(i-234,5) = ctdg.pres(121,i);
    micro_fake_b.months(i-234,5) = ctdg.months(i); %1200
end

% do OI on them - temperature first. basically doing 1D interpolation now

zc=171; % 171 for sal, 790 for temp

z_corr_func=@(z) exp(-(z(:)/zc).^2);

B_int.temp=micro_fake_b.temp(:,:);
B_int.sal=micro_fake_b.sal(:,:);
% B_temp_mean=nanmean(B_int.temp,2); % change this - want mean of the 4/5 subsamples at each time
% for i=1:108
%     B_int.temp(i,:)=B_int.temp(i,:)-B_temp_mean(i);
% end

D_int.temp=micro_fake_d.temp(:,:);
D_int.sal=micro_fake_d.sal(:,:);
% D_temp_mean=nanmean(D_int.temp,2);
% for i=1:108
%     D_int.temp(i,:)=D_int.temp(i,:)-D_temp_mean(i);
% end

% try subtracting off climatology
for i=1:108
    for j=1:12
        if ctdg.months(1,234+i)==j
            D_int.temp(i,1)=D_int.temp(i,1)-depth_t_clim(33,j);
            D_int.temp(i,2)=D_int.temp(i,2)-depth_t_clim(40,j);
            D_int.temp(i,3)=D_int.temp(i,3)-depth_t_clim(50,j);
            D_int.temp(i,4)=D_int.temp(i,4)-depth_t_clim(75,j);
            B_int.temp(i,1)=B_int.temp(i,1)-depth_t_clim(35,j);
            B_int.temp(i,2)=B_int.temp(i,2)-depth_t_clim(36,j);
            B_int.temp(i,3)=B_int.temp(i,3)-depth_t_clim(38,j);
            B_int.temp(i,4)=B_int.temp(i,4)-depth_t_clim(45,j);
            B_int.temp(i,5)=B_int.temp(i,5)-depth_t_clim(60,j);
        end
    end
end

for i=1:108
    for j=1:12
        if ctdg.months(1,234+i)==j
            D_int.sal(i,1)=D_int.sal(i,1)-depth_s_clim(33,j);
            D_int.sal(i,2)=D_int.sal(i,2)-depth_s_clim(40,j);
            D_int.sal(i,3)=D_int.sal(i,3)-depth_s_clim(50,j);
            D_int.sal(i,4)=D_int.sal(i,4)-depth_s_clim(75,j);
            B_int.sal(i,1)=B_int.sal(i,1)-depth_s_clim(35,j);
            B_int.sal(i,2)=B_int.sal(i,2)-depth_s_clim(36,j);
            B_int.sal(i,3)=B_int.sal(i,3)-depth_s_clim(38,j);
            B_int.sal(i,4)=B_int.sal(i,4)-depth_s_clim(45,j);
            B_int.sal(i,5)=B_int.sal(i,5)-depth_s_clim(60,j);
        end
    end
end

B_int.z=micro_b.depth(:,:);
D_int.z=micro_d.depth(:,:);

B_temp_var=nanvar(B_int.temp);
D_temp_var=nanvar(D_int.temp);

B_sal_var=nanvar(B_int.sal);
D_sal_var=nanvar(D_int.sal);
% 
% dx_obs=[B_dx;B_dx;B_dx;B_dx;B_dx;D_dx;D_dx;D_dx;D_dx];
% for i=1:108
%     dz_obs(:,i)=[B_int.z(i,:).';D_int.z(i,:).'];
%     temp_obs(:,i)=[B_int.temp(i,:).';D_int.temp(i,:).'];
% end
      
ratio_obs_b=noise_micro_b(:,2).'./B_sal_var; %change index of noise from 1 (temp) to 2 (sal)
ratio_obs_d=noise_micro_d(:,2).'./D_sal_var;

ratio_b=zeros(5,5);
for i=1:5
    for j=1:5
        if i==j
            ratio_b(i,j)=ratio_obs_b(i);
        end
    end
end

ratio_d=zeros(4,4);
for i=1:4
    for j=1:4
        if i==j
            ratio_d(i,j)=ratio_obs_d(i);
        end
    end
end

zgrid=0:20:5000;

% no NaN's - should be simpler than the velocity data
% do for mooring B set up first
clear anom_value
clear weight_corr
clear cross_corr
clear weights
for time=1:108 %only doing 1D - will need to modify
    for i=1:length(ratio_b)
        for j=1:length(zgrid)
            weight_corr(i,j)=z_corr_func(abs(zgrid(j)-B_int.z(i)));
        end
    end
    
    for i=1:length(ratio_b)
        for j=1:length(ratio_b)
            cross_corr(i,j)=z_corr_func(abs(B_int.z(i)-B_int.z(j)));
        end
    end

    for j=1:length(zgrid)
            weights(:,j)=(ratio_b+cross_corr)\weight_corr(:,j);
    end
    
    for j=1:length(zgrid)
            anom_value_fake_b(j,time)=weights(:,j).'*B_int.sal(time,:).'; %check to see if dimensions work?
    end
end

% now add the anom values to the climatology and see how well we did
for i=1:251
    for k=1:108
        for l=1:12
            if ctdg.months(i+234)==l
                fake_t_d(i,k)=anom_value_fake_b(i,k)+depth_t_clim(i,l);
                fake_t_b(i,k)=anom_value_fake_d(i,k)+depth_t_clim(i,l);
            end
        end % WHY ARE TEMPS SO HIGH. SOMETHING ISN'T WORKING
    end %pretend like all in jan
end

% same thing for salinity
for i=1:251
    for k=1:108
        for l=1:12
            if ctdg.months(i+234)==l
                fake_s_d(i,k)=anom_value_fake_b(i,k)+depth_s_clim(i,l);
                fake_s_b(i,k)=anom_value_fake_d(i,k)+depth_s_clim(i,l);
            end
        end % WHY ARE TEMPS SO HIGH. SOMETHING ISN'T WORKING
    end %pretend like all in jan
end

% depth vs pressure issue again - what pressure does each depth match with
% - want depth_grid that is closest to dz - WILL NEED TO AMEND BUT LETS GET
% A BASIC IDEA
ctdg.depth=gsw_depth_from_z(gsw_z_from_p(ctdg.pgrid,-34));
dz=(0:20:5000).';

for i=1:length(ctdg.depth)
    for j=1:length(dz)
        dif(i,j)=abs(ctdg.depth(i)-dz(j));
        % want to get index of min difference at each dz (251 of them)
    end
end
[M,I]=min(dif);

% find closest  value to 20 m depth bins

for i=1:108
    for j=1:251
        profile_dif.temp(j,i) = fake_s_d(j,i)-ctdg.sal(I(j),i+234);
        rms(j,i) = sqrt((fake_s_d(j,i)-ctdg.sal(I(j),i+234)).^2);
    end
end

for j=1:251
    ave_dif(j) = nanmean(profile_dif.temp(j,:));
    std_dif(j) = nanstd(profile_dif.temp(j,:));
    rms_ave(j) = nanmean(rms(j,:)); % not 100% sure this is right
end

%figure
% seems weird that it's not even good at the subsampled depths. make sure
% i'm doing it right......
% I think there is an issue with adding back in the mean. I should add back
% in the climatology of the subsampled depths, to all depths... but then
% only like 3 degree temp range between top and bottom. Hmm
% maybe should remove a climatology, then add climatology back. let's try
% that - ok, much better, now just fix depth/pres issue. also redo the T/S
% from the time series and graph time mean of those
figure
hold on
plot(ave_dif(:),0:20:5000)
plot(rms_ave(:),0:20:5000)
% plot(std_dif(:),overall_profile_sub.pres(:,1))
axis 'ij'
x=[0,0];
y=[0,5000];
plot(x,y,'k')
x1=[-.01,-.01];
x2=[.01,.01];
plot(x1,y,'k',x2,y,'k')
title('Mooring D interpolation sensitivity - salinity') % seems like a really high degree.... try lower
xlabel('Salinity')
ylabel('Depth (m)')
legend('Mean difference','RMS') %,'Standard deviation')
ylim([0 3700])

%% look at variability
for i=1:182
    for j=1:246
        std_johns(j,i)=nanstd(int_t_johns(j,i,:));
        std_OI(j,i)=nanstd(int_t_OI(j,i,:));
    end
end

dif=std_johns-std_OI;

figure
hold on
contourf(xgrid(1:246,:),zgrid(1:246,:),dif)
axis 'ij'
colorbar
cmocean('curl')
ylim([0 4500])
title('Standard dev OI temperature - Johns temperature')
xlabel('Distance offshore (m)')
ylabel('Depth (m)')
scatter(micro_distance,micro_depth,50,'k','filled','Visible','on')
caxis([-1.2 1.2])

%% correlation???

for i=1:182
    for j=1:246
        corr_johns_OI(j,i)=corr(int_t_johns(j,i,:),int_t_OI(j,i,:));
    end
end

% doesn't work, just look at some time series
for i=1:722
series_1(i)=int_t_johns(50,120,i);
series_2(i)=int_t_OI(50,120,i);
end

figure
hold on
plot(1:722,series_1)
plot(1:722,series_2)

%% velocity figures
%% figure of time mean velocities


distance_offshore = 0:500:1000.*sw_dist([coast_lat -34.0435],[coast_lon 27.8603],'km'); %grr... should be 184 long... but middle point shared?

for k=1:251
    depth(k) = 20*(k-1);
end

[distance_offshore_grid,depth_grid] = meshgrid(distance_offshore,depth); 

% add in topography
% nice color bar
my_map = [0,0,255;
50,255,255;
255,233,0;
253,200,0;
253,141,0;
251,48,0;
233,2,0;
201,1,0;
168,0,0;
144,0,0;
123,0,0;
101,0,0];

for i=1:12
    for j=1:3
        my_map2(i,j) = my_map(13-i,j);
    end
end

% a_pos = [27.595,-33.5583];
% b_pos = [27.6428,-33.6674];
% c_pos = [27.7152,-33.7996];
% d_pos = [27.8603,-34.0435];

instr_distance = 1000.*[sw_dist([coast_lat -33.5583],[coast_lon 27.595],'km'),sw_dist([coast_lat -33.6674],[coast_lon 27.6428],'km'),sw_dist([coast_lat -33.7996],[coast_lon 27.7152],'km'),sw_dist([coast_lat -34.0435],[coast_lon 27.8603],'km')];
    
instr_distance_rep = [instr_distance(1,1),instr_distance(1,2),instr_distance(1,2),instr_distance(1,3),instr_distance(1,3),instr_distance(1,3),instr_distance(1,3),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4)];
instr_depth = [250,650,1000,700,1000,1500,2000,700,1000,1500,2000,2500,3000];

micro_distance = [instr_distance(1,2),instr_distance(1,2),instr_distance(1,2),instr_distance(1,2),instr_distance(1,2),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4)];
micro_depth = [750,750,750,1000,1300,750,900,1100,1500];
%meshgrid(distance_offshore,depth)
figure
hold on
contourf(distance_offshore./10^3,depth,nanmean(v_vel,3),'ShowText','on')
axis 'ij'
colorbar
xlabel('Distance (km)')
ylabel('Depth (m)')
caxis([-2 .4])
colormap(my_map2./255)
scatter(instr_distance_rep./10^3,instr_depth,50,'k','filled','Visible','on') %why can't I just
%plot the scatter on top of the contour???
%scatter(micro_distance./10^3,micro_depth,50,'m','filled','Visible','on')

% definitely not right. let's try with the time mean subtracted instead
% is it possible I rotated wrong?
%% calculate transport, figure of transport time series

%% look at error matrix produced
