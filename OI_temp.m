% optimally interpolate T and S using similar method as OI_vel

%% OKAY let's redo this, using similar techniques to OI_vel - TEMP here
% load coas_lat, coast_lon, micro_b_noise, micro_d_noise,
% micro_b_d_despike_filter
%% create climatology

ctdg.months = month(datetime(ctdg.datenum,'ConvertFrom','datenum'));

for k=1:501
    for j=1:2415
        for i=1:12
            if ctdg.months(1,j) == i
                monthly_clim_t(k,j,i) = ctdg.temp(k,j);
                monthly_clim_s(k,j,i) = ctdg.sal(k,j);
            end
        end
    end
end

monthly_clim_t(monthly_clim_t==0)=NaN;
monthly_clim_s(monthly_clim_s==0)=NaN;

for i=1:12
    for j=1:501
        ave_temp(j,i) = nanmean(monthly_clim_t(j,:,i));
        number_measures(j,i) = sum(~isnan(monthly_clim_t(j,:,i)));
        ave_sal(j,i) = nanmean(monthly_clim_s(j,:,i));
    end
end

% Below 2000m, make ave_temp(:,i) be the same for each month
for i=1:12
    for j=1:200
        ave_temp2(j,i) = ave_temp(j,i);
        ave_sal2(j,i) = ave_sal(j,i);
    end
    for j=201:501
        ave_temp2(j,i) = nanmean(ave_temp(j,:));
        ave_sal2(j,i) = nanmean(ave_sal(j,:));
    end
end

% issue - in pressure coordinates
depth_grid=gsw_depth_from_z(gsw_z_from_p(ctdg.pgrid,-34));
dz=(0:20:5000).';

% change from 10 dbar pressure bins to 20 m depth bins
depth_t_clim=zeros(251,12);
depth_s_clim=zeros(251,12);
count=zeros(251,12);
for k=1:12
    for i=1:251
        for j=1:501
            if depth_grid(j)>(i-1)*20 && depth_grid(j)<i*20
                depth_t_clim(i,k)=depth_t_clim(i,k)+ave_temp2(j,k);
                count(i,k)=count(i,k)+1;
                depth_s_clim(i,k)=depth_s_clim(i,k)+ave_sal2(j,k);
            end
        end
    end
end

for k=1:12
    for i=1:251
        depth_t_clim(i,k)=depth_t_clim(i,k)./count(i,k);
        depth_s_clim(i,k)=depth_s_clim(i,k)./count(i,k);
    end
end

%%

ctdg.depth=gsw_depth_from_z(gsw_z_from_p(ctdg.pgrid,-34));

b_pos = [27.6428,-33.6674];
d_pos = [27.8603,-34.0435];

B_dx=1000*sw_dist([coast_lat b_pos(2)],[coast_lon b_pos(1)],'km');
D_dx=1000*sw_dist([coast_lat d_pos(2)],[coast_lon d_pos(1)],'km');

xc=130*1000; %amend to somewhere between 50, 53
zc=790; %amend to somewhere between 1900,2600

x_corr_func=@(x) exp(-(x(:)/xc).^2).*cos(pi.*x(:)./(2.*xc));
z_corr_func=@(z) exp(-(z(:)/zc).^2);

B_int.temp=micro_b.temp(:,:);
D_int.temp=micro_d.temp(6:727,:); 

% mean_b=nanmean(B_int.temp);
% mean_d=nanmean(D_int.temp);

% for i=1:722
%     B_int.temp(i,:)=B_int.temp(i,:)-mean_b;
%     D_int.temp(i,:)=D_int.temp(i,:)-mean_d;
% end
% obs=[B_int.temp,D_int.temp];
% space_mean_t=nanmean(obs,2);
% for i=1:722 % do space mean at each time instead
%     B_int.temp(i,:)=B_int.temp(i,:)-space_mean_t(i);
%     D_int.temp(i,:)=D_int.temp(i,:)-space_mean_t(i);
% end

%how did it still work despite this?
% want to subtract off climatological value - at each time step find clim
% pres closest to micro depth, subtract that off
B_int.z=micro_b.depth(:,:);
D_int.z=micro_d.depth(6:727,:);

dist=nan(722,501,5); %this section is to subtract clim mean
for i=1:722
    for j=1:length(ctdg.depth)
        for l=1:5
            dist(i,j,l)=abs(ctdg.depth(j,1)-B_int.z(i,l)); %why does this loop take so long? preallocate, try again - also, is it what I want???
        end
    end
end

clear M
clear I
for i=1:722
    [M(i,:),I(i,:)]=min(dist(i,:,:));
end

B_int.tempanom=nan(722,5);
B_int.months=month(micro_b.date);
for i=1:722
    for j=1:5
        for k=1:12
            if B_int.months(i,j)==k
                B_int.tempanom(i,j)=B_int.temp(i,j)-depth_t_clim(I(i,j),k);
            end
        end
    end
end

% same for D - subtract clim temp
dist=nan(722,501,4);
for i=1:722
    for j=1:length(ctdg.depth)
        for l=1:4
            dist(i,j,l)=abs(ctdg.depth(j,1)-D_int.z(i,l)); %why does this loop take so long? preallocate, try again - also, is it what I want???
        end
    end
end

clear M
clear I
for i=1:722
    [M(i,:),I(i,:)]=min(dist(i,:,:));
end

D_int.months=month(micro_d.date(6:727,:));
for i=1:722
    for j=1:4
        for k=1:12
            if D_int.months(i,j)==k
                D_int.tempanom(i,j)=D_int.temp(i,j)-depth_t_clim(I(i,j),k);
            end
        end
    end
end

B_temp_var=nanvar(B_int.temp);
D_temp_var=nanvar(D_int.temp);

dx_obs=[B_dx;B_dx;B_dx;B_dx;B_dx;D_dx;D_dx;D_dx;D_dx];
for i=1:722
    dz_obs(:,i)=[B_int.z(i,:).';D_int.z(i,:).'];
    temp_obs(:,i)=[B_int.tempanom(i,:).';D_int.tempanom(i,:).'];
end
        
var_obs=[B_temp_var.';D_temp_var.'];
noise_obs=[noise_micro_b(:,1);noise_micro_d(:,1)];

ratio_obs=noise_obs./var_obs;

ratio=zeros(9,9);
for i=1:9
    for j=1:9
        if i==j
            ratio(i,j)=ratio_obs(i);
        end
    end
end


x=0:500:1000*sw_dist([coast_lat d_pos(2)],[coast_lon d_pos(1)],'km');
z=0:20:5000; % then later will make anything under topography into NaN
for i=1:length(z)
    xgrid(i,:)=x(:);
end
for i=1:length(x)
    zgrid(:,i)=z(:);
end

% check to see B_int.temp values make sense
% check to see that anom + clim is what we are getting

% no NaN's - should be simpler than the velocity data
clear weight_corr
clear cross_corr
clear weights
clear anom_value
for time=1:722
    for i=1:length(ratio)
        for j=1:length(zgrid)
            for k=1:size(zgrid,2)
                weight_corr(i,j,k)=x_corr_func(abs(xgrid(j,k)-dx_obs(i)))*z_corr_func(abs(zgrid(j,k)-dz_obs(i)));
            end
        end
    end
    
    for i=1:length(ratio)
        for j=1:length(ratio)
            cross_corr(i,j)=x_corr_func(abs(dx_obs(i)-dx_obs(j)))*z_corr_func(abs(dz_obs(i)-dz_obs(j)));
        end
    end

    for j=1:length(zgrid)
        for k=1:size(zgrid,2)
            weights(:,j,k)=(ratio+cross_corr)\weight_corr(:,j,k);
        end
    end
    
    for j=1:length(zgrid)
        for k=1:size(zgrid,2)
            anom_value(j,k,time)=weights(:,j,k).'*temp_obs(:,time); %check to see if dimensions work?
        end
    end
end


%% same thing with salinity
% make sure the temp_obs is after subtracting the climatology!!!!

b_pos = [27.6428,-33.6674];
d_pos = [27.8603,-34.0435];

B_dx=1000*sw_dist([coast_lat b_pos(2)],[coast_lon b_pos(1)],'km');
D_dx=1000*sw_dist([coast_lat d_pos(2)],[coast_lon d_pos(1)],'km');

xc=77*1000; 
zc=171; 

x_corr_func=@(x) exp(-(x(:)/xc).^2).*cos(pi.*x(:)./(2.*xc));
z_corr_func=@(z) exp(-(z(:)/zc).^2);

B_int.temp=micro_b.sal(:,:);
D_int.temp=micro_d.sal(6:727,:);

clear obs
obs=[B_int.temp,D_int.temp];
space_mean_s=nanmean(obs,2);
for i=1:722
    B_int.temp(i,:)=B_int.temp(i,:)-space_mean_s(i);
    D_int.temp(i,:)=D_int.temp(i,:)-space_mean_s(i);
end
% B_s_mean=nanmean(B_int.temp); % this is the time mean
% D_s_mean=nanmean(D_int.temp);
% 
% for i=1:722
%     B_int.temp(i,:)=B_int.temp(i,:)-B_s_mean;
%     D_int.temp(i,:)=D_int.temp(i,:)-D_s_mean;
% end

% this section subtracts off clim mean
% dist=nan(722,501,5);
% for i=1:722
%     for j=1:length(ctdg.depth)
%         for l=1:5
%             dist(i,j,l)=abs(ctdg.depth(j,1)-B_int.z(i,l)); %why does this loop take so long? preallocate, try again - also, is it what I want???
%         end
%     end
% end
% 
% clear M
% clear I
% for i=1:722
%     [M(i,:),I(i,:)]=min(dist(i,:,:));
% end
% 
% B_int.months=month(micro_b.date);
% for i=1:722
%     for j=1:5
%         for k=1:12
%             if B_int.months(i,1)==k
%                 B_int.temp(i,j)=B_int.temp(i,j)-depth_s_clim(I(i,j),k);
%             end
%         end
%     end
% end
% 
% % same for D - subtract clim temp
% dist=nan(722,501,4);
% for i=1:722
%     for j=1:length(ctdg.depth)
%         for l=1:4
%             dist(i,j,l)=abs(ctdg.depth(j,1)-D_int.z(i,l)); %why does this loop take so long? preallocate, try again - also, is it what I want???
%         end
%     end
% end
% 
% clear M
% clear I
% for i=1:722
%     [M(i,:),I(i,:)]=min(dist(i,:,:));
% end
% 
% D_int.months=month(micro_d.date(6:727,:));
% for i=1:722
%     for j=1:4
%         for k=1:12
%             if D_int.months(i,j)==k
%                 D_int.temp(i,j)=D_int.temp(i,j)-depth_s_clim(I(i,j),k);
%             end
%         end
%     end
% end
B_int.z=micro_b.depth(:,:);
D_int.z=micro_d.depth(4:727,:);

B_temp_var=nanvar(B_int.temp);
D_temp_var=nanvar(D_int.temp);

dx_obs=[B_dx;B_dx;B_dx;B_dx;B_dx;D_dx;D_dx;D_dx;D_dx];
for i=1:722
    dz_obs(:,i)=[B_int.z(i,:).';D_int.z(i,:).'];
    temp_obs(:,i)=[B_int.temp(i,:).';D_int.temp(i,:).'];
end
        
var_obs=[B_temp_var.';D_temp_var.'];
noise_obs=[noise_micro_b(:,4);noise_micro_d(:,4)];

ratio_obs=noise_obs./var_obs;

ratio=zeros(9,9);
for i=1:9
    for j=1:9
        if i==j
            ratio(i,j)=ratio_obs(i);
        end
    end
end


x=0:500:1000*sw_dist([coast_lat d_pos(2)],[coast_lon d_pos(1)],'km');
z=0:20:5000; % then later will make anything under topography into NaN
for i=1:length(z)
    xgrid(i,:)=x(:);
end
for i=1:length(x)
    zgrid(:,i)=z(:);
end


% no NaN's - should be simpler than the velocity data
clear weight_corr
clear cross_corr
clear weights
clear anom_value
for time=1:722
    for i=1:length(ratio)
        for j=1:length(zgrid)
            for k=1:size(zgrid,2)
                weight_corr(i,j,k)=x_corr_func(abs(xgrid(j,k)-dx_obs(i)))*z_corr_func(abs(zgrid(j,k)-dz_obs(i)));
            end
        end
    end
    
    for i=1:length(ratio)
        for j=1:length(ratio)
            cross_corr(i,j)=x_corr_func(abs(dx_obs(i)-dx_obs(j)))*z_corr_func(abs(dz_obs(i)-dz_obs(j)));
        end
    end

    for j=1:length(zgrid)
        for k=1:size(zgrid,2)
            weights(:,j,k)=(ratio+cross_corr)\weight_corr(:,j,k);
        end
    end
    
    for j=1:length(zgrid)
        for k=1:size(zgrid,2)
            anom_value(j,k,time)=weights(:,j,k).'*temp_obs(:,time); %check to see if dimensions work?
        end
    end
end

%% add the anomaly values to a climatology
% % add anomalies to climatology
% % first need to get months for anomaly field
micro_months=month(micro_b.date);

for i=1:251
    for j=1:182
        for k=1:722
            for l=1:12
                if micro_months(k,1)==l
                    int_t(i,j,k)=anom_value_clim2(i,j,k)+depth_t_clim(i,l);
                    %int_s(i,j,k)=anom_value_sal3(i,j,k)+depth_s_clim(i,l);
                end
            end
        end
    end
end

%% make parts under topography NaN

% a_b=1000*sw_dist([-33.5583 -33.6674],[27.595 27.6428],'km');
% b_c=1000*sw_dist([-33.6674 -33.7996],[27.6428 27.7152],'km');
% c_d=1000*sw_dist([-33.7996 -34.0435],[27.7152 27.8603],'km');
% 
% a_b_lat = interp1([0,a_b],[-33.5583,-33.6674],[0:500:a_b]);
% a_b_lon = interp1([0,a_b],[27.595,27.6428],[0:500:a_b]);
% 
% b_c_lat = interp1([a_b,b_c+a_b],[-33.6674,-33.7996],[a_b:500:b_c+a_b]);
% b_c_lon = interp1([a_b,b_c+a_b],[27.6428,27.7152],[a_b:500:b_c+a_b]);
% 
% c_d_lat = interp1([b_c+a_b,c_d+a_b+b_c],[-33.7996,-34.0435],[b_c+a_b:500:c_d+a_b+b_c]);
% c_d_lon = interp1([b_c+a_b,c_d+a_b+b_c],[27.7152,27.8603],[b_c+a_b:500:c_d+a_b+b_c]);

overall_lat = interp1([0,D_dx],[coast_lat,d_pos(2)],[0:500:D_dx]); 
overall_lon = interp1([0,D_dx],[coast_lon,d_pos(1)],[0:500:D_dx]);

% now compare our interpolated points to ETOPO1 and find the closest match
% to each point
[elev,long,lat]=m_etopo2([27 28 -34 -33]);

dist = nan(182,61,61);
for i=1:182
    for j=1:61
        for k=1:61
            dist(i,j,k) = sw_dist([overall_lat(1,i) lat(j,k)],[overall_lon(1,i) long(j,k)],'km');
        end % why overall_lat, overall_lon shorter than my interpolated values
    end
end

clear M
clear I
closest = nan(182,1);
for i=1:182
    min_temp = dist(i,:,:);
    temp = min_temp(:);
    [M,I] = min(temp);
    closest(i,1) = I;
end

for i=1:182
    [I_row(i),I_col(i)] = ind2sub([61,61],closest(i));
end

%sw_dist([overall_lat(1,2) lat(I_row(1,2),I_col(1,2))],[overall_lon(1,2) long(I_row(1,2),I_col(1,2))],'km');
% yay it works. Now need to get the topography at each of those points, in
% order

topo = nan(182,1);
for i=1:182
    topo(i,1) = elev(I_row(i),I_col(i));
end
% doesn't match perfectly - high resolution topo data from ACT would be
% better
depth_int = nan(251,725,182);
for i=1:725
    for j=1:182
        for k=1:251
            depth_int(k,i,j) = 20*(k-1);
        end
    end
end

for i=1:722
    for j=1:182
        for k=1:251
            if depth_int(k,i,j) > abs(topo(j,1))
                int_s(k,j,i) = NaN;
                int_t(k,j,i) = NaN; % this deleted unphysical depths, still need to extrapolate inland...
            end
        end
    end
end

% I just did OI on whole domain - may be better to only use OI between
% measurement points then use no slip or something to coast
%% graph ave T, S
instr_distance = 1000.*[sw_dist([coast_lat -33.5583],[coast_lon 27.595],'km'),sw_dist([coast_lat -33.6674],[coast_lon 27.6428],'km'),sw_dist([coast_lat -33.7996],[coast_lon 27.7152],'km'),sw_dist([coast_lat -34.0435],[coast_lon 27.8603],'km')];
    
instr_distance_rep = [instr_distance(1,1),instr_distance(1,2),instr_distance(1,2),instr_distance(1,3),instr_distance(1,3),instr_distance(1,3),instr_distance(1,3),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4)];
instr_depth = [250,650,1000,700,1000,1500,2000,700,1000,1500,2000,2500,3000];

micro_distance = [instr_distance(1,2),instr_distance(1,2),instr_distance(1,2),instr_distance(1,2),instr_distance(1,2),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4),instr_distance(1,4)];
micro_depth = [750,750,750,1000,1300,750,900,1100,1500];

%dif=nanmean(int_t_OI_clim(1:246,:,:),3)-nanmean(int_t_johns,3);

figure
hold on
contourf(xgrid,zgrid,nanmean(int_t,3))
axis 'ij'
cmocean('thermal')
colorbar
ylim([0 4500])
title('OI temperature')
xlabel('Distance offshore (m)')
ylabel('Depth (m)')
scatter(micro_distance,micro_depth,50,'k','filled','Visible','on')

dif=nanmean(int_s_OI(1:246,:,:),3)-nanmean(int_s_johns,3);
figure
hold on
contourf(xgrid(1:246,:),zgrid(1:246,:),dif)
axis 'ij'
cmocean('haline')
colorbar
ylim([0 4500])
title('OI-Johns salinity')
xlabel('Distance offshore (m)')
ylabel('Depth (m)')
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


