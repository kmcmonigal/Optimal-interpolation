% Do OI on the U velocities anomalies - then will need to determine a
% background to add to it

% load filter sub rotate
a_pos = [27.595,-33.5583];
b_pos = [27.6428,-33.6674];
c_pos = [27.7152,-33.7996];
d_pos = [27.8603,-34.0435];

A_int.v=A.vtmp3_filt(:,1:725);
B_int.v=B.vtmp4_filt(:,:);
C_int.v=C.vtmp4_filt(:,2:726);
D_int.v=D.vtmp4_filt(:,3:727);

all_obs=[A_int.v;B_int.v;C_int.v;D_int.v]; % shit - need to rerun with [] (check)
mean_obs=nanmean(all_obs);
% No - want to remove mean of all instruments at each time step
% this does the space mean. let's do time mean
% for i=1:725
%     A_int.v(:,i)=A_int.v(:,i)-mean_obs(i);
%     B_int.v(:,i)=B_int.v(:,i)-mean_obs(i);
%     C_int.v(:,i)=C_int.v(:,i)-mean_obs(i);
%     D_int.v(:,i)=D_int.v(:,i)-mean_obs(i);
% end
A_mean=nanmean(A_int.v);
B_mean=nanmean(B_int.v);
C_mean=nanmean(C_int.v);
D_mean=nanmean(D_int.v);

for i=1:725
    A_int.v(:,i)=A_int.v(:,i)-A_mean(i); %time mean
    B_int.v(:,i)=B_int.v(:,i)-B_mean(i);
    C_int.v(:,i)=C_int.v(:,i)-C_mean(i);
    D_int.v(:,i)=D_int.v(:,i)-D_mean(i);
end


A_int.z=A.ztmp3(:,1:725);
B_int.z=B.ztmp4(:,:);
C_int.z=C.ztmp4(:,2:726);
D_int.z=D.ztmp4(:,3:727);

x=0:500:1000*sw_dist([coast_lat d_pos(2)],[coast_lon d_pos(1)],'km');
z=0:20:5000; % then later will make anything under topography into NaN
for i=1:length(z)
    xgrid(i,:)=x(:);
end
for i=1:length(x)
    zgrid(:,i)=z(:);
end

A_dx=1000*sw_dist([coast_lat a_pos(2)],[coast_lon a_pos(1)],'km');
B_dx=1000*sw_dist([coast_lat b_pos(2)],[coast_lon b_pos(1)],'km');
C_dx=1000*sw_dist([coast_lat c_pos(2)],[coast_lon c_pos(1)],'km');
D_dx=1000*sw_dist([coast_lat d_pos(2)],[coast_lon d_pos(1)],'km');

xc=50*1000; %amend to somewhere between 50, 53
zc=2200; %amend to somewhere between 1900,2600

x_corr_func=@(x) exp(-(x(:)/xc).^2).*cos(pi.*x(:)./(2.*xc));
z_corr_func=@(z) exp(-(z(:)/zc).^2);


% need to do Noise matrix,ratio with actual variance
for time=535:725
    clear Noise2
    clear u_obs_time
    clear dx
    clear weights
    clear cross_corr
    clear ratio
    clear weight_corr
    u_obs_time=[A_int.v(~isnan(A_int.v(:,time)),:);B_int.v(~isnan(B_int.v(:,time)),:);C_int.v(~isnan(C_int.v(:,time)),:);D_int.v(~isnan(D_int.v(:,time)),:)];
    dz_obs=[A_int.z(~isnan(A_int.v(:,time)),time);B_int.z(~isnan(B_int.v(:,time)),time);C_int.z(~isnan(C_int.v(:,time)),time);D_int.z(~isnan(D_int.v(:,time)),time)];
    u_obs=[A_int.v(~isnan(A_int.v(:,time)),time);B_int.v(~isnan(B_int.v(:,time)),time);C_int.v(~isnan(C_int.v(:,time)),time);D_int.v(~isnan(D_int.v(:,time)),time)];
    % setting up Noise2 matrix
    for i=1:sum(~isnan(A_int.v(:,time)))
        Noise2(i)=Noise(1);
    end
    for i=sum(~isnan(A_int.v(:,time)))+1:sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))-1
        Noise2(i)=Noise(2);
    end
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time))))=Noise(3);
    for i=sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+1:sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))-3
        Noise2(i)=Noise(4);
    end
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))-2)=Noise(5);
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))-1)=Noise(6);
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time))))=Noise(7);
    for i=sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+1:sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+sum(~isnan(D_int.v(:,time)))-5
        Noise2(i)=Noise(8);
    end
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+sum(~isnan(D_int.v(:,time)))-4)=Noise(9);
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+sum(~isnan(D_int.v(:,time)))-3)=Noise(10);
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+sum(~isnan(D_int.v(:,time)))-2)=Noise(11);
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+sum(~isnan(D_int.v(:,time)))-1)=Noise(12);
    Noise2(sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+sum(~isnan(D_int.v(:,time))))=Noise(13);
    % set up ratio matrix
    ratio=zeros(length(Noise2),length(Noise2));
    Noise2=Noise2/100;
    for i=1:length(Noise2)
        for j=1:length(Noise2)
            if i==j
                ratio(i,j)=Noise2(i)/nanvar(u_obs_time(i,:)); % seems bad that these values are so big. is this right? what was it on the one obs case?
            end
        end
    end
    % set up dx matrix
    for i=1:sum(~isnan(A_int.v(:,time)))
        dx(i)=A_dx;
    end
    for i=sum(~isnan(A_int.v(:,time)))+1:sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))
        dx(i)=B_dx;
    end
    for i=sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+1:sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))
        dx(i)=C_dx;
    end
    for i=sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+1:sum(~isnan(A_int.v(:,time)))+sum(~isnan(B_int.v(:,time)))+sum(~isnan(C_int.v(:,time)))+sum(~isnan(D_int.v(:,time)))
        dx(i)=D_dx;
    end
    % now get cross correlations between instruments, grid points
    for i=1:length(Noise2)
        for j=1:length(zgrid)
            for k=1:size(zgrid,2)
                weight_corr(i,j,k)=x_corr_func(abs(xgrid(j,k)-dx(i)))*z_corr_func(abs(zgrid(j,k)-dz_obs(i)));
            end
        end
    end
    % get cross corr between instruments - should this be theoretical or
    % "real" though
    for i=1:length(Noise2)
        for j=1:length(Noise2)
            cross_corr(i,j)=x_corr_func(abs(dx(i)-dx(j)))*z_corr_func(abs(dz_obs(i)-dz_obs(j)));
        end
    end
    % solve for weights
    for j=1:length(zgrid)
        for k=1:size(zgrid,2)
            weights(:,j,k)=(ratio+cross_corr)\weight_corr(:,j,k); % may need to do transform of weight_corr
        end
    end

    for j=1:length(zgrid)
        for k=1:size(zgrid,2)
            anom_value(j,k,time)=weights(:,j,k).'*u_obs; % mean already subtracted
        end
    end
    
end

%% do for u

A_int.u=A.utmp3_filt(:,1:725);
B_int.u=B.utmp4_filt(:,:);
C_int.u=C.utmp4_filt(:,2:726);
D_int.u=D.utmp4_filt(:,3:727);

all_obs=[A_int.u;B_int.u;C_int.u;D_int.u];
% mean_obs=nanmean(all_obs);
% No - want to remove mean of all instruments at each time step

% for i=1:725
%     A_int.u(:,i)=A_int.u(:,i)-mean_obs(i); %space mean - let's try time
%     B_int.u(:,i)=B_int.u(:,i)-mean_obs(i);
%     C_int.u(:,i)=C_int.u(:,i)-mean_obs(i);
%     D_int.u(:,i)=D_int.u(:,i)-mean_obs(i);
% end

A_mean_u=nanmean(A_int.u);
B_mean_u=nanmean(B_int.u);
C_mean_u=nanmean(C_int.u);
D_mean_u=nanmean(D_int.u);

for i=1:725
    A_int.u(:,i)=A_int.u(:,i)-A_mean_u(i);
    B_int.u(:,i)=B_int.u(:,i)-B_mean_u(i);
    C_int.u(:,i)=C_int.u(:,i)-C_mean_u(i);
    D_int.u(:,i)=D_int.u(:,i)-D_mean_u(i);
end

for time=1:725
    clear Noise2
    clear u_obs_time
    clear dx
    clear weights
    clear cross_corr
    clear ratio
    clear weight_corr
    u_obs_time=[A_int.u(~isnan(A_int.u(:,time)),:);B_int.u(~isnan(B_int.u(:,time)),:);C_int.u(~isnan(C_int.u(:,time)),:);D_int.u(~isnan(D_int.u(:,time)),:)];
    dz_obs=[A_int.z(~isnan(A_int.u(:,time)),time);B_int.z(~isnan(B_int.u(:,time)),time);C_int.z(~isnan(C_int.u(:,time)),time);D_int.z(~isnan(D_int.u(:,time)),time)];
    u_obs=[A_int.u(~isnan(A_int.u(:,time)),time);B_int.u(~isnan(B_int.u(:,time)),time);C_int.u(~isnan(C_int.u(:,time)),time);D_int.u(~isnan(D_int.u(:,time)),time)];
    % setting up Noise2 matrix
    for i=1:sum(~isnan(A_int.u(:,time)))
        Noise2(i)=Noise(1);
    end
    for i=sum(~isnan(A_int.u(:,time)))+1:sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))-1
        Noise2(i)=Noise(2);
    end
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time))))=Noise(3)
    for i=sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+1:sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))-3
        Noise2(i)=Noise(4);
    end
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))-2)=Noise(5);
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))-1)=Noise(6);
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time))))=Noise(7);
    for i=sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+1:sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+sum(~isnan(D_int.u(:,time)))-5
        Noise2(i)=Noise(8);
    end
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+sum(~isnan(D_int.u(:,time)))-4)=Noise(9);
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+sum(~isnan(D_int.u(:,time)))-3)=Noise(10);
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+sum(~isnan(D_int.u(:,time)))-2)=Noise(11);
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+sum(~isnan(D_int.u(:,time)))-1)=Noise(12);
    Noise2(sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+sum(~isnan(D_int.u(:,time))))=Noise(13);
    % set up ratio matrix
    ratio=zeros(length(Noise2),length(Noise2));
    Noise2=Noise2/100;
    for i=1:length(Noise2)
        for j=1:length(Noise2)
            if i==j
                ratio(i,j)=Noise2(i)/nanvar(u_obs_time(i,:)); % seems bad that these values are so big. is this right? what was it on the one obs case?
            end
        end
    end
    % set up dx matrix
    for i=1:sum(~isnan(A_int.u(:,time)))
        dx(i)=A_dx;
    end
    for i=sum(~isnan(A_int.u(:,time)))+1:sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))
        dx(i)=B_dx;
    end
    for i=sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+1:sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))
        dx(i)=C_dx;
    end
    for i=sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+1:sum(~isnan(A_int.u(:,time)))+sum(~isnan(B_int.u(:,time)))+sum(~isnan(C_int.u(:,time)))+sum(~isnan(D_int.u(:,time)))
        dx(i)=D_dx;
    end
    % now get cross correlations between instruments, grid points
    for i=1:length(Noise2)
        for j=1:length(zgrid)
            for k=1:size(zgrid,2)
                weight_corr(i,j,k)=x_corr_func(abs(xgrid(j,k)-dx(i)))*z_corr_func(abs(zgrid(j,k)-dz_obs(i)));
            end
        end
    end
    % get cross corr between instruments - should this be theoretical or
    % "real" though
    for i=1:length(Noise2)
        for j=1:length(Noise2)
            cross_corr(i,j)=x_corr_func(abs(dx(i)-dx(j)))*z_corr_func(abs(dz_obs(i)-dz_obs(j)));
        end
    end
    % solve for weights
    for j=1:length(zgrid)
        for k=1:size(zgrid,2)
            weights(:,j,k)=(ratio+cross_corr)\weight_corr(:,j,k); % may need to do transform of weight_corr
        end
    end

    for j=1:length(zgrid)
        for k=1:size(zgrid,2)
            anom_value(j,k,time)=weights(:,j,k).'*u_obs; % mean already subtracted
        end
    end
    
end

%% try adding back on the mean, then calculating transport, time mean velocity
% add back on the mean
for i=1:251
    for j=1:182
        for k=1:725
            u_vel(i,j,k)=anom_value_u(i,j,k)+mean_obs_u(1,k);
            v_vel(i,j,k)=anom_value_v(i,j,k)+mean_obs_v(1,k);
        end
    end
end

% make NaN if below topography

b_pos = [27.6428,-33.6674];
d_pos = [27.8603,-34.0435];

B_dx=1000*sw_dist([coast_lat b_pos(2)],[coast_lon b_pos(1)],'km');
D_dx=1000*sw_dist([coast_lat d_pos(2)],[coast_lon d_pos(1)],'km');

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

for i=1:725
    for j=1:182
        for k=1:251
            if depth_int(k,i,j) > abs(topo(j,1))
                u_vel(k,j,i) = NaN;
                v_vel(k,j,i) = NaN; % this deleted unphysical depths, still need to extrapolate inland...
            end
        end
    end
end

% rotate into along, cross stream
[across,along,dist,angle]=rotate_hydro([b_pos(2) d_pos(2)],[b_pos(1) d_pos(1)],u_vel,v_vel);

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