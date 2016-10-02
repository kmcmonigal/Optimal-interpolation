% optimally interpolate salinity
% first create a climatology to subtract 
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