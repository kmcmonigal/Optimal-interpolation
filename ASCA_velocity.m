%% Optimal interpolation for velocity data recovered on ASCA 042018 cruise
% first need to process data: must start and end at same time point, same
% sampling interval, spikes/gaps removed, 40 hour low pass filtered, with u
% and v corresponding to along, across current axis

% first, load our velocity data and organize it so that previous code works
veldir=dir('/Users/kayleenmcmonigal/Documents/Matlab/ASCA0618/mooring_data/ASCA2018_currents_prelim/*.mat');

for i=1:length(veldir)
    load(strcat(veldir(i).folder,'/',veldir(i).name));
end

A.u=A.u_adcp(:,A.sd>=C_adcp.sd(1));
A.v=A.v_adcp(:,A.sd>=C_adcp.sd(1));
A.z=A.z_adcp(:,A.sd>=C_adcp.sd(1));
A.sd2=A.sd(:,A.sd>=C_adcp.sd(1));

B.u=nan(51,13723);
B.v=nan(51,13723);
B.z=nan(51,13723);
tmp.u=B_adcp.u_adcp(:,B_adcp.sd>=C_adcp.sd(1));
tmp2.u=B_rcm.u(B_rcm.sd>=C_adcp.sd(1));
tmp.v=B_adcp.v_adcp(:,B_adcp.sd>=C_adcp.sd(1));
tmp2.v=B_rcm.v(B_rcm.sd>=C_adcp.sd(1));
tmp.z=B_adcp.z_adcp(:,B_adcp.sd>=C_adcp.sd(1));
tmp2.z=B_rcm.z(B_rcm.sd>=C_adcp.sd(1));
for i=1:50
    for j=1:length(tmp.u)
        B.u(i,j)=tmp.u(i,j);
        B.v(i,j)=tmp.v(i,j);
        B.z(i,j)=tmp.z(i,j);
    end
end
count=1;
for i=51
    for j=1:3:length(tmp2.u)-2
        B.u(i,count)=nanmean(tmp2.u(j:j+2));
        B.v(i,count)=nanmean(tmp2.v(j:j+2));
        B.z(i,count)=nanmean(tmp2.z(j:j+2));
        count=count+1;
    end
end

C.u=nan(53,18756);
C.v=nan(53,18756);
C.z=nan(53,18756);
tmp.u=C_adcp.u(:,C_adcp.sd>=C_adcp.sd(1));
tmp2.u=C_rcm.u{1,1}(C_rcm.sd{1,1}>=C_adcp.sd(1));
tmp3.u=C_rcm.u{1,2}(C_rcm.sd{1,2}>=C_adcp.sd(1));
tmp4.u=C_rcm.u{1,3}(C_rcm.sd{1,3}>=C_adcp.sd(1));
tmp.v=C_adcp.v(:,C_adcp.sd>=C_adcp.sd(1));
tmp2.v=C_rcm.v{1,1}(C_rcm.sd{1,1}>=C_adcp.sd(1));
tmp3.v=C_rcm.v{1,2}(C_rcm.sd{1,2}>=C_adcp.sd(1));
tmp4.v=C_rcm.v{1,3}(C_rcm.sd{1,3}>=C_adcp.sd(1));
tmp.z=C_adcp.z(:,C_adcp.sd>=C_adcp.sd(1));
tmp2.z=C_rcm.z{1,1}(C_rcm.sd{1,1}>=C_adcp.sd(1));
tmp3.z=C_rcm.z{1,2}(C_rcm.sd{1,2}>=C_adcp.sd(1));
tmp4.z=C_rcm.z{1,3}(C_rcm.sd{1,3}>=C_adcp.sd(1));
for i=1:50
    for j=1:length(tmp.u)
        C.u(i,j)=tmp.u(i,j);
        C.v(i,j)=tmp.v(i,j);
        C.z(i,j)=tmp.z(i,j);
    end
end
count=1;
for i=51
    for j=1:3:length(tmp2.u)-2
        C.u(i,count)=nanmean(tmp2.u(j:j+2));
        C.v(i,count)=nanmean(tmp2.v(j:j+2));
        C.z(i,count)=nanmean(tmp2.z(j:j+2));
        count=count+1;
    end
end
count=1;
for i=52
    for j=1:3:length(tmp3.u)-2
        C.u(i,count)=nanmean(tmp3.u(j:j+2));
        C.v(i,count)=nanmean(tmp3.v(j:j+2));
        C.z(i,count)=nanmean(tmp3.z(j:j+2));
        count=count+1;
    end
end
count=1;
for i=53
    for j=1:3:length(tmp4.u)-2
        C.u(i,count)=nanmean(tmp4.u(j:j+2));
        C.v(i,count)=nanmean(tmp4.v(j:j+2));
        C.z(i,count)=nanmean(tmp4.z(j:j+2));
        count=count+1;
    end
end

D.u=nan(53,18763);
D.v=nan(53,18763);
D.z=nan(53,18763);
tmp.u=D_adcp.u(:,D_adcp.sd>=C_adcp.sd(1));
tmp2.u=D_rcm.u{1,1}(D_rcm.sd{1,1}>=C_adcp.sd(1));
tmp3.u=D_rcm.u{1,2}(D_rcm.sd{1,2}>=C_adcp.sd(1));
tmp4.u=D_rcm.u{1,3}(D_rcm.sd{1,3}>=C_adcp.sd(1));
tmp.v=D_adcp.v(:,D_adcp.sd>=C_adcp.sd(1));
tmp2.v=D_rcm.v{1,1}(D_rcm.sd{1,1}>=C_adcp.sd(1));
tmp3.v=D_rcm.v{1,2}(D_rcm.sd{1,2}>=C_adcp.sd(1));
tmp4.v=D_rcm.v{1,3}(D_rcm.sd{1,3}>=C_adcp.sd(1));
tmp.z=D_adcp.z(:,D_adcp.sd>=C_adcp.sd(1));
tmp2.z=D_rcm.z{1,1}(D_rcm.sd{1,1}>=C_adcp.sd(1));
tmp3.z=D_rcm.z{1,2}(D_rcm.sd{1,2}>=C_adcp.sd(1));
tmp4.z=D_rcm.z{1,3}(D_rcm.sd{1,3}>=C_adcp.sd(1));
for i=1:50
    for j=1:length(tmp.u)
        D.u(i,j)=tmp.u(i,j);
        D.v(i,j)=tmp.v(i,j);
        D.z(i,j)=tmp.z(i,j);
    end
end
count=1;
for i=51
    for j=1:3:length(tmp2.u)-2
        D.u(i,count)=nanmean(tmp2.u(j:j+2));
        D.v(i,count)=nanmean(tmp2.v(j:j+2));
        D.z(i,count)=nanmean(tmp2.z(j:j+2));
        count=count+1;
    end
end
count=1;
for i=52
    for j=1:3:length(tmp3.u)-2
        D.u(i,count)=nanmean(tmp3.u(j:j+2));
        D.v(i,count)=nanmean(tmp3.v(j:j+2));
        D.z(i,count)=nanmean(tmp3.z(j:j+2));
        count=count+1;
    end
end
count=1;
for i=53
    for j=1:3:length(tmp4.u)-2
        D.u(i,count)=nanmean(tmp4.u(j:j+2));
        D.v(i,count)=nanmean(tmp4.v(j:j+2));
        D.z(i,count)=nanmean(tmp4.z(j:j+2));
        count=count+1;
    end
end

E.u=nan(54,56343/3);
E.v=nan(54,56343/3);
E.z=nan(54,56343/3);
tmp.u=E_adcp.u(:,E_adcp.sd>=C_adcp.sd(1));
tmp2.u=E_nrtk.u{1,1}(E_nrtk.sd{1,1}>=C_adcp.sd(1));
tmp3.u=E_nrtk.u{1,2}(E_nrtk.sd{1,2}>=C_adcp.sd(1));
tmp4.u=E_nrtk.u{1,3}(E_nrtk.sd{1,3}>=C_adcp.sd(1));
tmp5.u=E_nrtk.u{1,4}(E_nrtk.sd{1,4}>=C_adcp.sd(1));
tmp.v=E_adcp.v(:,E_adcp.sd>=C_adcp.sd(1));
tmp2.v=E_nrtk.v{1,1}(E_nrtk.sd{1,1}>=C_adcp.sd(1));
tmp3.v=E_nrtk.v{1,2}(E_nrtk.sd{1,2}>=C_adcp.sd(1));
tmp4.v=E_nrtk.v{1,3}(E_nrtk.sd{1,3}>=C_adcp.sd(1));
tmp5.v=E_nrtk.v{1,4}(E_nrtk.sd{1,4}>=C_adcp.sd(1));
tmp.z=E_adcp.z(:,E_adcp.sd>=C_adcp.sd(1));
tmp2.z=E_nrtk.z{1,1}(E_nrtk.sd{1,1}>=C_adcp.sd(1));
tmp3.z=E_nrtk.z{1,2}(E_nrtk.sd{1,2}>=C_adcp.sd(1));
tmp4.z=E_nrtk.z{1,3}(E_nrtk.sd{1,3}>=C_adcp.sd(1));
tmp5.z=E_nrtk.z{1,4}(E_nrtk.sd{1,4}>=C_adcp.sd(1));
for i=1:50
    for j=1:length(tmp.u)
        E.u(i,j)=tmp.u(i,j);
        E.v(i,j)=tmp.v(i,j);
        E.z(i,j)=tmp.z(i,j);
    end
end
count=1;
for i=51
    for j=1:3:length(tmp2.u)-2
        E.u(i,count)=nanmean(tmp2.u(j:j+2));
        E.v(i,count)=nanmean(tmp2.v(j:j+2));
        E.z(i,count)=nanmean(tmp2.z(j:j+2));
        count=count+1;
    end
end
count=1;
for i=52
    for j=1:3:length(tmp3.u)-2
        E.u(i,count)=nanmean(tmp3.u(j:j+2));
        E.v(i,count)=nanmean(tmp3.v(j:j+2));
        E.z(i,count)=nanmean(tmp3.z(j:j+2));
        count=count+1;
    end
end
count=1;
for i=53
    for j=1:3:length(tmp4.u)-2
        E.u(i,count)=nanmean(tmp4.u(j:j+2));
        E.v(i,count)=nanmean(tmp4.v(j:j+2));
        E.z(i,count)=nanmean(tmp4.z(j:j+2));
        count=count+1;
    end
end
count=1;
for i=54
    for j=1:3:length(tmp5.u)-2
        E.u(i,count)=nanmean(tmp5.u(j:j+2));
        E.v(i,count)=nanmean(tmp5.v(j:j+2));
        E.z(i,count)=nanmean(tmp5.z(j:j+2));
        count=count+1;
    end
end

F.u=nan(55,55980/3);
F.v=nan(55,55980/3);
F.z=nan(55,55980/3);
tmp.u=F_adcp.u(:,F_adcp.sd>=C_adcp.sd(1));
tmp2.u=F_nrtk.u{1,1}(F_nrtk.sd{1,1}>=C_adcp.sd(1));
tmp3.u=F_nrtk.u{1,2}(F_nrtk.sd{1,2}>=C_adcp.sd(1));
tmp4.u=F_nrtk.u{1,3}(F_nrtk.sd{1,3}>=C_adcp.sd(1));
tmp5.u=F_nrtk.u{1,4}(F_nrtk.sd{1,4}>=C_adcp.sd(1));
tmp6.u=F_nrtk.u{1,5}(F_nrtk.sd{1,5}>=C_adcp.sd(1));
tmp.v=F_adcp.v(:,F_adcp.sd>=C_adcp.sd(1));
tmp2.v=F_nrtk.v{1,1}(F_nrtk.sd{1,1}>=C_adcp.sd(1));
tmp3.v=F_nrtk.v{1,2}(F_nrtk.sd{1,2}>=C_adcp.sd(1));
tmp4.v=F_nrtk.v{1,3}(F_nrtk.sd{1,3}>=C_adcp.sd(1));
tmp5.v=F_nrtk.v{1,4}(F_nrtk.sd{1,4}>=C_adcp.sd(1));
tmp6.v=F_nrtk.v{1,5}(F_nrtk.sd{1,5}>=C_adcp.sd(1));
tmp.z=F_adcp.z(:,F_adcp.sd>=C_adcp.sd(1));
tmp2.z=F_nrtk.z{1,1}(F_nrtk.sd{1,1}>=C_adcp.sd(1));
tmp3.z=F_nrtk.z{1,2}(F_nrtk.sd{1,2}>=C_adcp.sd(1));
tmp4.z=F_nrtk.z{1,3}(F_nrtk.sd{1,3}>=C_adcp.sd(1));
tmp5.z=F_nrtk.z{1,4}(F_nrtk.sd{1,4}>=C_adcp.sd(1));
tmp6.z=F_nrtk.z{1,5}(F_nrtk.sd{1,5}>=C_adcp.sd(1));
for i=1:50
    for j=1:length(tmp.u)
        F.u(i,j)=tmp.u(i,j);
        F.v(i,j)=tmp.v(i,j);
        F.z(i,j)=tmp.z(i,j);
    end
end
count=1;
for i=51
    for j=1:3:length(tmp2.u)-2
        F.u(i,count)=nanmean(tmp2.u(j:j+2));
        F.v(i,count)=nanmean(tmp2.v(j:j+2));
        F.z(i,count)=nanmean(tmp2.z(j:j+2));
        count=count+1;
    end
end
count=1;
for i=52
    for j=1:3:length(tmp3.u)-2
        F.u(i,count)=nanmean(tmp3.u(j:j+2));
        F.v(i,count)=nanmean(tmp3.v(j:j+2));
        F.z(i,count)=nanmean(tmp3.z(j:j+2));
        count=count+1;
    end
end
count=1;
for i=53
    for j=1:3:length(tmp4.u)-2
        F.u(i,count)=nanmean(tmp4.u(j:j+2));
        F.v(i,count)=nanmean(tmp4.v(j:j+2));
        F.z(i,count)=nanmean(tmp4.z(j:j+2));
        count=count+1;
    end
end
count=1;
for i=54
    for j=1:3:length(tmp5.u)-2
        F.u(i,count)=nanmean(tmp5.u(j:j+2));
        F.v(i,count)=nanmean(tmp5.v(j:j+2));
        F.z(i,count)=nanmean(tmp5.z(j:j+2));
        count=count+1;
    end
end
count=1;
for i=55
    for j=1:3:length(tmp6.u)-2
        F.u(i,count)=nanmean(tmp6.u(j:j+2));
        F.v(i,count)=nanmean(tmp6.v(j:j+2));
        F.z(i,count)=nanmean(tmp6.z(j:j+2));
        count=count+1;
    end
end


G.u=nan(55,55980/3);
G.v=nan(55,55980/3);
G.z=nan(55,55980/3);
tmp.u=G_adcp.u(:,G_adcp.sd>=C_adcp.sd(1));
tmp2.u=G_nrtk.u{1,1}(G_nrtk.sd{1,1}>=C_adcp.sd(1));
tmp3.u=G_nrtk.u{1,2}(G_nrtk.sd{1,2}>=C_adcp.sd(1));
tmp4.u=G_nrtk.u{1,3}(G_nrtk.sd{1,3}>=C_adcp.sd(1));
tmp5.u=G_nrtk.u{1,4}(G_nrtk.sd{1,4}>=C_adcp.sd(1));
tmp6.u=G_nrtk.u{1,5}(G_nrtk.sd{1,5}>=C_adcp.sd(1));
tmp.v=G_adcp.v(:,G_adcp.sd>=C_adcp.sd(1));
tmp2.v=G_nrtk.v{1,1}(G_nrtk.sd{1,1}>=C_adcp.sd(1));
tmp3.v=G_nrtk.v{1,2}(G_nrtk.sd{1,2}>=C_adcp.sd(1));
tmp4.v=G_nrtk.v{1,3}(G_nrtk.sd{1,3}>=C_adcp.sd(1));
tmp5.v=G_nrtk.v{1,4}(G_nrtk.sd{1,4}>=C_adcp.sd(1));
tmp6.v=G_nrtk.v{1,5}(G_nrtk.sd{1,5}>=C_adcp.sd(1));
tmp.z=G_adcp.z(:,G_adcp.sd>=C_adcp.sd(1));
tmp2.z=G_nrtk.z{1,1}(G_nrtk.sd{1,1}>=C_adcp.sd(1));
tmp3.z=G_nrtk.z{1,2}(G_nrtk.sd{1,2}>=C_adcp.sd(1));
tmp4.z=G_nrtk.z{1,3}(G_nrtk.sd{1,3}>=C_adcp.sd(1));
tmp5.z=G_nrtk.z{1,4}(G_nrtk.sd{1,4}>=C_adcp.sd(1));
tmp6.z=G_nrtk.z{1,5}(G_nrtk.sd{1,5}>=C_adcp.sd(1));
for i=1:50
    for j=1:length(tmp.u)
        G.u(i,j)=tmp.u(i,j);
        G.v(i,j)=tmp.v(i,j);
        G.z(i,j)=tmp.z(i,j);
    end
end
count=1;
for i=51
    for j=1:3:length(tmp2.u)-2
        G.u(i,count)=nanmean(tmp2.u(j:j+2));
        G.v(i,count)=nanmean(tmp2.v(j:j+2));
        G.z(i,count)=nanmean(tmp2.z(j:j+2));
        count=count+1;
    end
end
count=1;
for i=52
    for j=1:3:length(tmp3.u)-2
        G.u(i,count)=nanmean(tmp3.u(j:j+2));
        G.v(i,count)=nanmean(tmp3.v(j:j+2));
        G.z(i,count)=nanmean(tmp3.z(j:j+2));
        count=count+1;
    end
end
count=1;
for i=53
    for j=1:3:length(tmp4.u)-2
        G.u(i,count)=nanmean(tmp4.u(j:j+2));
        G.v(i,count)=nanmean(tmp4.v(j:j+2));
        G.z(i,count)=nanmean(tmp4.z(j:j+2));
        count=count+1;
    end
end
count=1;
for i=54
    for j=1:3:length(tmp5.u)-2
        G.u(i,count)=nanmean(tmp5.u(j:j+2));
        G.v(i,count)=nanmean(tmp5.v(j:j+2));
        G.z(i,count)=nanmean(tmp5.z(j:j+2));
        count=count+1;
    end
end
count=1;
for i=55
    for j=1:3:length(tmp6.u)-2
        G.u(i,count)=nanmean(tmp6.u(j:j+2));
        G.v(i,count)=nanmean(tmp6.v(j:j+2));
        G.z(i,count)=nanmean(tmp6.z(j:j+2));
        count=count+1;
    end
end

% across=-A.u_adcp(30,:)*sind(64)-A.v_adcp(30,:)*cosd(64); % seems correct to get cross array velocity
% along=A.u_adcp(30,:)*cosd(64)-A.v_adcp(30,:)*sind(64);

% 40 hour low pass filter records,
% 20 hour subsample
d1 = designfilt('lowpassiir','FilterOrder',12, ...
    'HalfPowerFrequency',1/20,'DesignMethod','butter');

A.gaps=isnan(A.u);
B.gaps=isnan(B.u);
C.gaps=isnan(C.u);
D.gaps=isnan(D.u);
E.gaps=isnan(E.u);
F.gaps=isnan(F.u);
G.gaps=isnan(G.u);

A.u(isnan(A.u))=0;A.v(isnan(A.v))=0;A.z(isnan(A.z))=0;
B.u(isnan(B.u))=0;B.v(isnan(B.v))=0;B.z(isnan(B.z))=0;
C.u(isnan(C.u))=0;C.v(isnan(C.v))=0;C.z(isnan(C.z))=0;
D.u(isnan(D.u))=0;D.v(isnan(D.v))=0;D.z(isnan(D.z))=0;
E.u(isnan(E.u))=0;E.v(isnan(E.v))=0;E.z(isnan(E.z))=0;
F.u(isnan(F.u))=0;F.v(isnan(F.v))=0;F.z(isnan(F.z))=0;
G.u(isnan(G.u))=0;G.v(isnan(G.v))=0;G.z(isnan(G.z))=0;

for i=1:size(A.u,1)
    A.u2(i,:)=filtfilt(d1,A.u(i,:)); 
    A.v2(i,:)=filtfilt(d1,A.v(i,:));
    A.z2(i,:)=filtfilt(d1,A.z(i,:));
end
for i=1:size(B.u,1)
    B.u2(i,:)=filtfilt(d1,B.u(i,:));
    B.v2(i,:)=filtfilt(d1,B.v(i,:));
    B.z2(i,:)=filtfilt(d1,B.z(i,:));
end
for i=1:size(C.u,1)
    C.u2(i,:)=filtfilt(d1,C.u(i,:));
    C.v2(i,:)=filtfilt(d1,C.v(i,:));
    C.z2(i,:)=filtfilt(d1,C.z(i,:));
end
for i=1:size(D.u,1)
    D.u2(i,:)=filtfilt(d1,D.u(i,:));
    D.v2(i,:)=filtfilt(d1,D.v(i,:));
    D.z2(i,:)=filtfilt(d1,D.z(i,:));
end
for i=1:size(E.u,1)
    E.u2(i,:)=filtfilt(d1,E.u(i,:));
    E.v2(i,:)=filtfilt(d1,E.v(i,:));
    E.z2(i,:)=filtfilt(d1,E.z(i,:));
end
for i=1:size(F.u,1)
    F.u2(i,:)=filtfilt(d1,F.u(i,:));
    F.v2(i,:)=filtfilt(d1,F.v(i,:));
    F.z2(i,:)=filtfilt(d1,F.z(i,:));
end
for i=1:size(G.u,1)
    G.u2(i,:)=filtfilt(d1,G.u(i,:));
    G.v2(i,:)=filtfilt(d1,G.v(i,:));
    G.z2(i,:)=filtfilt(d1,G.z(i,:));
end
A.u2(A.gaps)=nan;A.v2(A.gaps)=nan;A.z(A.gaps)=nan;
B.u2(B.gaps)=nan;B.v2(B.gaps)=nan;B.z(B.gaps)=nan;
C.u2(C.gaps)=nan;C.v2(C.gaps)=nan;C.z(C.gaps)=nan;
D.u2(D.gaps)=nan;D.v2(D.gaps)=nan;D.z(D.gaps)=nan;
E.u2(E.gaps)=nan;E.v2(E.gaps)=nan;E.z(E.gaps)=nan;
F.u2(F.gaps)=nan;F.v2(F.gaps)=nan;F.z(F.gaps)=nan;
G.u2(G.gaps)=nan;G.v2(G.gaps)=nan;G.z(G.gaps)=nan;
A.u3=A.u2(:,1:20:end);A.v3=A.v2(:,1:20:end);A.z3=A.z2(:,1:20:end);
B.u3=B.u2(:,1:20:end);B.v3=B.v2(:,1:20:end);B.z3=B.z2(:,1:20:end);
C.u3=C.u2(:,1:20:end);C.v3=C.v2(:,1:20:end);C.z3=C.z2(:,1:20:end);
D.u3=D.u2(:,1:20:end);D.v3=D.v2(:,1:20:end);D.z3=D.z2(:,1:20:end);
E.u3=E.u2(:,1:20:end);E.v3=E.v2(:,1:20:end);E.z3=E.z2(:,1:20:end);
F.u3=F.u2(:,1:20:end);F.v3=F.v2(:,1:20:end);F.z3=F.z2(:,1:20:end);
G.u3=G.u2(:,1:20:end);G.v3=G.v2(:,1:20:end);G.z3=G.z2(:,1:20:end);

% make them all as long as the longest record
A.u=nan(size(A.u3,1),940);A.v=nan(size(A.u3,1),940);A.z=nan(size(A.u3,1),940);
B.u=nan(size(B.u3,1),940);B.v=nan(size(B.u3,1),940);B.z=nan(size(B.u3,1),940);
C.u=nan(size(C.u3,1),940);C.v=nan(size(C.u3,1),940);C.z=nan(size(C.u3,1),940);
D.u=nan(size(D.u3,1),940);D.v=nan(size(D.u3,1),940);D.z=nan(size(D.u3,1),940);
E.u=nan(size(E.u3,1),940);E.v=nan(size(E.u3,1),940);E.z=nan(size(E.u3,1),940);
F.u=nan(size(F.u3,1),940);F.v=nan(size(F.u3,1),940);F.z=nan(size(F.u3,1),940);
G.u=nan(size(G.u3,1),940);G.v=nan(size(G.u3,1),940);G.z=nan(size(G.u3,1),940);

for i=1:length(A.u3)
    A.u(:,i)=A.u3(:,i);
    A.v(:,i)=A.v3(:,i);
    A.z(:,i)=A.z3(:,i);
end

for i=1:length(B.u3)
    B.u(:,i)=B.u3(:,i);
    B.v(:,i)=B.v3(:,i);
    B.z(:,i)=B.z3(:,i);
end

for i=1:length(C.u3)
    C.u(:,i)=C.u3(:,i);
    C.v(:,i)=C.v3(:,i);
    C.z(:,i)=C.z3(:,i);
end

for i=1:length(D.u3)
    D.u(:,i)=D.u3(:,i);
    D.v(:,i)=D.v3(:,i);
    D.z(:,i)=D.z3(:,i);
end

for i=1:length(E.u3)
    E.u(:,i)=E.u3(:,i);
    E.v(:,i)=E.v3(:,i);
    E.z(:,i)=E.z3(:,i);
end

for i=1:length(F.u3)
    F.u(:,i)=F.u3(:,i);
    F.v(:,i)=F.v3(:,i);
    F.z(:,i)=F.z3(:,i);
end

for i=1:length(G.u3)
    G.u(:,i)=G.u3(:,i);
    G.v(:,i)=G.v3(:,i);
    G.z(:,i)=G.z3(:,i);
end

% delete things we don't need to keep it cleaner
clear B_adcp B_rcm C_adcp C_rcm D_adcp D_rcm E_adcp E_nrtk F_adcp F_nrtk G_adcp G_nrtk tmp tmp2 tmp3 tmp4 tmp5 tmp6 A.u2 A.u3 A.v2 A.v3 A.z2 A.z3
clear B.u2 B.v2 B.z2 B.u3 B.v3 B.z3 C.u2 C.u3 C.v2 C.v3 C.z2 C.z3 D.u2 D.u3 D.v2 D.v3 D.z2 D.z3 E.u2 E.v2 E.z2 E.u3 E.v3 E.z3 F.u2 F.u3 F.v2 F.v3 F.z2
clear F.z3 G.u2 G.u3 G.v3 G.v2 G.z2 G.z3 count

%% mapping begins
% make some choices
use_background=0; % subtract a "background field" and map anomalies vs mapping field itself
Lx=63*1000; % horizontal decorrelation length scale
Lz=1913; % vertical decorrelation length scale
lx=.2*63*1000; % small scale decorrelation length scales 
lz=.2*1913;
N=.07; % noise ratio

% correlation functions
Xc=@(x) exp(-(x(:)/Lx).^2).*cos(pi.*x(:)./(2.*Lx))+exp(-(x(:)/lx).^2).*cos(pi.*x(:)./(2.*lx));
Zc=@(z) exp(-(z(:)/Lz).^2)+exp(-(z(:)/lz).^2);

% instrument locations - check that these are the locations we want to use
a_pos = [27.595,-33.5583];
b_pos = [27.6428,-33.6674];
c_pos = [27.7152,-33.7996];
d_pos = [27.8603,-34.0435];
e_pos = [28.0316,-34.29];
f_pos = [28.17,-34.5380];
g_pos = [28.3453,-34.8213];

% map with position of mooring a as left bound, mooring g as right bound
B_dx=1000*sw_dist([a_pos(2) b_pos(2)],[a_pos(1) b_pos(1)],'km');
C_dx=1000*sw_dist([a_pos(2) c_pos(2)],[a_pos(1) c_pos(1)],'km');
D_dx=1000*sw_dist([a_pos(2) d_pos(2)],[a_pos(1) d_pos(1)],'km');
E_dx=1000*sw_dist([a_pos(2) e_pos(2)],[a_pos(1) e_pos(1)],'km');
F_dx=1000*sw_dist([a_pos(2) f_pos(2)],[a_pos(1) f_pos(1)],'km');
G_dx=1000*sw_dist([a_pos(2) g_pos(2)],[a_pos(1) g_pos(1)],'km');

% make a matrix of the distance of each observation from A, will need to
% fill in
dist=zeros(1,size(A.u,1)+size(B.u,1)+size(C.u,1)+size(D.u,1)+size(E.u,1)+size(F.u,1)+size(G.u,1));
for i=1:size(A.u,1)
    dist(i)=0;
end
for i=size(A.u,1)+1:size(A.u,1)+size(B.u,1)
    dist(i)=B_dx;
end
for i=size(A.u,1)+size(B.u,1)+1:size(A.u,1)+size(B.u,1)+size(C.u,1)
    dist(i)=C_dx;
end
for i=size(A.u,1)+size(B.u,1)+size(C.u,1)+1:size(A.u,1)+size(B.u,1)+size(C.u,1)+size(D.u,1)
    dist(i)=D_dx;
end
for i=size(A.u,1)+size(B.u,1)+size(C.u,1)+size(D.u,1)+1:size(A.u,1)+size(B.u,1)+size(C.u,1)+size(D.u,1)+size(E.u,1)
    dist(i)=E_dx;
end
for i=size(A.u,1)+size(B.u,1)+size(C.u,1)+size(D.u,1)+size(E.u,1)+1:size(A.u,1)+size(B.u,1)+size(C.u,1)+size(D.u,1)+size(E.u,1)+size(F.u,1)
    dist(i)=F_dx;
end
for i=size(A.u,1)+size(B.u,1)+size(C.u,1)+size(D.u,1)+size(E.u,1)+size(F.u,1)+1:size(A.u,1)+size(B.u,1)+size(C.u,1)+size(D.u,1)+size(E.u,1)+size(F.u,1)+size(G.u,1)
    dist(i)=G_dx;
end

% make a grid to map data onto
x=0:500:G_dx;
z=0:20:4500; % then later will make anything under topography nan

[xgrid,zgrid]=meshgrid(x,z);

% set up a noise matrix by finding variace of each measurement. this is the matrix B in 4.2.10 of Daley ("background error covariance matrix")

noisev=N./nanvar([A.v;B.v;C.v;D.v;E.v;F.v;G.v].');
noiseu=N./nanvar([A.u;B.u;C.u;D.u;E.u;F.v;G.u].'); % noise:signal ratio is instrument noise (set above)/variance of measurements
noise=(noiseu+noisev)/2; % take the mean since they are very similar

z_mean=[nanmean(A.z,2);nanmean(B.z,2);nanmean(C.z,2);nanmean(D.z,2);nanmean(E.z,2);nanmean(F.z,2);nanmean(G.z,2)];

u=[nanmean(A.u,2);nanmean(B.u,2);nanmean(C.u,2);nanmean(D.u,2);nanmean(E.u,2);nanmean(F.u,2);nanmean(G.u,2)]; 
v=[nanmean(A.v,2);nanmean(B.v,2);nanmean(C.v,2);nanmean(D.v,2);nanmean(E.v,2);nanmean(F.v,2);nanmean(G.v,2)];

mean_u=0;%nanmean(nanmean(u)); % we could assume some background field but let's not for now
mean_v=0;%nanmean(nanmean(v));

u_anom=u-mean_u;
v_anom=v-mean_v;
gaps=isnan(u_anom);
u_anom2=u_anom(gaps==0);
v_anom2=v_anom(gaps==0);
z_mean2=z_mean(gaps==0);
noise2=noise(gaps==0);
dist2=dist(gaps==0);
noise3=diag(noise2);

% this creates a matrix which is the error covariance between each instrument and each grid point (B_i in 4.2.10 of Daley)
weight_corr=nan(length(noise2),size(zgrid,1),size(zgrid,2));
for j=1:length(noise2)
    for k=1:size(zgrid,1)
        for l=1:size(zgrid,2)
            weight_corr(j,k,l)=Xc(abs(xgrid(k,l)-dist2(j)))*Zc(abs(zgrid(k,l)-z_mean2(j)));
        end
    end
end

% get cross correlation between instruments. this is the matrix O in 4.2.10 of Daley ("observation error covariance matrix")
cross_corr=nan(length(noise2),length(noise2));
for j=1:length(noise2)
    for k=1:length(noise2)
        cross_corr(j,k)=Xc(abs(dist2(j)-dist2(k)))*Zc(abs(z_mean2(j)-z_mean2(k)));
    end
end

% solve for weights (W_i in Daley)
weights=nan(size(weight_corr,1),size(zgrid,1),size(zgrid,2));
for j=1:size(zgrid,1)
    for k=1:size(zgrid,2)
        weights(:,j,k)=((noise3+cross_corr)\weight_corr(:,j,k));
    end
end

% use the weights that we just solved for to get the u and v velocity at each grid point
for j=1:size(zgrid,1)
    for k=1:size(zgrid,2)
        vel_mean_u(j,k)=(weights(:,j,k).'*u_anom2)+mean_u;
        vel_mean_v(j,k)=(weights(:,j,k).'*v_anom2)+mean_v;
    end
end

%% now map anomalies relative to that time mean
% full weights
tic
for i=1:940 % map over each time step
    % put all the obervations from that point in time into a matrix
    u=[A.u(:,i);B.u(:,i);C.u(:,i);D.u(:,i);E.u(:,i);F.u(:,i);G.u(:,i)]; 
    v=[A.v(:,i);B.v(:,i);C.v(:,i);D.v(:,i);E.v(:,i);F.v(:,i);G.v(:,i)];% cpies34(:),i);cpies45(:,i)];
    
    % depth of each of those observations (no CPIES for now)
    z=[A.z(:,i);B.z(:,i);C.z(:,i);D.z(:,i);E.z(:,i);F.z(:,i);G.z(:,i)];% dpth34_(~isnan(cpies34_(:,i)));dpth45_(~isnan(cpies45_(:,i)))];
    
    % see if any of the observations are nan, get rid of those
    gaps=[isnan(A.u(:,i));isnan(B.u(:,i));isnan(C.u(:,i));isnan(D.u(:,i));isnan(E.u(:,i));isnan(F.u(:,i));isnan(G.u(:,i))];
    
    u2=u(gaps==0); % see if more efficient way to do this
    v2=v(gaps==0);
    z2=z(gaps==0);
    dist2=dist(gaps==0);
    noisev2=noisev(gaps==0);
    noiseu2=noiseu(gaps==0);

    % make diagonal noise matrix
    noise=diag((noiseu2+noisev2)/2); % just take mean of u, v noise since they are very similar
    
    vel_mean_anomu=nan(1,size(z2,1));
    vel_mean_anomv=nan(1,size(z2,1));
    for j=1:size(z2,1)
        for k=1:size(zgrid,1)
            dz_dist(j,k)=abs(zgrid(k,1)-z2(j));
        end
        [Mz,Iz] = min(dz_dist(j,:));
        for k=1:size(xgrid,2)
            grid_dist(j,k)=abs(xgrid(1,k)-dist2(j));
        end
        [M,I] = min(grid_dist(j,:));
        vel_mean_anomu(j)=u2(j)-vel_mean_u(Iz,I);
        vel_mean_anomv(j)=v2(j)-vel_mean_v(Iz,I);
    end
   
    
    % now get cross correlations between instruments, grid points
    weight_corr2=nan(length(noise),size(zgrid,1),size(zgrid,2));
    for j=1:length(noise)
        for k=1:size(zgrid,1)
            for l=1:size(zgrid,2)
                weight_corr2(j,k,l)=Xc(abs(xgrid(k,l)-dist2(j)))*Zc(abs(zgrid(k,l)-z2(j)));
            end
        end
    end
    
    % get cross corr between instruments
    cross_corr2=nan(length(noise),length(noise));
    for j=1:length(noise)
        for k=1:length(noise)
            cross_corr2(j,k)=Xc(abs(dist2(j)-dist2(k)))*Zc(abs(z2(j)-z2(k)));
        end
    end
    
    % solve for weights
    weights2=nan(size(weight_corr2,1),size(zgrid,1),size(zgrid,2));
    for j=1:size(zgrid,1)
        for k=1:size(zgrid,2)
            weights2(:,j,k)=((noise+cross_corr2)\weight_corr2(:,j,k)); 
            temp_u(j,k)=(weights2(:,j,k).'*vel_mean_anomu.');
            temp_v(j,k)=(weights2(:,j,k).'*vel_mean_anomv.');
        end
    end
    vel_total_u(:,:,i)=temp_u+vel_mean_u;%+vel_mean2u;
    vel_total_v(:,:,i)=temp_v+vel_mean_v;%+vel_mean2v;
    
    % repeat time varying on small scale - don't do this for now as it doubles the time to run

%     for j=1:size(dz_obs,1)
%         for k=1:size(zgrid,1)
%             dz_dist(j,k)=abs(zgrid(k,1)-dz_obs(j));
%         end
%         [Mz,Iz] = min(dz_dist(j,:));
%         for k=1:size(xgrid,2)
%             grid_dist(j,k)=abs(xgrid(1,k)-dx(j));
%         end
%         [M,I] = min(grid_dist(j,:));
%         temp=tem(Iz,I);
%         ac_obs_anom(j)=ac_obs(j)-temp;
%     end
%     
%     cross_corr=nan(length(Noise2),length(Noise2));
%     for j=1:length(Noise2)
%         for k=1:length(Noise2)
%             cross_corr(j,k)=Xc(abs(dx(j)-dx(k)))*Zc(abs(dz_obs(j)-dz_obs(k)));
%         end
%     end
%     
%     % weights are different than before
%     
%     weights_ac=nan(size(weight_corr,1),size(zgrid,1),size(zgrid,2));
%     for j=1:size(zgrid,1)
%         for k=1:size(zgrid,2)
%             weights_ac(:,j,k)=(ratio_ac+cross_corr)\weight_corr(:,j,k);
%         end
%     end
%     
%     for j=1:size(zgrid,1)
%         for k=1:size(zgrid,2)
%             temp2(j,k)=weights_ac(:,j,k).'*ac_obs_anom.';
%         end
%     end
%     int_ac2(:,:,time)=temp2;

    disp([i toc])
end

save('vel_mapping','vel_total_u','vel_total_v','xgrid','zgrid','Lx','lx','Lz','lz')
