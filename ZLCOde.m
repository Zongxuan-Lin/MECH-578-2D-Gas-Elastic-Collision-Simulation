clc;
clear all;
%% given parameters and table setup
t=0; %time start at 0
ancol=0;
boxsize=21; %box size
num=(boxsize-1)^2;%A total of 400 molecules in the simulation
moleculed=0.1; %diameter and mass of molecule
molecuelm=1;
m2=molecuelm*2;
vo=1; %initial velocity for all molecules;
Etot=0.5*molecuelm*vo^2;
bconst=1.0; %fictional boltzman constant in the problem
temp=1/2;%initial temperature
ncol=200000;%number of collisions
volrec=zeros(ncol+1,num,2);
posrec=zeros(ncol+1,num,2);
spdrec=zeros(ncol+1,num,1);
inidir=rand(1,num)*2*pi;%initialized direction of travel
pos=zeros(num,2);%position table of the particles
vel=zeros(num,2);%velocity table of the particles
speed=zeros(num,1);
ancols=zeros(ncol+1,1);
time=zeros(ncol+1,1);
dts=zeros(ncol,1);
Cnum=0;
afb=(num+1)*ones(num,1);%affected ball id in collision
%% initial Position and velocity Setup
for i=1:num
    speed(i)=1;
    x=fix(i/(boxsize-1))+1;
    y=mod(i,20);
    if y==0
        y=20; x=x-1;
    end
    pos(i,1)=x;pos(i,2)=y;%initialize position
    vel(i,1)=vo*cos(inidir(i));vel(i,2)=vo*sin(inidir(i));%initialize velocity
end
volrec(1,:,:)=vel; posrec(1,:,:)=pos; spdrec(1,:,:)=speed;
%% First collision time calculation
comb=nchoosek(num,2);%total number of pairs to consider
colliBall=zeros(comb,3);%potential collision between particles table, update between every collision
colliWall=zeros(num,3);%potential collision between particles and wall table, update between every collision
temppos=1; %navigator for collision table

%first run of collision between ball simulation
for i=1:(num-1)
    for j=(i+1):num
        colliBall(temppos,1)=i;%initializing pair id
        colliBall(temppos,2)=j;
        rel_pos=pos(i,:)-pos(j,:);
        rel_vel=vel(i,:)-vel(j,:);
        mark1=dot(rel_pos, rel_vel);%Dot Product of rel position and rel velocity
        velsq=norm(rel_vel)^2;
        possq=norm(rel_pos)^2;
        mark2=mark1^2-velsq*(possq-moleculed^2);
        if mark1<0 && mark2>=0
            th1=(-mark1+sqrt(mark2))/velsq; th2=(-mark1-sqrt(mark2))/velsq;
            colliBall(temppos,3)=compare(th1,th2);
        else
            colliBall(temppos,3)=inf;
        end
        temppos=temppos+1;
    end
end
%wall 1 left  0
%wall 2 right     21
%wall 3 bottom    0
%wall 4 top   21
for i=1:num
    colliWall(i,1)=i;
    posx=pos(i,1);  posy=pos(i,2);  velx=vel(i,1);  vely=vel(i,2);
    tx=0;ty=0;wx=0;wy=0;
    if (velx>0)
        tx=(boxsize-posx-moleculed/2)/velx;wx=2;
    elseif(velx<0)
        tx=(0-posx+moleculed/2)/velx;wx=1;
    else
        tx=inf;wx=5;
    end
    if (vely>0)
        ty=(boxsize-posy-moleculed/2)/vely;wy=4;
    elseif(vely<0)
        ty=(0-posy+moleculed/2)/vely;wy=3;
    else
        ty=inf;wy=5;
    end
    if(tx<0)
        tx=0;
    end
    if(ty<0)
        ty=0;
    end
    if (tx<ty)
        colliWall(i,2)=wx;colliWall(i,3)=tx;
    elseif(tx>ty)
        colliWall(i,2)=wy; colliWall(i,3)=ty;
    elseif(tx~=inf && tx==ty)
        colliWall(i,2)=6;colliWall(i,3)=ty;%code 6 collision with both wall together
    else
        colliWall(i,2)=5;colliWall(i,3)=inf;
    end
end

%% first post-collision Velocity Analysis 
%find the first collision happend
%2d elastic collision between balls calc credit to fiflulu
   temporarystorage1=colliBall(:,3);temporarystorage2=colliWall(:,3);
   [mint_ball, minBallInd] = min(temporarystorage1);
   [mint_wall, minwallInd] = min(temporarystorage2);
   ancol=ancol+1;
%% Collision with Wall Velocity Analysis 
if mint_wall<mint_ball%hit wall first
    wallid=colliWall(minwallInd,2);deltat=mint_wall;
    t=t+deltat;
    dts(1)=deltat;
    time(2)=t;
    pos=pos+vel*deltat;
    if wallid <=2%x direction
        vel(minwallInd,1)=-vel(minwallInd,1);
    elseif wallid>2 && wallid <=4%y direction
        vel(minwallInd,2)=-vel(minwallInd,2);
    elseif wallid==6%in the corner
        vel(minwallInd,1)=-vel(minwallInd,1);
        vel(minwallInd,2)=-vel(minwallInd,2);
    end
    afb(1)=minwallInd;Cnum=1;
    for i=1:comb
        xx=colliBall(i,1);yx=colliBall(i,2);
        dist=norm(pos(xx,:)-pos(yx,:));
        if dist<moleculed
            Cnum=Cnum+2;
            afb(Cnum-1)=xx;afb(Cnum)=yx;
            r12unit=(pos(yx,:)-pos(xx,:))/(norm(pos(yx,:)-pos(xx,:)));
            cmnterm=dot((vel(xx,:)-vel(yx,:)),r12unit)*r12unit;
            vb1=vel(xx,:)-cmnterm;vel(xx,:)=vb1;
            vb2=vel(yx,:)+cmnterm;vel(yx,:)=vb2;
            speed(xx)=norm(vel(xx,:));speed(yx)=norm(vel(yx,:));
            posxx=pos(xx,:)+vb1*(moleculed-dist)/speed(xx);pos(xx,:)=posxx;
            posyx=pos(yx,:)+vb2*(moleculed-dist)/speed(yx);pos(yx,:)=posyx;
            ancol=ancol+1;
        end
    end
    volrec(2,:,:)=vel;posrec(2,:,:)=pos;spdrec(2,:,:)=speed;
    colliWall(:,3)=colliWall(:,3)-deltat;colliWall(:,3)=colliWall(:,3)-deltat;

%% Collision with Ball Velocity Analysis 
elseif mint_ball<mint_wall  %hit ball first
    ball1=colliBall(minBallInd,1); ball2=colliBall(minBallInd,2);
    deltat=mint_ball; t=t+deltat;dts(1)=deltat;
    time(2)=t;
    % Simple procedure found in Halie
    r12unit=(pos(ball2,:)-pos(ball1,:))/(norm(pos(ball2,:)-pos(ball1,:)));
    cmnterm=dot((vel(ball1,:)-vel(ball2,:)),r12unit)*r12unit;
    vb1=vel(ball1,:)-cmnterm;vel(ball1,:)=vb1;
    vb2=vel(ball2,:)+cmnterm;vel(ball2,:)=vb2;
    speed(ball1)=norm(vel(ball1,:));speed(ball2)=norm(vel(ball2,:));
    pos=pos+vel*deltat;
    afb(1)=ball1;afb(2)=ball2;Cnum=2;
    for i=1:comb
        xx=colliBall(i,1);yx=colliBall(i,2);
        dist=norm(pos(xx,:)-pos(yx,:));
        if dist<moleculed
            Cnum=Cnum+2;
            afb(Cnum-1)=xx;afb(Cnum)=yx;
            r12unit=(pos(yx,:)-pos(xx,:))/(norm(pos(yx,:)-pos(xx,:)));
            cmnterm=dot((vel(xx,:)-vel(yx,:)),r12unit)*r12unit;
            vb1=vel(xx,:)-cmnterm;vel(xx,:)=vb1;
            vb2=vel(yx,:)+cmnterm;vel(yx,:)=vb2;
            speed(xx)=norm(vel(xx,:));speed(yx)=norm(vel(yx,:));
            posxx=pos(xx,:)+vb1*(moleculed-dist)/speed(xx);pos(xx,:)=posxx;
            posyx=pos(yx,:)+vb2*(moleculed-dist)/speed(yx);pos(yx,:)=posyx;
            ancol=ancol+1;
        end
    end
    volrec(2,:,:)=vel;posrec(2,:,:)=pos;spdrec(2,:,:)=speed;
    colliWall(:,3)=colliWall(:,3)-deltat;colliWall(:,3)=colliWall(:,3)-deltat;

%% Collision with Wall and Ball Simultaneous Velocity Analysis 
else%hit wall and ball simultaneously
    wallid=colliWall(minwallInd,2);
    deltat=colliWall(minwallInd,3);
    t=t+deltat;time(2)=t;dts(1)=deltat;
    if wallid <=2
        vel(minwallInd,1)=-vel(minwallInd,1);
    elseif wallid>2 && wallid <=4
        vel(minwallInd,2)=-vel(minwallInd,2);
    elseif wallid==6%in the corner
        vel(minwallInd,1)=-vel(minwallInd,1);
        vel(minwallInd,2)=-vel(minwallInd,2);
    end

    ball1=colliBall(minBallInd,1); ball2=colliBall(minBallInd,2);
    % Simple procedure found in Halie
    r12unit=(pos(ball2,:)-pos(ball1,:))/(norm(pos(ball2,:)-pos(ball1,:)));
    cmnterm=dot((vel(ball1,:)-vel(ball2,:)),r12unit)*r12unit;
    vb1=vel(ball1,:)-cmnterm;vel(ball1,:)=vb1;
    vb2=vel(ball2,:)+cmnterm;vel(ball2,:)=vb2;
    speed(ball1)=norm(vel(ball1,:));speed(ball2)=norm(vel(ball2,:));
    pos=pos+vel*deltat;
    afb(1)=ball1;afb(2)=ball2;afb(3)=minwallInd;Cnum=3;
    for i=1:comb
        xx=colliBall(i,1);yx=colliBall(i,2);
        dist=norm(pos(xx,:)-pos(yx,:));
        if dist<moleculed
            Cnum=Cnum+2;
            afb(Cnum-1)=xx;afb(Cnum)=yx;
            r12unit=(pos(yx,:)-pos(xx,:))/(norm(pos(yx,:)-pos(xx,:)));
            cmnterm=dot((vel(xx,:)-vel(yx,:)),r12unit)*r12unit;
            vb1=vel(xx,:)-cmnterm;vel(xx,:)=vb1;
            vb2=vel(yx,:)+cmnterm;vel(yx,:)=vb2;
            speed(xx)=norm(vel(xx,:));speed(yx)=norm(vel(yx,:));
            posxx=pos(xx,:)+vb1*(moleculed-dist)/speed(xx);pos(xx,:)=posxx;
            posyx=pos(yx,:)+vb2*(moleculed-dist)/speed(yx);pos(yx,:)=posyx;
            ancol=ancol+1;
        end
    end
    volrec(2,:,:)=vel;posrec(2,:,:)=pos;spdrec(2,:,:)=speed;
    colliWall(:,3)=colliWall(:,3)-deltat;colliWall(:,3)=colliWall(:,3)-deltat;

end
ancols(2)=ancol;
for Simnum=3:(ncol+1)
    for i=1:comb
        x=colliBall(i,1);y=colliBall(i,2);
        c1= concern(x,afb,Cnum);c2=concern(y,afb,Cnum);
        if (c1 || c2)
            rel_pos=pos(x,:)-pos(y,:);
            rel_vel=vel(x,:)-vel(y,:);
            mark1=dot(rel_pos, rel_vel);%Dot Product of rel position and rel velocity
            velsq=dot(rel_vel, rel_vel);
            possq=dot(rel_pos, rel_pos);
            mark2=mark1^2-velsq*(possq-moleculed^2);
            if mark1<0 && mark2>=0
                th1=(-mark1+sqrt(mark2))/velsq; th2=(-mark1-sqrt(mark2))/velsq;
                colliBall(i,3)=compare(th1,th2);
            else
                colliBall(i,3)=inf;
            end
        end
    end
    Cnum=0;
    afb=(num+1)*ones(num,1);
    for j=1:num
            posx=pos(j,1);  posy=pos(j,2);  velx=vel(j,1);  vely=vel(j,2);
            tx=0;ty=0;wx=0;wy=0;
            if (velx>0)
                tx=(boxsize-posx-moleculed/2)/velx;wx=2;
            elseif(velx<0)
                tx=(0-posx+moleculed/2)/velx;wx=1;
            elseif(velx==0)
                tx=inf;wx=5;
            end
            if (vely>0)
                ty=(boxsize-posy-moleculed/2)/vely;wy=4;
            elseif(vely<0)
                ty=(0-posy+moleculed/2)/vely;wy=3;
            elseif(vely==0)
                ty=inf;wy=5;
            end
            if(tx<0)
                tx=0;
            end
            if(ty<0)
                ty=0;
            end
            if (tx<ty)
                colliWall(j,2)=wx;colliWall(j,3)=tx;
            elseif(tx>ty)
                colliWall(j,2)=wy; colliWall(j,3)=ty;
            elseif(tx~=inf && tx==ty)
                colliWall(j,2)=6;colliWall(j,3)=ty;%code 6 collision with both wall together
            else
                colliWall(j,2)=5;colliWall(j,3)=inf;
            end
   end
   temporarystorage1=colliBall(:,3);temporarystorage2=colliWall(:,3);
   [mint_ball, minBallInd] = min(temporarystorage1);
   [mint_wall, minwallInd] = min(temporarystorage2);
   ancol=ancol+1;
   if mint_wall<mint_ball%hit wall first
       wallid=colliWall(minwallInd,2);deltat=mint_wall;t=t+deltat;
       time(Simnum)=t;dts(Simnum-1)=deltat;
       pos=pos+vel*deltat;
   if wallid <=2%x direction
       vel(minwallInd,1)=-vel(minwallInd,1);
   elseif wallid>2 && wallid <=4%y direction
       vel(minwallInd,2)=-vel(minwallInd,2);
   elseif wallid==6%in the corner
       vel(minwallInd,1)=-vel(minwallInd,1);
       vel(minwallInd,2)=-vel(minwallInd,2);
   end
   afb(1)=minwallInd;Cnum=1;
   for i=1:comb
       xx=colliBall(i,1);yx=colliBall(i,2);
       dist=norm(pos(xx,:)-pos(yx,:));
       if dist<moleculed
           Cnum=Cnum+2;
           afb(Cnum-1)=xx;afb(Cnum)=yx;
           r12unit=(pos(yx,:)-pos(xx,:))/(norm(pos(yx,:)-pos(xx,:)));
           cmnterm=dot((vel(xx,:)-vel(yx,:)),r12unit)*r12unit;
           vb1=vel(xx,:)-cmnterm;vel(xx,:)=vb1;
           vb2=vel(yx,:)+cmnterm;vel(yx,:)=vb2;
           speed(xx)=norm(vel(xx,:));speed(yx)=norm(vel(yx,:));
           posxx=pos(xx,:)+vb1*(moleculed-dist)/speed(xx);pos(xx,:)=posxx;
           posyx=pos(yx,:)+vb2*(moleculed-dist)/speed(yx);pos(yx,:)=posyx;
           ancol=ancol+1;
       end
   end
   volrec(Simnum,:,:)=vel;posrec(Simnum,:,:)=pos;spdrec(Simnum,:,:)=speed;
   colliWall(:,3)=colliWall(:,3)-deltat;colliWall(:,3)=colliWall(:,3)-deltat;
%% Collision with Ball Velocity Analysis 
   elseif mint_ball<mint_wall  %hit ball first
       ball1=colliBall(minBallInd,1); ball2=colliBall(minBallInd,2);
       deltat=mint_ball; t=t+deltat;dts(Simnum-1)=deltat;
       time(Simnum)=t;
       % Simple procedure found in Halie
       r12unit=(pos(ball2,:)-pos(ball1,:))/(norm(pos(ball2,:)-pos(ball1,:)));
       cmnterm=dot((vel(ball1,:)-vel(ball2,:)),r12unit)*r12unit;
       vb1=vel(ball1,:)-cmnterm;vel(ball1,:)=vb1;
       vb2=vel(ball2,:)+cmnterm;vel(ball2,:)=vb2;
       speed(ball1)=norm(vel(ball1,:));speed(ball2)=norm(vel(ball2,:));
       pos=pos+vel*deltat;
        afb(1)=ball1;afb(2)=ball2;Cnum=2;
        for i=1:comb
            xx=colliBall(i,1);yx=colliBall(i,2);
            dist=norm(pos(xx,:)-pos(yx,:));
            if dist<moleculed
                Cnum=Cnum+2;
                afb(Cnum-1)=xx;afb(Cnum)=yx;
                r12unit=(pos(yx,:)-pos(xx,:))/(norm(pos(yx,:)-pos(xx,:)));
                cmnterm=dot((vel(xx,:)-vel(yx,:)),r12unit)*r12unit;
                vb1=vel(xx,:)-cmnterm;vel(xx,:)=vb1;
                vb2=vel(yx,:)+cmnterm;vel(yx,:)=vb2;
                speed(xx)=norm(vel(xx,:));speed(yx)=norm(vel(yx,:));
                posxx=pos(xx,:)+vb1*(moleculed-dist)/speed(xx);pos(xx,:)=posxx;
                posyx=pos(yx,:)+vb2*(moleculed-dist)/speed(yx);pos(yx,:)=posyx;
                ancol=ancol+1;
            end
        end
       volrec(Simnum,:,:)=vel;posrec(Simnum,:,:)=pos;spdrec(Simnum,:,:)=speed;
       colliWall(:,3)=colliWall(:,3)-deltat;colliWall(:,3)=colliWall(:,3)-deltat;

%% Collision with Wall and Ball Simultaneous Velocity Analysis 
    else%hit wall and ball simultaneously
        wallid=colliWall(minwallInd,2);
        deltat=colliWall(minwallInd,3);
        t=t+deltat;time(Simnum)=t;dts(Simnum-1)=deltat;
        if wallid <=2
            vel(minwallInd,1)=-vel(minwallInd,1);
        elseif wallid>2 && wallid <=4
            vel(minwallInd,2)=-vel(minwallInd,2);
        elseif wallid==6%in the corner
            vel(minwallInd,1)=-vel(minwallInd,1);
            vel(minwallInd,2)=-vel(minwallInd,2);
        end
        ball1=colliBall(minBallInd,1); ball2=colliBall(minBallInd,2);
        % Simple procedure found in Halie
        r12unit=(pos(ball2,:)-pos(ball1,:))/(norm(pos(ball2,:)-pos(ball1,:)));
        cmnterm=dot((vel(ball1,:)-vel(ball2,:)),r12unit)*r12unit;
        vb1=vel(ball1,:)-cmnterm;vel(ball1,:)=vb1;
        vb2=vel(ball2,:)+cmnterm;vel(ball2,:)=vb2;
        speed(ball1)=norm(vel(ball1,:));speed(ball2)=norm(vel(ball2,:));
        pos=pos+vel*deltat;
        afb(1)=ball1;afb(2)=ball2;afb(3)=minwallInd;Cnum=3;
        for i=1:comb
            xx=colliBall(i,1);yx=colliBall(i,2);
            dist=norm(pos(xx,:)-pos(yx,:));
            if dist<moleculed
                Cnum=Cnum+2;
                afb(Cnum-1)=xx;afb(Cnum)=yx;
                r12unit=(pos(yx,:)-pos(xx,:))/(norm(pos(yx,:)-pos(xx,:)));
                cmnterm=dot((vel(xx,:)-vel(yx,:)),r12unit)*r12unit;
                vb1=vel(xx,:)-cmnterm;vel(xx,:)=vb1;
                vb2=vel(yx,:)+cmnterm;vel(yx,:)=vb2;
                speed(xx)=norm(vel(xx,:));speed(yx)=norm(vel(yx,:));
                posxx=pos(xx,:)+vb1*(moleculed-dist)/speed(xx);pos(xx,:)=posxx;
                posyx=pos(yx,:)+vb2*(moleculed-dist)/speed(yx);pos(yx,:)=posyx;
                ancol=ancol+1;
            end
        end
        volrec(Simnum,:,:)=vel;posrec(Simnum,:,:)=pos;spdrec(Simnum,:,:)=speed;
        colliWall(:,3)=colliWall(:,3)-deltat;colliWall(:,3)=colliWall(:,3)-deltat;
   end
   ancols(Simnum)=ancol;
end
        
%% Plot Figures and Save as images
min(dts)
temp=1/2;
figure('Color', 'white','units','normalized','outerposition',[0 0 1 1])
p1=subplot(2,2,1);
h1 = plot(posrec(ncol+1,:,1),posrec(ncol+1,:,2),...
	'o','Color','red','MarkerFaceColor','red','MarkerSize',2);
xlim([0 21]);ylim([0 21]);interval=0.05;
edges=(-5:interval:5);edgeb=(0:interval:5);
coe1=(sqrt(molecuelm/(2*pi*temp*bconst)));%probability histogram vx vy
coe3=-(molecuelm/(2*temp*bconst));%exp stuff
coe4=(molecuelm/(temp*bconst));%probability histogram speed
curve1=zeros(length(edges),1);
for i=1:length(edges)
    curve1(i)=coe1*exp(coe3*edges(i)^2);
end
curve3=zeros(length(edgeb),1);
for i=1:length(edgeb)
    curve3(i)=coe4*edgeb(i)*exp(coe3*edgeb(i)^2);
end
simtitle=strcat('Simulation at time   ::   ', string(time(ncol+1)),'    Number of collisions   ::', string(ancol));
p1.XLabel.String = simtitle;
p2=subplot(2,2,2);
g =histogram(volrec(ncol+1,:,1),edges,'Normalization','pdf'); hold on;
cp1=plot(edges, curve1); hold on;
p2.XLabel.String = 'Vx (u/t)';p2.YLabel.String = 'Distribution';
p3=subplot(2,2,3);
k1 =histogram(volrec(ncol+1,:,2),edges,'Normalization','pdf');hold on;
cp2=plot(edges, curve1); hold on;
p3.XLabel.String = 'Vy (u/t)';p3.YLabel.String = 'Distribution';
p4=subplot(2,2,4);
l1=histogram(spdrec(ncol+1,:,1),edgeb,'Normalization','pdf');hold on;
cp3=plot(edgeb, curve3); hold on;
p4.XLabel.String = 'Speed (u/t)';p4.YLabel.String = 'Distribution';
%%
inframe=0;
for i=1:num
    if(pos(i,1)>=0 && pos(i,1)<=21 && pos(i,2)>=0 && pos(i,2)<=21)
        inframe=inframe+1;
    end
end
inframe
energy=(1/2)*sum(sum(speed.^2,2))
figure('Color', 'white','units','normalized','outerposition',[0 0 1 1])
p1=subplot(1,1,1);
h = plot(posrec(1,:,1),posrec(1,:,2),...
	'o','Color','red','MarkerFaceColor','red','MarkerSize',3);
xlim([0 21]);ylim([0 21]);interval=0.05;
simtitle=strcat('Simulation at time   ::   ', string(time(1)),'    Number of collisions   ::', string(0));
p1.XLabel.String = simtitle;
for i=2:ncol+1
    pause(0.1);
    h.XData=posrec(i,:,1);
    h.YData=posrec(i,:,2);
    simtitle=strcat('Simulation at time  ::   ', string(time(i)),'    Number of collisions   ::', string(i));
    p1.XLabel.String = simtitle;
end
%%
figure('Color', 'white','units','normalized','outerposition',[0 0 1 1])
p1=subplot(2,2,1);
h = plot(posrec(1,:,1),posrec(1,:,2),...
	'o','Color','red','MarkerFaceColor','red','MarkerSize',2);
xlim([0 21]);ylim([0 21]);interval=0.05;
simtitle=strcat('Simulation at time   ::   ', string(time(1)),'    Number of collisions   ::', string(ancols(1)));
p1.XLabel.String = simtitle;
p2=subplot(2,2,2);
g =histogram(volrec(1,:,1),edges,'Normalization','pdf');hold on;
cp1=plot(edges, curve1); hold on;
p2.XLabel.String = 'Vx (u/t)';p2.YLabel.String = 'Distribution';
p3=subplot(2,2,3);
k =histogram(volrec(1,:,2),edges,'Normalization','pdf');hold on;
cp2=plot(edges, curve1); hold on;
p3.XLabel.String = 'Vy (u/t)';p3.YLabel.String = 'Distribution';
p4=subplot(2,2,4);
l=histogram(spdrec(1,:,1),edgeb,'Normalization','pdf');hold on;
cp3=plot(edgeb, curve3); hold on;
p4.XLabel.String = 'Speed (u/t)';p4.YLabel.String = 'Distribution';
figtitle='Initial Condition';
saveas(gcf,figtitle,'png')
for i=2:ncol+1
    pause(0.002);
    h.XData=posrec(i,:,1);
    h.YData=posrec(i,:,2);
    simtitle=strcat('Simulation at time  ::   ', string(time(i)),'    Number of collisions   ::', string(ancols(i)));
    p1.XLabel.String = simtitle;
    subplot(2,2,2)
    hold off;
    g =histogram(volrec(i,:,1),edges,'Normalization','pdf');hold on;
    cp1=plot(edges, curve1); hold on;
    p2.XLabel.String = 'Vx (u/t)';
    p2.YLabel.String = 'Distribution';
    subplot(2,2,3)
    hold off;
    k =histogram(volrec(i,:,2),edges,'Normalization','pdf');hold on;
    cp2=plot(edges, curve1); hold on;
    p3.XLabel.String = 'Vy (u/t)';
    p3.YLabel.String = 'Distribution';
    subplot(2,2,4)
    hold off;
    l=histogram(spdrec(i,:,1),edgeb,'Normalization','pdf');hold on;
    cp3=plot(edgeb, curve3); hold on;
    p4.XLabel.String = 'Speed (u/t)';
    p4.YLabel.String = 'Distribution';
    if mod((i-1),200)==0
        figtitle=strcat('collision number',string(i-1));
        saveas(gcf,figtitle,'png')
    end
end
figtitle=strcat('collision number',string(ancols(i)));
saveas(gcf,figtitle,'png')

function toupdate = concern(i,afb,Cnum)
    toupdate=false;
    for j=1:Cnum
        if(i==afb(j))
            toupdate=true;
            break
        end
    end
end
function approtime=compare(time1,time2)
    if (time1>=0 && time2>=0)
        approtime=min(time1,time2);
    elseif(time1<0 && time2>=0)
        approtime=time2;
    elseif(time2<0 && time1>=0)
        approtime=time1;
    else
        approtime=inf;
    end
end
