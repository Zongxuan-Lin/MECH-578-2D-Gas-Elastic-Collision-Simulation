clc;clear all;close all;
t=0; %time start at 0
ancol=0;
boxsize=21; %box size
num=(boxsize-1)^2;%A total of 400 molecules in the simulation
bconst=1.0; %fictional boltzman constant in the problem
temp=1/2;%initial temperature
ncol=10;%number of collisions
moleculed=0.1; %diameter and mass of molecule
molecuelm=1;
m2=molecuelm*2;
rangepair=[-5 -3 0;5 7 10];
bolfac=4;%size of sample scaling factor for boltz

intervals3=triangle([-2 0 2], num*bolfac);%for boltzmann distribution function
intervals4=triangle([-2 0 2], num*bolfac);%for boltzmann distribution function
%velbol=impulse(2,num*bolfac);%boltzman velocity 
%velbol=settemp(velbol,temp,num*bolfac);%boltzman velocity 
velbol=radinrange(intervals3,intervals4,num*bolfac);
%% Boltzmann Equation Distribution
vstd=sqrt(bconst*temp/molecuelm);
ntcit=20; %number of montecalo sim for each collision
maxiniV=max(abs(velbol),[],'all');
numsec=30;%number of steps for vel discretization in each direction
dom=3*vstd;%initialize dimension of velocity space to 3*std which is 99% confidence zone
if maxiniV>dom
    dom=maxiniv;
end
stepsize=dom*2/numsec;
velrec=linspace(-dom,dom,numsec+1);
[X,Y] = meshgrid(velrec);
xpop=zeros(1,numsec+1);ypop=zeros(1,numsec+1);
N=num*bolfac;
Nuv=zeros(numsec+1);
for i =1:N
    coor1=round((velbol(i,:)+dom)/stepsize);
    Nuv(coor1(1),coor1(2))=Nuv(coor1(1),coor1(2))+1;
    for k = 1:numsec+1
         if (velbol(i,1)>=(velrec(k)-stepsize/2)) && (velbol(i,1)< (velrec(k)+stepsize/2))
            xpop(k) = xpop(k)+1;
         end
         if (velbol(i,2)>=(velrec(k)-stepsize/2)) && (velbol(i,2)< (velrec(k)+stepsize/2))
            ypop(k) = ypop(k)+1;%distribution statistics
        end
    end
end
f_u = xpop./(N) ;%distribution probability
f_v = ypop./(N) ;
fuv=Nuv/N;
surf(X,Y,fuv);
cinit=zeros(numsec+1);
size1=numsec+1;
dt=1;
for times=1:50
    cinit=zeros(numsec+1);
    for i=1:numsec+1
        for j=1:numsec+1
            term1=colli(i,j,fuv, size1,dom,ntcit,stepsize,velrec,moleculed);
            term1=term1*dom^2*moleculed;
            cinit(i,j)=term1;
        end
    end
    fuv=fuv+cinit*dt;
    fuv=fuv*1/sum(fuv,'all');
contourf(X,Y,fuv);

end
for times=1:50
    cinit=zeros(numsec+1);
    for i=1:numsec+1
        for j=1:numsec+1
            term1=colli(i,j,fuv, size1,dom,ntcit,stepsize,velrec,moleculed);
            term1=term1*dom^2*moleculed;
            cinit(i,j)=term1;
        end
    end
    fuv=fuv+cinit*dt;
    fuv=fuv*1/sum(fuv,'all');
contourf(X,Y,fuv);

end
for times=1:50
    cinit=zeros(numsec+1);
    for i=1:numsec+1
        for j=1:numsec+1
            term1=colli(i,j,fuv, size1,dom,ntcit,stepsize,velrec,moleculed);
            term1=term1*dom^2*moleculed;
            cinit(i,j)=term1;
        end
    end
    fuv=fuv+cinit*dt;
    fuv=fuv*1/sum(fuv,'all');
contourf(X,Y,fuv);

end
1
function term1=colli(i,j,fuv, size1,dom,ntcit,stepsize,velrec,moleculed)
fav=fuv(i,j);
va=[velrec(i); velrec(j)];
term1=0;
db=linspace(-moleculed,moleculed,ntcit+1);
for m=1:size1
    for n=1:size1
        if m~=i && n~=j
            fbv=fuv(m,n);
            vb=[velrec(m); velrec(n)];
            ctrv=(va+vb)/2;
            for k=1:ntcit
                b=db(k);
                direc=asin(b/moleculed);
                kgn=[cos(direc);sin(direc)];
                t1=sqrt(1-(b/moleculed)^2);t2=b/moleculed;
                g=vb-va;
                %vasl=va+kgn.*(kgn'*g);
                vbsl=t2*[t2 -t1;t1 t2]*g+vb;%v1 aki
                vasl=t1*[t1 t2;-t2 t1]*g+vb;%v2 aki
                cora=(vasl+dom)/stepsize;
                ca=round(cora);
                if ca(1)>length(velrec)
                    ca(1)=length(velrec);
                elseif ca(1)<=1
                    ca(1)=1;
                end
                if ca(2)>length(velrec)
                    ca(2)=length(velrec);
                elseif ca(2)<=1
                    ca(2)=1;
                end
                corc=(ctrv+dom)/stepsize;

                %vbsl=vb-kgn.*(kgn'*g);
                corb=(vbsl+dom)/stepsize;
                cb=round(corb);
                if cb(1)>length(velrec)
                    cb(1)=length(velrec);
                elseif cb(1)<=1
                    cb(1)=1;
                end
                if cb(2)>length(velrec)
                    cb(2)=length(velrec);
                elseif cb(2)<=1
                    cb(2)=1;
                end
                curres=(fav*fbv-fuv(ca(1),ca(2))*fuv(ca(1),cb(2)))*stepsize^2;
                 [ina, outa]=relpos(corc,cora,velrec);
                 [inb, outb]=relpos(corc,corb,velrec);
%                 e1=(velrec(ina(1))^2+velrec(ina(2))^2+velrec(inb(1))^2+velrec(inb(2))^2);
%                 e2=(velrec(outa(1))^2+velrec(outa(2))^2+velrec(outb(1))^2+velrec(outb(2))^2);
%                 e0=sum(vasl.^2)+sum(vbsl.^2);
%                 rv=(e0-e1)/(e2-e1);
                 dv1=fuv(ina(1),ina(2))*fuv(inb(1),inb(2));
                 dv2=fuv(outa(1),outa(2))*fuv(outb(1),outb(2));
                dv0=fbv*fav;
                    rv=0.5;
                 fafabsl=0.5*dv1+0.5*dv2;
%                 rvstar=rv*dv1/(rv*dv1+(1-rv)*dv2);
                 curres=(-1+(1-rv)*((ina(1)==i && ina(2)==j)+(inb(1)==i && inb(2)==j))+rv*((outa(1)==i && outa(2)==j)+(outb(1)==i && outb(2)==j)))*(dv0-fafabsl);
%                 curres=((1-rvstar)*dv1+rvstar*dv2-dv0)+...
%                 ((ina(1)==i && ina(2)==j)+(inb(1)==i && inb(2)==j))*((1-rv)*dv0-(1-rvstar)*dv1)+...
%                 ((outa(1)==i && outa(2)==j)+(outb(1)==i && outb(2)==j))*(rv*dv0-rvstar*dv2);
                curres=curres*2*moleculed*norm(g)/ntcit;%*stepsize*stepsize;
                if ~isnan(curres)
                    term1=term1+curres;
                end
            end
        end
    end
end
end



function [in, out]=relpos(ctrv,vsl,velrec)
ro=vsl-ctrv;
a=floor(vsl);ra=a-ctrv;dra=sqrt(sum(ra.^2))-sqrt(sum(ro.^2));
b=ceil(vsl);rb=b-ctrv;drb=sqrt(sum(rb.^2))-sqrt(sum(ro.^2));
c=[a(1);b(2)];rc=c-ctrv;drc=sqrt(sum(rc.^2))-sqrt(sum(ro.^2));
d=[a(2);b(1)];rd=d-ctrv;drd=sqrt(sum(rd.^2))-sqrt(sum(ro.^2));
alco=[a b c d];
ad=[sqrt(sum(ra.^2)) sqrt(sum(rb.^2)) sqrt(sum(rc.^2)) sqrt(sum(rd.^2))];
cor=[1 2;1 3;1 4;2 3;2 4;3 4];
judge=[dra*drb dra*drc dra*drd drb*drc drb*drd drc*drd];
[M,I] = min(judge); 
ind=cor(I,:);
d1=ad(ind(1));d2=ad(ind(2));
if d1<=d2
    in=alco(:,ind(1));
    if in(1)>length(velrec)
        in(1)=length(velrec);
    elseif in(1)<=1
        in(1)=1;
    end
    if in(2)>length(velrec)
        in(2)=length(velrec);
    elseif in(2)<=1
        in(2)=1;
    end
    out=alco(:,ind(2));
    if out(1)>length(velrec)
        out(1)=length(velrec);
    elseif out(1)<=1
        out(1)=1;
    end
    if out(2)>length(velrec)
        out(2)=length(velrec);
    elseif out(2)<=1
        out(2)=1;
    end
else
    in=alco(:,ind(2));
    out=alco(:,ind(1));
    if in(1)>length(velrec)
        in(1)=length(velrec);
    elseif in(1)<=1
        in(1)=1;
    end
    if in(2)>length(velrec)
        in(2)=length(velrec);
    elseif in(2)<=1
        in(2)=1;
    end
    if out(1)>length(velrec)
        out(1)=length(velrec);
    elseif out(1)<=1
        out(1)=1;
    end
    if out(2)>length(velrec)
        out(2)=length(velrec);
    elseif out(2)<=1
        out(2)=1;
    end
end
end




%% Different Distribution Initialization
%determine range of velocity randomization for 
function intervals=tophat(rangepair, num)
%range pair n*m m=2, then it is tophats, each interval is ranked from
%if m=3, the third parameter should be population in the subrange
%smallest to largest
intervals=zeros(num,4);
if size(rangepair,2)==2 && size(rangepair,1)>1
    tempr=linspace(0,sum(rangepair(:,2)-rangepair(:,1)),num+1);
    intervals(:,1)=tempr(1:num)+rangepair(1,1);intervals(:,2)=tempr(2:num+1)+rangepair(1,1);
    for i=2:size(rangepair,1)
        decision=(intervals(:,2)>rangepair(i-1,2))+(intervals(:,1)>rangepair(i-1,2));
        [~, idx] = max((decision==1));
        temp=intervals(idx,2);intervals(idx,2)=rangepair(i-1,2);
        intervals(idx,3)=rangepair(i,1);
        intervals(idx,4)=rangepair(i,1)+temp-rangepair(i-1,2);
        intervals(:,1)=intervals(:,1)+(decision==2)*(rangepair(i,1)-rangepair(i-1,2));
        intervals(:,2)=intervals(:,2)+(decision==2)*(rangepair(i,1)-rangepair(i-1,2));
    end
elseif size(rangepair,2)==2 && size(rangepair,1)==1
    tempr=linspace(0,(rangepair(2)-rangepair(1)),num+1);
    intervals(:,1)=tempr(1:num)+rangepair(1);intervals(:,2)=tempr(2:num+1)+rangepair(1);
elseif size(rangepair,2)==3
    step=1;
    for i=1:size(rangepair,1)
        tempr=linspace(rangepair(i,1),rangepair(i,2),rangepair(i,3)+1);
        intervals(step:rangepair(i,3)+step-1,1)=tempr(1:rangepair(i,3));
        intervals(step:rangepair(i,3)+step-1,2)=tempr(2:rangepair(i,3)+1);
        step=rangepair(i,3)+step;
    end
end
end
function intervals=triaglint(rangepair,tempr1,height)
%range pair[a,b,c] not for right angle triangles
    base1=(rangepair(2)-rangepair(1));base2=(rangepair(3)-rangepair(2));
    tempr=sqrt(tempr1*2*base1/height)+rangepair(1);
    num=length(tempr1)-1;
    intervals=zeros(num,2);
    intervals(:,1)=tempr(1:num);intervals(:,2)=tempr(2:num+1);
    decision=(intervals(:,2)>rangepair(2))+(intervals(:,1)>rangepair(2));
    [~, idx] = max((decision==1));
    intervals(idx,2)=rangepair(2);
    leftarea=base1*height/2;
    decision2=(tempr1-leftarea)>0;
    tempr2=decision2.*base2*height/2-(tempr1-leftarea).*decision2;
    tempr2=-sqrt(tempr2*2*base2/height)+rangepair(3)*decision2;
    intervals2=zeros(num,2);
    intervals2(:,1)=tempr2(1:num);intervals2(:,2)=tempr2(2:num+1);
    intervals(idx,2)=intervals2(idx,2);intervals2(idx,2)=0;intervals2(idx,1)=0;
    intervals(idx+1:num,:)=0;
    intervals=intervals+intervals2;
    intervals(num,2)=rangepair(3);
end
function intervals=thint(rangepair,tempr,height)
    num=length(tempr)-1;
    tempr=tempr/height+rangepair(1);
    intervals=zeros(num,2);
    intervals(:,1)=tempr(1:num);intervals(:,2)=tempr(2:num+1);
end
function intervals=triangle(rangepair, num) 
%range pair n*m m=3, then it is tophats, each interval is ranked from
%smallest to largest like y coordinate of the 3 points from smallest to
%largest
%if m=4 the forth coordinatewould be the population in that range
intervals=zeros(num,4);
if size(rangepair,2)==3 && size(rangepair,1)==1
    height=num*2/(rangepair(3)-rangepair(1));
    base1=(rangepair(2)-rangepair(1));base2=(rangepair(3)-rangepair(2));
    tempr=linspace(0,400,num+1);
    tempr=sqrt(tempr*2*base1/height)+rangepair(1);
    intervals(:,1)=tempr(1:num);intervals(:,2)=tempr(2:num+1);
    decision=(intervals(:,2)>rangepair(2))+(intervals(:,1)>rangepair(2));
    [~, idx] = max((decision==1));
    intervals(idx,2)=rangepair(2);
    stepsleft=num-idx+1;
    tempr2=linspace(stepsleft,0,stepsleft+1);
    tempr2=-sqrt(tempr2*2*base2/height)+rangepair(3);
    intervals2=zeros(num,4);
    intervals2(idx:num,1)=tempr2(1:length(tempr2)-1);intervals2(idx:num,2)=tempr2(1:length(tempr2)-1);
    intervals(idx,2)=intervals2(idx,2);intervals2(idx,2)=0;intervals2(idx,1)=0;
    intervals(idx+1:num,:)=0;intervals=intervals+intervals2;
elseif size(rangepair,2)==3 && size(rangepair,1)>1
    height=num*2/sum(rangepair(:,3)-rangepair(:,1));
    areas=height*(rangepair(:,3)-rangepair(:,1))/2;
    loc=1;
    st=1;
    for i=1:length(areas)
        atc=linspace(0,floor(areas(i)),floor(areas(i))+1);
        if st~=1
            atc=atc+st;
            atc=[0 atc];
        end
        if i~=length(areas)  && (floor(areas(i)))~=areas(i)
            atc=[atc areas(i)];
        end
        st=1-(areas(i)-floor(areas(i))-st);
        if st==2
            st=1;
        end
        subset=triaglint(rangepair(i,:),atc,height);
        if loc==1
            intervals(loc:length(subset),1)=subset(:,1);
            intervals(loc:length(subset),2)=subset(:,2);
            loc=length(subset);
        else
            if st~=1
                intervals(loc,3)=subset(1,1);
                intervals(loc,4)=subset(1,2);
                intervals(loc+1:length(subset)+loc-1,1)=subset(2:length(subset),1);
                intervals(loc+1:length(subset)+loc-1,2)=subset(2:length(subset),2);
            else
                intervals(loc+1:length(subset)+loc,1)=subset(1:length(subset),1);
                intervals(loc+1:length(subset)+loc,2)=subset(1:length(subset),2);
            end
        end
    end
    st=0;
elseif size(rangepair,2)==4 % the forth parameter is the number of particle in the range
        step=1;
    for i=1:size(rangepair,1)
        tempr=0:1:rangepair(i,4);
        height=rangepair(i,4)*2/(rangepair(i,3)-rangepair(i,1));
        subset=triaglint(rangepair(i,1:3),tempr,height);
        intervals(step:rangepair(i,4)+step-1,1)=subset(:,1);
        intervals(step:rangepair(i,4)+step-1,2)=subset(:,2);
        step=rangepair(i,4)+step;
    end
end
end
function intervals=trapezoid(rangepair, num)
%range pair n*m m=4, then it is tophats, each interval is ranked from
%smallest to largest like y coordinate of the 4 points from smallest to
%largest
%if m=5 the fifth coordinatewould be the population in that range
intervals=zeros(num,4);    
if size(rangepair,2)==4 && size(rangepair,1)==1
    height=num/((rangepair(4)-rangepair(1)+rangepair(3)-rangepair(2))/2);
    a1=height*(rangepair(2)-rangepair(1))/2;
    a2=height*(rangepair(3)-rangepair(2));
    a3=height*(rangepair(4)-rangepair(3))/2;
    step=1;
    st=0;
    if a1~=0
        atc=linspace(0,floor(a1),floor(a1)+1);
        if (floor(a1))~=a1
            atc=[atc a1];
        end
        st=1-(a1-floor(a1)-st);
        subset1=rtriaglint(rangepair(1:2),atc,height);
        intervals(step:length(subset1)+step-1,1)=subset1(:,1);
        intervals(step:length(subset1)+step-1,2)=subset1(:,2);
        step=length(subset1);
    end
    atc=linspace(0,floor(a2),floor(a2)+1);
    if st==2
        st=1;
    end
    if st~=1
        atc=atc+st;
        atc=[0 atc];
    end
    if (floor(a2))~=a2
        atc=[atc a2];
    end
    subset2=thint(rangepair(2:3),atc,height);
    if st~=1
        intervals(step,2)=subset2(1,2);
        intervals(step+1:length(subset2)+step-1,1)=subset2(2:length(subset2),1);
        intervals(step+1:length(subset2)+step-1,2)=subset2(2:length(subset2),2);
        
    else
        st=1;
        intervals(step+1:length(subset2)+step,1)=subset2(1:length(subset2),1);
        intervals(step+1:length(subset2)+step,2)=subset2(1:length(subset2),2);
    end
    step=length(subset2)+step;
    st=1-(a2-floor(a2)-st);
    if st==2
        st=1;
    end
    if a3~=0
        atc=linspace(0,floor(a3),floor(a3)+1);
        if st~=1
            atc=atc+st;
            atc=[0 atc];
        end
        if (floor(a3))~=a3
            atc=[atc a3];
        end
        subset3=ltrianglint(rangepair(3:4),atc,height);
        if st~=1
            intervals(step,2)=subset3(1,2);
            intervals(step+1:length(subset3)+step-1,1)=subset3(2:length(subset3),1);
            intervals(step+1:length(subset3)+step-1,2)=subset3(2:length(subset3),2);
        else
            st=1;
            intervals(step+1:length(subset3)+step,1)=subset3(1:length(subset3),1);
            intervals(step+1:length(subset3)+step,2)=subset3(1:length(subset3),2);
        end
    end
    step=length(subset3)+step;
    st=1-(a3-floor(a3)-st);
elseif size(rangepair,2)==4 && size(rangepair,1)>1
        height=num/sum((rangepair(:,4)-rangepair(:,1)+rangepair(:,3)-rangepair(:,2))/2);
        st=1;step=1;
        for i=1:size(rangepair,1) 
            a1=height*(rangepair(i,2)-rangepair(i,1))/2;
            a2=height*(rangepair(i,3)-rangepair(i,2));
            a3=height*(rangepair(i,4)-rangepair(i,3))/2;
            if a1~=0
                atc=linspace(0,floor(a1),floor(a1)+1);
                if (floor(a1))~=a1
                    atc=[atc a1];
                end
                subset1=rtriaglint(rangepair(i,1:2),atc,height);
                if i==1
                    intervals(step:length(subset1)+step-1,1)=subset1(:,1);
                    intervals(step:length(subset1)+step-1,2)=subset1(:,2);
                    step=length(subset1);
                else
                    if st~=1
                        intervals(step,3)=subset1(1,1);intervals(step,4)=subset1(1,2);
                        intervals(step+1:length(subset1)+step-1,1)=subset1(2:length(subset1),1);
                        intervals(step+1:length(subset1)+step-1,2)=subset1(2:length(subset1),2);
                    else
                        intervals(step+1:length(subset1)+step,1)=subset1(1:length(subset1),1);
                        intervals(step+1:length(subset1)+step,2)=subset1(1:length(subset1),2);
                    end
                    step=length(subset1)+step;
                end
            end
            st=1-(a1-floor(a1)-st);
            atc=linspace(0,floor(a2),floor(a2)+1);
            if st==2
                st=1;
            end
            if st~=1
                atc=atc+st;
                atc=[0 atc];
            end
            if (floor(a2))~=a2
                atc=[atc a2];
            end
            subset2=thint(rangepair(i,2:3),atc,height);
            if st~=1
                intervals(step,2)=subset2(1,2);
                intervals(step+1:length(subset2)+step-1,1)=subset2(2:length(subset2),1);
                intervals(step+1:length(subset2)+step-1,2)=subset2(2:length(subset2),2);
            else
                intervals(step+1:length(subset2)+step,1)=subset2(1:length(subset2),1);
                intervals(step+1:length(subset2)+step,2)=subset2(1:length(subset2),2);
            end
            step=length(subset2)+step;
            st=1-(a2-floor(a2)-st);
            if st==2
                st=1;
            end
            if a3~=0
                atc=linspace(0,floor(a3),floor(a3)+1);
                if st~=1
                    atc=atc+st;
                    atc=[0 atc];
                end
                if (floor(a3))~=a3
                    atc=[atc a3];
                end
                subset3=ltrianglint(rangepair(i,3:4),atc,height);
                if st~=1
                    intervals(step,2)=subset3(1,2);
                    intervals(step+1:length(subset3)+step-1,1)=subset3(2:length(subset3),1);
                    intervals(step+1:length(subset3)+step-1,2)=subset3(2:length(subset3),2);
                else
                    st=1;
                    intervals(step+1:length(subset3)+step,1)=subset3(1:length(subset3),1);
                    intervals(step+1:length(subset3)+step,2)=subset3(1:length(subset3),2);
                end
            end
            step=length(subset3)+step;
            st=1-(a3-floor(a3)-st);
            if st==2
                st=1;
            end
        end
elseif size(rangepair,2)==5
        height=num/sum((rangepair(:,4)-rangepair(:,1)+rangepair(:,3)-rangepair(:,2))/2);
        step=1;
        for i=1:size(rangepair,1) 
            st=1;
            height=rangepair(i,5)/((rangepair(i,4)-rangepair(i,1)+rangepair(i,3)-rangepair(i,2))/2);
            a1=height*(rangepair(i,2)-rangepair(i,1))/2;
            a2=height*(rangepair(i,3)-rangepair(i,2));
            a3=height*(rangepair(i,4)-rangepair(i,3))/2;
            if a1~=0
                atc=linspace(0,floor(a1),floor(a1)+1);
                if (floor(a1))~=a1
                    atc=[atc a1];
                end
                subset1=rtriaglint(rangepair(i,1:2),atc,height);
                if i==1
                    intervals(step:length(subset1)+step-1,1)=subset1(:,1);
                    intervals(step:length(subset1)+step-1,2)=subset1(:,2);
                    step=length(subset1);
                else
                    intervals(step+1:length(subset1)+step,1)=subset1(1:length(subset1),1);
                    intervals(step+1:length(subset1)+step,2)=subset1(1:length(subset1),2);
                    step=length(subset1)+step;
                end
            end
            st=1-(a1-floor(a1)-st);
            atc=linspace(0,floor(a2),floor(a2)+1);
            if st==2
                st=1;
            end
            if st~=1
                atc=atc+st;
                atc=[0 atc];
            end
            if (floor(a2))~=a2
                atc=[atc a2];
            end
            subset2=thint(rangepair(i,2:3),atc,height);
            if st~=1
                intervals(step,2)=subset2(1,2);
                intervals(step+1:length(subset2)+step-1,1)=subset2(2:length(subset2),1);
                intervals(step+1:length(subset2)+step-1,2)=subset2(2:length(subset2),2);
            else
                intervals(step+1:length(subset2)+step,1)=subset2(1:length(subset2),1);
                intervals(step+1:length(subset2)+step,2)=subset2(1:length(subset2),2);
            end
            step=length(subset2)+step;
            st=1-(a2-floor(a2)-st);
            if st==2
                st=1;
            end
            if a3~=0
                atc=linspace(0,floor(a3),floor(a3)+1);
                if st~=1
                    atc=atc+st;
                    atc=[0 atc];
                end
                if (floor(a3))~=a3
                    atc=[atc a3];
                end
                subset3=ltrianglint(rangepair(i,3:4),atc,height);
                if st~=1
                    intervals(step,2)=subset3(1,2);
                    intervals(step+1:length(subset3)+step-1,1)=subset3(2:length(subset3),1);
                    intervals(step+1:length(subset3)+step-1,2)=subset3(2:length(subset3),2);
                else
                    st=1;
                    intervals(step+1:length(subset3)+step,1)=subset3(1:length(subset3),1);
                    intervals(step+1:length(subset3)+step,2)=subset3(1:length(subset3),2);
                end
            end
            step=length(subset3)+step;
        end        
end
end
function intervals=rtriaglint(rangepair,tempr1,height)
%range pair [a,b]for triangles with right angle on the right
    base1=(rangepair(2)-rangepair(1));
    tempr=sqrt(tempr1*2*base1/height)+rangepair(1);
    num=length(tempr1)-1;
    intervals=zeros(num,2);
    intervals(:,1)=tempr(1:num);intervals(:,2)=tempr(2:num+1);
end
function intervals=ltrianglint(rangepair,tempr1,height)
%range pair [a b] for triangle with right angle on the left
    base2=(rangepair(2)-rangepair(1));
    totarea=base2*height/2;
    tempr2=totarea-tempr1;
    num=length(tempr1)-1;
    tempr2=-sqrt(tempr2*2*base2/height)+rangepair(2);
    intervals=zeros(num,2);
    intervals(:,1)=tempr2(1:num);intervals(:,2)=tempr2(2:num+1);
end
function xvvel=impulse(rangepair, num)
% rangepair [a,b]*n  is the speed b is the population with the speed
% if range pair is just a number,  then all population should have the same
% speed
% with initial speed directly output velocitys
xvvel=zeros(num,2);
spd=zeros(num,1);
inidir=rand(1,num)*2*pi;
if size(rangepair,1)==2
    step=1;
    for i=1:size(rangepair,2)
        spd(step:rangepair(i,2)+step-1)=ones(1,rangepair(i,2))*rangepair(i,1);
        step=rangepair(i,2)+step;
    end
    xvvel(:,1)=spd.*cos(inidir)';xvvel(:,2)=spd.*sin(inidir)';
elseif size(rangepair,1)==1
    xvvel(:,1)=rangepair*cos(inidir);xvvel(:,2)=rangepair*sin(inidir);
end
end
function xvvel=randomdis(range, num)
%range [a b]
inidir=rand(1,num)*2*pi;
spd=rand(num,1)*(range(2)-range(1))+range(1);
xvvel(:,1)=spd.*cos(inidir)';xvvel(:,2)=spd.*sin(inidir)';    
end
function xvvel=settemp(vel,temp,num)
xvvel=vel;
sumvel=sum(vel.^2,"all");
curtemp=sumvel/num;
factor=sqrt(temp/curtemp);
xvvel=vel*factor;
end
function xvvel=radinrange(rangex,rangey,num)
xloc=randperm(num);
yloc=randperm(num);
vxpro=rand(num,1);
vypro=rand(num,1);
vxr1=rangex(:,2)-rangex(:,1);
vxr2=rangex(:,4)-rangex(:,3);
vxrtot=vxr2+vxr1;
dvx=vxrtot.*vxpro;
ifrx=(dvx<=vxr1);
vx=ifrx.*(dvx+rangex(:,1))+(ifrx==0).*(dvx-vxr1+rangex(:,3));
vyr1=rangey(:,2)-rangey(:,1);
vyr2=rangey(:,4)-rangey(:,3);
vyrtot=vyr2+vyr1;
dvy=vyrtot.*vxpro;
ifry=(dvy<=vyr1);
vy=ifry.*(dvy+rangey(:,1))+(ifry==0).*(dvy-vyr1+rangey(:,3));
xvvel=zeros(num,2);
for i=1:num
    xvvel(i,1)=vx(xloc(i));
    xvvel(i,2)=vy(yloc(i));    
end
end
