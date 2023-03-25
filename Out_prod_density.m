% Outer product of two univariate Gaussian
clear;
close all;

%----------Gamma Mixture-----------%
w_m1=[3/5 2/5]; % weights of the 1st mixture density
Mu1=[2 2];
Sig1=[3 5];
N1=1e4; % particle number
X1=[];  % save particle value
for i=1:N1
    u=rand;
    r=w_m1(1);
    for j=1:length(w_m1);
        if u<r
            if j==1 
                X1=[X1 gamrnd(Mu1(j),Sig1(j))-10];
            else
                X1=[X1 10-gamrnd(Mu1(j),Sig1(j))];
            end
            break;
        else
            r=r+w_m1(j+1);
        end
    end
end
[f1,xi1] = ksdensity(X1); 
  

%---------------Skewed Normal Mixture-------------%
w_m2=[3/4 1/4]; % weights of the 1st mixture density
loc2=[3 -3];
scale2=[1 3];
shape2=[5 -6];
N2=1e4; % particle number
X2=[];  % save particle value
for i=1:N2
    u=rand;
    r=w_m2(1);
    for j=1:length(w_m2);
        if u<r             
            X2=[X2 rsn(1,loc2(j),scale2(j),shape2(j))]; %   rskt(1,loc2(j),scale2(j),shape2(j))];            
            break;
        else
            r=r+w_m2(j+1);
        end
    end
end
[f2,xi2] = ksdensity(X2); 


%----------------------Student's T---------------------%
N3=1e4; % particle number
X3=[];  % save particle value
for i=1:N3
    X3=[X3 mvtrnd(9,4)];
end
[f3,xi3] = ksdensity(X3); 


%---------------Beta+Normal Mixture-------------%
w_m4=[1/2 1/2];
N4=1e4; % particle number
X4=[];  % save particle value
for i=1:N4
    u=rand;
    r=w_m4(1);
    for j=1:length(w_m4);
        if u<r   
            if j==1
               X4=[X4 betarnd(3,3)-3];
            else
               X4=[X4 randn]; 
            end
            break;
        else
            r=r+w_m4(j+1);
        end
    end
end
[f4,xi4] = ksdensity(X4); 
 

%---------------Expontential Mixture-------------%
w_m5=[1/2 1/2];
N5=1e4; % particle number
X5=[];  % save particle value
for i=1:N5
    u=rand;
    r=w_m5(1);
    for j=1:length(w_m5);
        if u<r   
            if j==1
               X5=[X5 exprnd(1)];
            else
               X5=[X5 -exprnd(1)]; 
            end
            break;
        else
            r=r+w_m5(j+1);
        end
    end
end
[f5,xi5] = ksdensity(X5); 
 

%----------------------Skewed Normal---------------------%
N6=1e4; % particle number
X6=[];  % save particle value
for i=1:N6
    X6=[X6 rsn(1,0,8,-3)];
end
[f6,xi6] = ksdensity(X6); 


%----------Normal Mixture-----------%
w_m7=[1/8 1/4 5/8]; % weights of the 1st mixture density
Mu7=[-10 0 7];
Sig7=[.1 .15 .2];
N7=1e4; % particle number
X7=[];  % save particle value
for i=1:N7
    u=rand;
    r=w_m7(1);
    for j=1:length(w_m7);
        if u<r
            X7=[X7 mvnrnd(Mu7(j),Sig7(j))];
            break;
        else
            r=r+w_m7(j+1);
        end
    end
end
[f7,xi7] = ksdensity(X7); 

figure,
subplot(3,3,1);plot(xi1,f1);
subplot(3,3,2);plot(xi2,f2); 
subplot(3,3,3);plot(xi3,f3); 
subplot(3,3,4);plot(xi4,f4);
subplot(3,3,5);plot(xi5,f5);
subplot(3,3,6);plot(xi6,f6); 
subplot(3,3,7);plot(xi7,f7); 
% subplot(7,7,1);plot(xi1,f1);
% subplot(7,7,9);plot(xi2,f2); 
% subplot(7,7,17);plot(xi3,f3); 
% subplot(7,7,25);plot(xi4,f4);
% subplot(7,7,33);plot(xi5,f5);
% subplot(7,7,41);plot(xi6,f6); 
% subplot(7,7,49);plot(xi7,f7); 
