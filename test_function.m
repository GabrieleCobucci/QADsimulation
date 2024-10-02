clear all
clc

%Ensemble parameters
m = 8; %Number of states in the ensemble
d = 8; %Physical dimension
r = 5; %QAD

%Physical identity
id = eye(d);

%MUBs
omega = exp(2*pi*i/d);

%MUBs
for l = 0 : d-1
    mub{l+1,1} = id(:,l+1);
end

for j = 0 : d-1
    for l = 0 : d-1
        mub{l+1,j+2} = 0;
        for k = 0 : d-1
            mub{l+1,j+2} = mub{l+1,j+2} + 1/sqrt(d)*omega^(k*l + j*k^2)*mub{k+1,1};
        end
    end
end

%Pure states for the ensemble
for x = 1 : m-1
    psi{x} = mub{x,1}*mub{x,1}';
end

%Unitaries
for y = 1 : 2
    U{y} = 0;
    for k = 1 : d
        U{y} = U{y} + mub{k,y}*mub{k,1}';
    end
end

psi{m} = mub{1,2}*mub{1,2}';

v = QADsimulation(d,r,m,2,U,psi);