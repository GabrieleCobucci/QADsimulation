function v = QADsimulation(d,r,m,nUnit,F,psi)

%-------------------------------------------------------------------------%
%This function computes the simulation of a d-dimensional quantum ensemble
%of states subject to white noise with QAD = r.

%Inputs:
% - d: physical dimension of each state in the ensemble
% - r: Quantum Absolute Dimension for the simulation
% - m: number of states in the ensemble
% - nUnit: number of unitaries to use in the simulation
% - F{y}, with y = 1,...,nUnit: unitaries from the computational basis to
% the basis for the simulation
% - psi{x}, with x = 1,...,m: pure target states for the simulation of
% rho{x} = v*psi{x} + (1-v)/d*id

%Output:
% - v: maximum visibility for a QAD = r simulation of the ensemble
%-------------------------------------------------------------------------%

%Identity
id = eye(d);

% Select subspaces

%Number of simulation boxes
nUnit = 2; %Number of unitaries
nSub = nchoosek(d,r); %Number of subspaces
subsets = nchoosek([1:d],r); %Subspaces

nLambda = nUnit*nSub; %Number of boxes

%Projectors into the subspaces
for y = 1 : nUnit
    for mu = 1 : nSub
        subarr = zeros(1,d);
        subarr(1,subsets(mu,:)) = 1; %Fixing to 1 the interested subset
        proj{y,mu} = F{y}*diag(subarr)*F{y}'; %Projectors
    end
end

%% SDP

%Constraints
C = [];

%Variables
v = sdpvar(1);

%Simulation coefficients
sumq = 0;
for y = 1 : nUnit
    for mu = 1 : nSub
        q{y,mu} = sdpvar(1);
        C = [C, q{y,mu} >= 0];
        sumq = sumq + q{y,mu};
    end
end
C = [C, sumq == 1];

%Simulation states
for x = 1 : m
    for y = 1 : nUnit
        for mu = 1 : nSub
            sigma{x,y,mu} = sdpvar(d,d,'hermitian','complex');
            C = [C, sigma{x,y,mu} >= 0, trace(sigma{x,y,mu}) == q{y,mu}];
        end
    end
end

%Projector constraints
for x = 1 :  m
    for y = 1 : nUnit
        for mu = 1 : nSub
            projcon = proj{y,mu}*sigma{x,y,mu}*proj{y,mu};
            C = [C, projcon == sigma{x,y,mu}];
        end
    end
end

%White noise constraint
for x = 1 : m
    s = 0;
    for y = 1 : nUnit
        for mu = 1 : nSub
            s = s + sigma{x,y,mu};
        end
    end
    C = [C, v*psi{x} + (1-v)/d*id == s];
end

%SolveSDP
disp('Options')
ops=sdpsettings('solver','mosek', 'cachesolvers', 1);
diagnostic=solvesdp(C,-v,ops)

v=double(v)

end