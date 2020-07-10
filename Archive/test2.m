clear
clc
close all

    stateP(1,:) = {[],[],[],0,[],[],[],[]};
    stateT(1,:) = {[],[],[],0,[],[],[],[]};
    stateh(1,:) = {[],[],[],0,[],[],[],[]};
    
for ii = 1:1000
    stateP(1,2) = {ii*10};
    stateT(1,3) = {ii};
    stateh(1,8) = {ii*100};
    [stateP] = unFAIR3(stateP,1);
    [stateT] = unFAIR3(stateT,1);
    [stateh] = unFAIR3(stateh,1);
    [~,P(ii,1),T(ii,1),~,~,~,~,h(ii,1)] = stateP{1,:};
    [~,P(ii,2),T(ii,2),~,~,~,~,h(ii,2)] = stateT{1,:};
    [~,P(ii,3),T(ii,3),~,~,~,~,h(ii,3)] = stateh{1,:};
    stateP(1,:) = {[],[],[],0,[],[],[],[]};
    stateT(1,:) = {[],[],[],0,[],[],[],[]};
    stateh(1,:) = {[],[],[],0,[],[],[],[]};
end

% figure
% plot(T(:,1))
% figure
% plot(T(:,2))
% figure
% plot(T(:,3))

figure
plot(P(:,1))
figure
plot(P(:,2))
figure
plot(P(:,3))


figure
plot(h(:,1))
figure
plot(h(:,2))
figure
plot(h(:,3))