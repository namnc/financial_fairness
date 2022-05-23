%https://tomlankhorst.nl/matlab-to-latex-with-matlab2tikz/
addpath('./matlab2tikz/src/');
addpath('./aboxplot/');

ladderdat=readmatrix('DW_Ladder_4.csv');
lockeddat=readmatrix('DW_LockedLadder_4P_2r.csv');
planteddat=readmatrix('DW_PlantedLadder_4P_2r.csv');
amortizeddat=readmatrix('DW_AmortizedLadder_4P.csv');
bitcoindat=readmatrix('block_gen_mean_sd_variance.csv');

q=10000;
drate = 238;
rrate = -log(1+drate/10000);
myear = 365*24*60;

ladderdim=size(ladderdat);
lockeddim=size(lockeddat);
planteddim=size(planteddat);
amortizeddim=size(amortizeddat);
bitcoindim=size(bitcoindat);

n=ladderdim(1,2)-1;

deposits=zeros(5,n);
%ladderds=zeros(1,n);
%lockedds=zeros(1,n);
%plantedds=zeros(1,n);

for i=1:n
    %ladderds(1,i)=sum(abs(ladderdat(:,i+1)))/2;
    %lockedds(1,i)=sum(abs(lockeddat(:,i+1)))/2;
    %plantedds(1,i)=sum(abs(planteddat(:,i+1)))/2;
    deposits(2,i)=sum(abs(ladderdat(:,i+1)))/2;
    deposits(4,i)=sum(abs(lockeddat(:,i+1)))/2;
    deposits(5,i)=sum(abs(planteddat(:,i+1)))/2;
    deposits(3,i)=sum(abs(amortizeddat(:,i+1)))/2;
    deposits(1,i)= n-1;
end

deposits2 = zeros(5,4);
for i=1:5
    deposits2(i,1) = deposits(i,1);
    deposits2(i,2) = deposits(i,10);
    deposits2(i,3) = deposits(i,25);
    deposits2(i,4) = deposits(i,n);
end

figure
tda = bar(deposits2,'FaceColor','flat');
set(tda(1),'FaceColor','k');
set(tda(2),'FaceColor','g');
set(tda(3),'FaceColor','b');
set(tda(4),'FaceColor','r');
legend('P1','P10','P25','P55');
title('');
ylabel('Total Deposit Amount (x times q)');
xlabel('Penalty Protocols (55 parties, 2 stages)');
xticks([1 2 3 4 5]);
xticklabels({'\MLMech','\LMech','\ALMech','\LLMech','\PLMech'});
%matlab2tikz('TDA.tex');

% figure
% plot(1:n,deposits(5,:),'y',1:n,deposits(1,:),'g',1:n,deposits(2,:),'r',1:n,deposits(3,:),'b');
% ylim([0,max(max(deposits(:,:)))+5]);
% xlim([1, n]);
% legend('MultiLock','Ladder','LockedLadder/AmortizedLadder','PlantedLadder');
% title('Total Deposit Amount');
% ylabel('Total Deposit Amount (x times q)');
% xlabel('Parties');
% xticks([1 2 3 4]);
% xticklabels({'P1','P2','P3','P4'});

% matlab2tikz('TDA.tex');

rounds=max([ladderdim(1,1), lockeddim(1,1), planteddim(1,1)]);
trials=bitcoindim(1,1);

Round=zeros(rounds,trials);

sixblocks = 60;

for t=1:trials
    pdBitcoin = truncate(makedist('Normal','mu',bitcoindat(t,2),'sigma',bitcoindat(t,3)),bitcoindat(t,3),bitcoindat(t,3)*6);
    Round(1,t)= 60;%sixblocks(pdBitcoin);
    for r=2:rounds
        Round(r,t)=Round(r-1,t)+60;%sixblocks(pdBitcoin);
    end
end

RoundBeginDeposits=zeros(5,n);
RoundEndWithdraws=zeros(5,n);
RoundLocks=zeros(5,n);

% Ladder
NPVL=zeros(n,trials);
TimeDepositsL=zeros(n,trials);
for i=1:n
    for r=1:ladderdim(1)
        NPVL(i,:) = NPVL(i,:) - ladderdat(r,i+1).*exp(rrate.*Round(r,:)/myear);
        if (ladderdat(r,i+1) ~= 0)
            if(RoundBeginDeposits(1,i) == 0 || r < RoundBeginDeposits(1,i))
                RoundBeginDeposits(1,i) = r;
            end
            if (r > RoundEndWithdraws(1,i))
                RoundEndWithdraws(1,i) = r;
            end
        end
    end
    TimeDepositsL(i,:)=-Round(RoundBeginDeposits(1,i))+Round(RoundEndWithdraws(1,i));
end

% Locked Ladder
NPVLL=zeros(n,trials);
TimeDepositsLL=zeros(n,trials);
for i=1:n
    for r=1:lockeddim(1)
        NPVLL(i,:) = NPVLL(i,:) - lockeddat(r,i+1).*exp(rrate.*Round(r,:)/myear);
        if (lockeddat(r,i+1) ~= 0)
            if(RoundBeginDeposits(2,i) == 0 || r < RoundBeginDeposits(2,i))
                RoundBeginDeposits(2,i) = r;
            end
            if (r > RoundEndWithdraws(2,i))
                RoundEndWithdraws(2,i) = r;
            end
        end
    end
    TimeDepositsLL(i,:)=-Round(RoundBeginDeposits(2,i))+Round(RoundEndWithdraws(2,i));
end

% Planted Ladder
NPVPL=zeros(n,trials);
TimeDepositsPL=zeros(n,trials);
for i=1:n
    for r=1:planteddim(1)
        NPVPL(i,:) = NPVPL(i,:) - planteddat(r,i+1).*exp(rrate.*Round(r,:)/myear);
        if (planteddat(r,i+1) ~= 0)
            if(RoundBeginDeposits(3,i) == 0 || r < RoundBeginDeposits(3,i))
                RoundBeginDeposits(3,i) = r;
            end
            if (r > RoundEndWithdraws(3,i))
                RoundEndWithdraws(3,i) = r;
            end
        end
    end
    TimeDepositsPL(i,:)=-Round(RoundBeginDeposits(3,i))+Round(RoundEndWithdraws(3,i));
end



% AmortizedLadder
NPVAL=zeros(n,trials);
TimeDepositsAL=zeros(n,trials);
for i=1:n
    for r=1:amortizeddim(1)
        NPVAL(i,:) = NPVAL(i,:) - amortizeddat(r,i+1).*exp(rrate.*Round(r,:)/myear);
        if (amortizeddat(r,i+1) ~= 0)
            if(RoundBeginDeposits(4,i) == 0 || r < RoundBeginDeposits(4,i))
                RoundBeginDeposits(4,i) = r;
            end
            if (r > RoundEndWithdraws(4,i))
                RoundEndWithdraws(4,i) = r;
            end
        end
    end
    TimeDepositsAL(i,:)=-Round(RoundBeginDeposits(4,i))+Round(RoundEndWithdraws(4,i));
end

%Multi Lock
for i=1:n
    RoundBeginDeposits(5,i)=1;
    RoundEndWithdraws(5,i)=2;
end
NPVML=zeros(n,trials);
for i=1:n
    NPVML(i,:) = -(n-1).*exp(rrate.*Round(1,:)/myear) + (n-1).*exp(rrate.*Round(2,:)/myear);
end

RoundLocks=RoundEndWithdraws-RoundBeginDeposits;

% figure
% plot(1:n,RoundLocks(4,:),'y',1:n,RoundLocks(1,:),'g',1:n,RoundLocks(2,:),'r',1:n,RoundLocks(3,:),'b');
% ylim([0,max(max(RoundLocks(:,:)))+5]);
% xlim([1,n]);
% legend('MultiLock/AmortizedLadder','Ladder','LockedLadder','PlantedLadder');
% title('Max Lock Time Window');
% ylabel('Lock Time Window (rounds)');
% xlabel('Parties');
% xticks([1 2 3 4]);
% xticklabels({'P1','P2','P3','P4'});
% % xlabel('');
% 
% matlab2tikz('MLTW.tex');

RoundLocks3 = zeros(5,4);
for i=1:5
    RoundLocks3(i,1) = RoundLocks(i,1);
    RoundLocks3(i,2) = RoundLocks(i,10);
    RoundLocks3(i,3) = RoundLocks(i,25);
    RoundLocks3(i,4) = RoundLocks(i,n);
end

RoundLocks2 = zeros(5,4);
RoundLocks2(1,:) = (RoundLocks3(5,:));
RoundLocks2(2,:) = (RoundLocks3(1,:));
RoundLocks2(3,:) = (RoundLocks3(4,:));
RoundLocks2(4,:) = (RoundLocks3(2,:));
RoundLocks2(5,:) = (RoundLocks3(3,:));



% RoundLocks2 = fliplr(RoundLocks2);

figure
mltw = bar(RoundLocks2,'FaceColor','flat');
set(mltw(1),'FaceColor','k');
set(mltw(2),'FaceColor','g');
set(mltw(3),'FaceColor','b');
set(mltw(4),'FaceColor','r');
legend('P1','P10','P25','P55');
% title('Max Lock Time Window');
ylabel('Lock Time Window (rounds)');
xlabel('Penalty Protocols (55 parties, 2 stages)');
ylim([0,max(max(RoundLocks3))+5]);
xticks([1 2 3 4 5]);
% xticklabels({'\PLMech','\LLMech','\ALMech','\LMech','\MLMech'});
xticklabels({'\MLMech','\LMech','\ALMech','\LLMech','\PLMech'});
% view([90 90]);
%matlab2tikz('MLTW.tex');

NPVB = zeros(5,4);
NPVB(1,1) = NPVML(1,1);
NPVB(1,2) = NPVML(10,1);
NPVB(1,3) = NPVML(25,1);
NPVB(1,4) = NPVML(n,1);
NPVB(2,1) = NPVL(1,1);
NPVB(2,2) = NPVL(10,1);
NPVB(2,3) = NPVL(25,1);
NPVB(2,4) = NPVL(n,1);
NPVB(3,1) = NPVAL(1,1);
NPVB(3,2) = NPVAL(10,1);
NPVB(3,3) = NPVAL(25,1);
NPVB(3,4) = NPVAL(n,1);
NPVB(4,1) = NPVLL(1,1);
NPVB(4,2) = NPVLL(10,1);
NPVB(4,3) = NPVLL(25,1);
NPVB(4,4) = NPVLL(n,1);
NPVB(5,1) = NPVPL(1,1);
NPVB(5,2) = NPVPL(10,1);
NPVB(5,3) = NPVPL(25,1);
NPVB(5,4) = NPVPL(n,1);

figure
ocn = bar(-q*NPVB,'FaceColor','flat');
% ylim([0,max(max(NPVB))+5]);
set(ocn(1),'FaceColor','k');
set(ocn(2),'FaceColor','g');
set(ocn(3),'FaceColor','b');
set(ocn(4),'FaceColor','r');
legend('P1','P10','P25','P55');
title('');
ylabel('Net present cost of participation (bps)');
xlabel('Penalty Protocols (55 parties, 2 stages)');
xticks([1 2 3 4 5]);
xticklabels({'\MLMech','\LMech','\ALMech','\LLMech','\PLMech'});
%matlab2tikz('OCN.tex');

NPVQ = NPVB*q;

% NPVL2 = zeros(4,180);
% NPVAL2 = zeros(4,180);
% NPVML2 = zeros(4,180);
% NPVLL2 = zeros(4,180);
% NPVPL2 = zeros(4,180);
% 
% NPVL2(1,:) = NPVL(1,:);
% NPVL2(2,:) = NPVL(10,:);
% NPVL2(3,:) = NPVL(25,:);
% NPVL2(4,:) = NPVL(n,:);
% 
% NPVAL2(1,:) = NPVAL(1,:);
% NPVAL2(2,:) = NPVAL(10,:);
% NPVAL2(3,:) = NPVAL(25,:);
% NPVAL2(4,:) = NPVAL(n,:);
% 
% NPVML2(1,:) = NPVML(1,:);
% NPVML2(2,:) = NPVML(10,:);
% NPVML2(3,:) = NPVML(25,:);
% NPVML2(4,:) = NPVML(n,:);
% 
% NPVLL2(1,:) = NPVLL(1,:);
% NPVLL2(2,:) = NPVLL(10,:);
% NPVLL2(3,:) = NPVLL(25,:);
% NPVLL2(4,:) = NPVLL(n,:);
% 
% NPVPL2(1,:) = NPVPL(1,:);
% NPVPL2(2,:) = NPVPL(10,:);
% NPVPL2(3,:) = NPVPL(25,:);
% NPVPL2(4,:) = NPVPL(n,:);
% 
% NPVNR=[NPVL2;NPVAL2;NPVML2];
% figure
% boxplot(q.*NPVNR',{reshape(repmat('A':'C',4,1),12,1) repmat((1:4)',3,1)},'factorgap',10,'color','kgbr');
% legend(findobj(gca,'Tag','Box'),'P1','P10','P25','P55');
% xticks([3 8 13]);
% xticklabels({'\LMech','\ALMech','\MLMech'});
% title('');
% ylabel('Net present cost of participation (bps)');
% xlabel('Non-Reactive Penalty Protocols (4 parties)');
% matlab2tikz('OCNR.tex');
% 
% NPVR=[NPVLL2;NPVPL2];
% figure
% boxplot(q.*NPVR',{reshape(repmat('A':'B',4,1),8,1) repmat((1:4)',2,1)},'factorgap',10,'color','kgbr');
% legend(findobj(gca,'Tag','Box'),'P1','P10','P25','P55');
% xticks([2.5 7.5]);
% xticklabels({'\LLMech','\PLMech'});
% title('');
% ylabel('Net present cost of participation (bps)');
% xlabel('Reactive Penalty Protocols (4 parties, 2 stages)');
% matlab2tikz('OCR.tex');

% NPVNR1 = cat(2, NPVL(1,:)', NPVAL(1,:)', NPVML(1,:)');
% NPVNR2 = cat(2, NPVL(2,:)', NPVAL(2,:)', NPVML(2,:)');
% NPVNR3 = cat(2, NPVL(3,:)', NPVAL(3,:)', NPVML(3,:)');
% NPVNR4 = cat(2, NPVL(4,:)', NPVAL(4,:)', NPVML(4,:)');
% NPVNR = {q.*NPVNR1; q.*NPVNR2; q.*NPVNR3; q.*NPVNR4};
% 
% figure
% aboxplot(NPVNR,'labels',{'\LMech','\ALMech','\MLMech'},'Colormap',[1 1 1; 0 1 0 ; 0 0 1 ;1 0 0]); % Advanced box plot
% legend('P1','P2','P3','P4');
% title('');
% ylabel('Net present cost of participation (bps)');
% xlabel('Non-Reactive Penalty Protocols (4 parties)');
% 
% matlab2tikz('OCNR.tex');
% 
% NPVR1 = cat(2, NPVLL(1,:)', NPVPL(1,:)');
% NPVR2 = cat(2, NPVLL(2,:)', NPVPL(2,:)');
% NPVR3 = cat(2, NPVLL(3,:)', NPVPL(3,:)');
% NPVR4 = cat(2, NPVLL(4,:)', NPVPL(4,:)');
% NPVR = {q.*NPVR1; q.*NPVR2; q.*NPVR3; q.*NPVR4};
% 
% figure
% aboxplot(NPVR,'labels',{'\LLMech','\PLMech'},'Colormap',[1 1 1; 0 1 0 ; 0 0 1 ;1 0 0]); % Advanced box plot
% ylabel('Net present cost of participation (bps)');
% title('');
% xlabel('Reactive Penalty Protocols (4 parties, 2 stages)');
% legend('P1','P2','P3','P4');
% matlab2tikz('OCR.tex');

% figure
% title('Oppotunity costs in base points (MultiLock)');
% boxplot(q.*NPVML','Labels',{'P1','P2','P3','P4'});
% ylim([min(min(q.*NPVML(:,:)))-0.5,0]);
% % xlim([1,n]);
% ylabel('Opportunity Costs % (MultiLock)');
% xlabel('Parties');
% 
% matlab2tikz('OCML.tex');
% 
% figure
% title('Oppotunity costs in base points (Amortized Ladder)');
% boxplot(q.*NPVAL','Labels',{'P1','P2','P3','P4'});
% ylim([min(min(q.*NPVAL(:,:)))-0.5,0]);
% % xlim([1,n]);
% ylabel('Opportunity Costs % (AmortizedLadder)');
% xlabel('Parties');
% 
% matlab2tikz('OCAL.tex');
% 
% figure
% title('Oppotunity costs in base points (Locked Ladder)');
% boxplot(q.*NPVLL','Labels',{'P1','P2','P3','P4'});
% ylim([min(min(q.*NPVLL(:,:)))-5,0]);
% % xlim([1,n]);
% ylabel('Opportunity Costs % (Locked Ladder)');
% xlabel('Parties');
% 
% matlab2tikz('OCLL.tex');
% 
% figure
% title('Oppotunity costs in base points (Planted Ladder)');
% boxplot(q.*NPVPL','Labels',{'P1','P2','P3','P4'});
% ylim([min(min(q.*NPVPL(:,:)))-5,0]);
% % xlim([1,n]);
% ylabel('Opportunity Costs % (Planted Ladder)');
% xlabel('Parties');
% 
% matlab2tikz('OCPL.tex');
% 
% NPVC=zeros(6,trials);
% NPVC(1,:)=(NPVML(1,:));
% NPVC(2,:)=(NPVL(1,:));
% NPVC(3,:)=(NPVAL(1,:));
% NPVC(4,:)=(NPVML(4,:));
% NPVC(5,:)=(NPVL(4,:));
% NPVC(6,:)=(NPVAL(4,:));
% 
% 
% figure
% title('Oppotunity costs in base points');
% boxplot((q.*NPVC'),'Labels',{'ML P1','L P1', 'AL P1', 'ML P4','L P4', 'AL P4'});
% ylim([min(min(q.*NPVC(:,:)))-0.1,0.1]);
% % xlim([1,n]);
% ylabel('Opportunity Costs %');
% xlabel('Protocols x Parties');
% % xlabel(['\MLMech';'\LMech';'\LLMech';'\PLMech']);
% 
% matlab2tikz('OCC.tex');

% NPVCL=zeros(2,trials);
% NPVCL(1,:)=(NPVML(4,:));
% NPVCL(2,:)=(NPVL(4,:));
% NPVCL(3,:)=(NPVLL(4,:));
% NPVCL(4,:)=(NPVPL(4,:));
% NPVCL = log10(q.*NPVCL);

% figure
% title('Oppotunity costs in base points (ALL, P4)');
% boxplot((NPVCL'),'Labels',{'\MLMech','\LMech','\LLMech','\PLMech'});
% ylim([min(min(NPVCL(:,:)))-0.1,0.1]);
% % xlim([1,n]);
% ylabel('Opportunity Costs % (log)');
% % xlabel(['\MLMech';'\LMech';'\LLMech';'\PLMech']);
% 
% NPVCF=zeros(4,trials);
% NPVCF(1,:)=NPVML(1,:);
% NPVCF(2,:)=NPVL(1,:);
% NPVCF(3,:)=NPVLL(1,:);
% NPVCF(4,:)=NPVPL(1,:);
% NPVCF = (q.*NPVCF);
% 
% figure
% title('Oppotunity costs in base points (ALL, P1)');
% boxplot(NPVCF','Labels',{'\MLMech','\LMech','\LLMech','\PLMech'});
% ylim([min(min(NPVCF(:,:)))-0.1,0.1]);
% % xlim([1,n]);
% ylabel('Opportunity Costs % (log)');
% % xlabel(['\MLMech';'\LMech';'\LLMech';'\PLMech']);
%%%%%%%%%%%%%%% OLD FIGS %%%%%%%%%%%%%%%

% figure
% plot(1:trials,TimeDeposit(n+1,:),'g',1:trials,TimeDeposit(1,:),'r',1:trials,TimeDeposit(4,:),'b');
% ylim([0,max(max(TimeDeposit(:,:)))]);
% legend('Our Protocol','Party 1','Party 4');
% title('Time your money is in escrow');
% ylabel('Times in Minutes');
% xlabel('Transactions');
% 
% figure
% boxplot(TimeDeposit');
% ylim([0,max(max(TimeDeposit(:,:)))]);
% ylabel('Times in Minutes');
% xlabel('Parties (5 is ours)');

% figure
% plot(1:trials,q.*NPV(n+1,:),'g',1:trials,q.*NPV(1,:),'r',1:trials,q.*NPV(4,:),'b');
% ylim([min(min(q.*NPV(:,:))),0]);
% legend('Our Protocol','Party 1','Party 4');
% title('Opportunity cost for money in escrow');
% ylabel('Money');
% xlabel('Transactions');
% 
% figure
% title('Oppotunity costs in base points');
% boxplot(q.*NPV');
% ylim([min(min(q.*NPV(:,:))),0]);
% ylabel('Opportunity Costs %');
% xlabel('Parties (5 is ours)');
