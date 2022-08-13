function plotLesionAfterLearning(btAll2d,plotRange)
%%
altMethod = {'mean','median','geomean'};
% global cenMethod edges_RT edges_HT edges_RelT smo_win;
cenMethod = altMethod{2}; % each subject's central estimate
grandCen = 'mean';
grandErr = 'sem';

edges_RT = 0:0.025:0.6; % Reaction Time (only correct)
edges_HT = 0:0.05:2.5; % Hold Time (all trials)
edges_RelT = 0:0.025:1; % Realease Time (correct + late trials)

smo_win = 8; % smoothdata('gaussian'), 
cmpWin = 3; % cmpWin sessions before/after surgery to compare
grpName = ["Lesion";"Sham"];
DataOut = [];
%% Data processing
btAll2d_use = btAll2d(:,plotRange);
% session by session, trial by trial: 2 packaging method
[SBS,TBT] = packData(btAll2d_use);

sess_pre = unique(SBS.Session(cellfun(@(x) ~isempty(x),strfind(SBS.Group,'Pre'))));
sess_post = unique(SBS.Session(cellfun(@(x) ~isempty(x),strfind(SBS.Group,'Post'))));
% subject * group
TBTsg = [estTBT_3FPs(TBT(cellfun(@(x) ~isempty(x),strfind(TBT.Group,'Lesion-Pre')) & ismember(TBT.Session,sess_pre(end-cmpWin+1:end)),:));...
    estTBT_3FPs(TBT(cellfun(@(x) ~isempty(x),strfind(TBT.Group,'Lesion-Post')) & ismember(TBT.Session,sess_post(1:cmpWin)),:));...
    estTBT_3FPs(TBT(cellfun(@(x) ~isempty(x),strfind(TBT.Group,'Sham-Pre')) & ismember(TBT.Session,sess_pre(end-cmpWin+1:end)),:));...
    estTBT_3FPs(TBT(cellfun(@(x) ~isempty(x),strfind(TBT.Group,'Sham-Post')) & ismember(TBT.Session,sess_post(1:cmpWin)),:))];

% group
TBTg = grpstats(removevars(TBTsg,{'Subject','Task','Type'}),'Group',{grandCen,grandErr});

DataOut.TBTsg = TBTsg;
DataOut.TBTg = TBTg; 

% between subject group
SBSbtw = grpstats(addvars(removevars(SBS,{'Subject','Date','Task','Group','Type'}),...
    erase(SBS.Group,{'-Pre','-Post'}),'NewVariableNames','Group'),...
    {'Group','Session'},{'mean',grandErr});
% subject * group
SBSsg = grpstats(removevars(SBS,{'Session','Date','Task','Type'}),{'Subject','Group'},{grandCen,grandErr});

xedges = struct;
xedges.RT = movmean(edges_RT,2,'Endpoints','discard');
xedges.HT = movmean(edges_HT,2,'Endpoints','discard');
xedges.RelT = movmean(edges_RelT,2,'Endpoints','discard');

DataOut.SBSbtw = SBSbtw;
DataOut.SBSsg = SBSsg;
DataOut.xedges = xedges;


%% Statistics
p = struct;
% performance, pre v.s. post
corLesBefore = TBTsg.Cor(strcmp(TBTsg.Group,'Lesion-Pre'))';
corLesAfter = TBTsg.Cor(strcmp(TBTsg.Group,'Lesion-Post'))';
corShamBefore = TBTsg.Cor(strcmp(TBTsg.Group,'Sham-Pre'))';
corShamAfter = TBTsg.Cor(strcmp(TBTsg.Group,'Sham-Post'))';
p.corLe = signrank(corLesBefore,corLesAfter);
p.corSh = signrank(corShamBefore,corShamAfter);
fprintf('Correct Pre-Post signrank test, Lesion p=%.2f, Sham p=%.2f\n',p.corLe,p.corSh);
[~,p.corLeT] = ttest(corLesBefore,corLesAfter);
[~,p.corShT] = ttest(corShamBefore,corShamAfter);
fprintf('Correct Pre-Post paired-ttest test, Lesion p=%.2f, Sham p=%.2f\n',p.corLeT,p.corShT);

preLesBefore = TBTsg.Pre(strcmp(TBTsg.Group,'Lesion-Pre'))';
preLesAfter = TBTsg.Pre(strcmp(TBTsg.Group,'Lesion-Post'))';
preShamBefore = TBTsg.Pre(strcmp(TBTsg.Group,'Sham-Pre'))';
preShamAfter = TBTsg.Pre(strcmp(TBTsg.Group,'Sham-Post'))';
p.preLe = signrank(preLesBefore,preLesAfter);
p.preSh = signrank(preShamBefore,preShamAfter);
fprintf('Premature Pre-Post signrank test, Lesion p=%.2f, Sham p=%.2f\n',p.preLe,p.preSh);
[~,p.preLeT] = ttest(preLesBefore,preLesAfter);
[~,p.preShT] = ttest(preShamBefore,preShamAfter);
fprintf('Premature Pre-Post paired-ttest test, Lesion p=%.2f, Sham p=%.2f\n',p.preLeT,p.preShT);

lateLesBefore = TBTsg.Late(strcmp(TBTsg.Group,'Lesion-Pre'))';
lateLesAfter = TBTsg.Late(strcmp(TBTsg.Group,'Lesion-Post'))';
lateShamBefore = TBTsg.Late(strcmp(TBTsg.Group,'Sham-Pre'))';
lateShamAfter = TBTsg.Late(strcmp(TBTsg.Group,'Sham-Post'))';
p.lateLe = signrank(lateLesBefore,lateLesAfter);
p.lateSh = signrank(lateShamBefore,lateShamAfter);
fprintf('Late Pre-Post signrank test, Lesion p=%.2f, Sham p=%.2f\n',p.lateLe,p.lateSh);
[~,p.lateLeT] = ttest(lateLesBefore,lateLesAfter);
[~,p.lateShT] = ttest(lateShamBefore,lateShamAfter);
fprintf('Late Pre-Post paired-ttest test, Lesion p=%.2f, Sham p=%.2f\n',p.lateLeT,p.lateShT);

% Late: 3FPs × Pre/Post friedman test in Lesion
lateSLesBefore = TBTsg.Late_S(strcmp(TBTsg.Group,'Lesion-Pre'))';
lateSLesAfter = TBTsg.Late_S(strcmp(TBTsg.Group,'Lesion-Post'))';
lateMLesBefore = TBTsg.Late_M(strcmp(TBTsg.Group,'Lesion-Pre'))';
lateMLesAfter = TBTsg.Late_M(strcmp(TBTsg.Group,'Lesion-Post'))';
lateLLesBefore = TBTsg.Late_L(strcmp(TBTsg.Group,'Lesion-Pre'))';
lateLLesAfter = TBTsg.Late_L(strcmp(TBTsg.Group,'Lesion-Post'))';
late3FPs = [lateSLesBefore',lateSLesAfter';...
    lateMLesBefore',lateMLesAfter';...
    lateLLesBefore',lateLLesAfter'];
[p.late3FPsLe,tb1,stats_late3FPs] = friedman(late3FPs,length(lateSLesBefore),'off');
% c = multcompare(stats_late3FPs);
fprintf('Late 3FPs×Pre/Post friedman test in Lesion, Pre/Post effect p=%.2f\n',p.late3FPsLe);

save('VarsToPlot.mat','TBT','TBTsg','TBTg','SBS','SBSbtw','SBSsg','xedges','p');
%% Plot
load('VarsToPlot.mat');

cTab20 = [0.0901960784313726,0.466666666666667,0.701960784313725;0.682352941176471,0.780392156862745,0.901960784313726;0.960784313725490,0.498039215686275,0.137254901960784;0.988235294117647,0.729411764705882,0.470588235294118;0.152941176470588,0.631372549019608,0.278431372549020;0.611764705882353,0.811764705882353,0.533333333333333;0.843137254901961,0.149019607843137,0.172549019607843;0.964705882352941,0.588235294117647,0.592156862745098;0.564705882352941,0.403921568627451,0.674509803921569;0.768627450980392,0.690196078431373,0.827450980392157;0.549019607843137,0.337254901960784,0.290196078431373;0.768627450980392,0.607843137254902,0.576470588235294;0.847058823529412,0.474509803921569,0.698039215686275;0.956862745098039,0.709803921568628,0.807843137254902;0.501960784313726,0.501960784313726,0.501960784313726;0.780392156862745,0.780392156862745,0.776470588235294;0.737254901960784,0.745098039215686,0.196078431372549;0.854901960784314,0.862745098039216,0.549019607843137;0.113725490196078,0.737254901960784,0.803921568627451;0.627450980392157,0.843137254901961,0.890196078431373];
cRed = cTab20(7,:);
cRed2 = cTab20(8,:);
cGreen = cTab20(5,:);
cGreen2 = cTab20(6,:);
cBlue = cTab20(1,:);
cBlue2 = cTab20(2,:);
cGray = cTab20(15,:);
cGray2 = cTab20(16,:);
cOrange = cTab20(3,:);
cOrange2 = cTab20(4,:);


hf = figure(44); clf(hf,'reset');
set(hf, 'name', 'Lesion effect', 'units', 'centimeters', 'position', [1 1 15.8 24.2],...
    'PaperPositionMode', 'auto','renderer','painter'); % 生科论文要求版面裁掉边距还剩，宽度15.8cm,高度24.2cm

size1 = [4,4*0.618];
size2 = [4*0.618,4*0.618];
size3 = [2,4*0.618];
ys = [1 4.7 8.6 12.5 16.3,20.2]; % yStart
xs = [1 6.1 11.2]; % xStart
xs2 = [1 4.25 7.5 10.75];

% PLOT x:sessions, y:%, color:cor/pre/late, Lesion
ha11 = axes;
set(ha11, 'units', 'centimeters', 'position', [xs(1) ys(1) size1], 'nextplot', 'add','tickDir', 'out',...
    'fontsize',7,'fontname','Arial');
lenPre = length(sess_pre)+1;
lenPost = length(sess_post)-1;

% % plot(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Dark(SBSbtw.Group==grpName(1)),...
% %     'linewidth',1.5,'color','k');
% plot(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Cor(SBSbtw.Group==grpName(1)),...
%     'linewidth',1.5,'color',cGreen);
% plot(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Pre(SBSbtw.Group==grpName(1)),...
%     'linewidth',1.5,'color',cRed);
% plot(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Late(SBSbtw.Group==grpName(1)),...
%     'linewidth',1.5,'color',cGray);
% % errorbar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Dark(SBSbtw.Group==grpName(1)),...
% %     SBSbtw.sem_Dark(SBSbtw.Group==grpName(1)),'linewidth',0.8,'color','k');
% errorbar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Cor(SBSbtw.Group==grpName(1)),...
%     SBSbtw.sem_Cor(SBSbtw.Group==grpName(1)),'linewidth',0.8,'color',cGreen,'CapSize',3);
% errorbar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Pre(SBSbtw.Group==grpName(1)),...
%     SBSbtw.sem_Pre(SBSbtw.Group==grpName(1)),'linewidth',0.8,'color',cRed,'CapSize',3);
% errorbar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Late(SBSbtw.Group==grpName(1)),...
%     SBSbtw.sem_Late(SBSbtw.Group==grpName(1)),'linewidth',0.8,'color',cGray,'CapSize',3);
shadedErrorBar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Cor(SBSbtw.Group==grpName(1)),...
    SBSbtw.sem_Cor(SBSbtw.Group==grpName(1)),...
    'lineProps',{'o-','linewidth',1.5,'color',cGreen,'markerSize',4,'markerFaceColor',cGreen,'markerEdgeColor','none'});
shadedErrorBar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Pre(SBSbtw.Group==grpName(1)),...
    SBSbtw.sem_Pre(SBSbtw.Group==grpName(1)),...
    'lineProps',{'o-','linewidth',1.5,'color',cRed,'markerSize',4,'markerFaceColor',cRed,'markerEdgeColor','none'});
shadedErrorBar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Late(SBSbtw.Group==grpName(1)),...
    SBSbtw.sem_Late(SBSbtw.Group==grpName(1)),...
    'lineProps',{'o-','linewidth',1.5,'color',cGray,'markerSize',4,'markerFaceColor',cGray,'markerEdgeColor','none'});
plot([-0.5,-0.5],[0,1],'k','linewidth',0.6);

le1 = legend({'Correct','Premature','Late'},'Fontsize',7,'units','centimeters','Position',[xs(3)+size2(1)+0.628,ys(2)+0.1,1,1]);% [4.7,2.8,1,1]
le1.ItemTokenSize = [12,22];
le1.Position = le1.Position + [0.025 0.045 0 0];
legend('boxoff');
% set(lemark(1:3:end),'markersize',5);
xlim([-lenPre+0.5,lenPost+0.5]);ylim([0,1]);
set(gca,'xtick',[-7,-5,-3,-1,0,2,4,6,8,10],'xticklabel',{'-7','-5','-3','-1','1','3','5','7','9','11'});
xlabel('Sessions','Fontsize',8,'FontName','Arial');
ylabel('Probability','Fontsize',8,'FontName','Arial');
title(grpName(1),'Fontsize',9,'FontName','Arial');

% PLOT x:sessions, y:%, color:cor/pre/late, Sham
ha12 = axes;
set(ha12, 'units', 'centimeters', 'position', [xs(2) ys(1) size1], 'nextplot', 'add','tickDir', 'out',...
    'fontsize',7,'fontname','Arial');
% % plot(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Dark(SBSbtw.Group==grpName(2)),...
% %     'linewidth',1.5,'color','k');
% plot(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Cor(SBSbtw.Group==grpName(2)),...
%     '-','linewidth',1.5,'color',cGreen);
% plot(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Pre(SBSbtw.Group==grpName(2)),...
%     '-','linewidth',1.5,'color',cRed);
% plot(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Late(SBSbtw.Group==grpName(2)),...
%     '-','linewidth',1.5,'color',cGray);
% % errorbar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Dark(SBSbtw.Group==grpName(2)),...
% %     SBSbtw.sem_Dark(SBSbtw.Group==grpName(2)),'linewidth',0.8,'color','k');
% errorbar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Cor(SBSbtw.Group==grpName(2)),...
%     SBSbtw.sem_Cor(SBSbtw.Group==grpName(2)),'linewidth',0.8,'color',cGreen,'CapSize',3);
% errorbar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Pre(SBSbtw.Group==grpName(2)),...
%     SBSbtw.sem_Pre(SBSbtw.Group==grpName(2)),'linewidth',0.8,'color',cRed,'CapSize',3);
% errorbar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Late(SBSbtw.Group==grpName(2)),...
%     SBSbtw.sem_Late(SBSbtw.Group==grpName(2)),'linewidth',0.8,'color',cGray,'CapSize',3);
shadedErrorBar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Cor(SBSbtw.Group==grpName(2)),...
    SBSbtw.sem_Cor(SBSbtw.Group==grpName(2)),...
    'lineProps',{'o-','linewidth',1.5,'color',cGreen,'markerSize',4,'markerFaceColor',cGreen,'markerEdgeColor','none'});
shadedErrorBar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Pre(SBSbtw.Group==grpName(2)),...
    SBSbtw.sem_Pre(SBSbtw.Group==grpName(2)),...
    'lineProps',{'o-','linewidth',1.5,'color',cRed,'markerSize',4,'markerFaceColor',cRed,'markerEdgeColor','none'});
shadedErrorBar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Late(SBSbtw.Group==grpName(2)),...
    SBSbtw.sem_Late(SBSbtw.Group==grpName(2)),...
    'lineProps',{'o-','linewidth',1.5,'color',cGray,'markerSize',4,'markerFaceColor',cGray,'markerEdgeColor','none'});
plot([-0.5,-0.5],[0,1],'k','linewidth',0.6);

xlim([-lenPre+0.5,lenPost+0.5]);ylim([0,1]);
set(gca,'xtick',[-7,-5,-3,-1,0,2,4,6,8,10],'xticklabel',{'-7','-5','-3','-1','1','3','5','7','9','11'});
xlabel('Sessions','Fontsize',8,'FontName','Arial');
ylabel('Probability','Fontsize',8,'FontName','Arial');
title(grpName(2),'Fontsize',9,'FontName','Arial');

% PLOT x:cor/pre/late * Pre/Post * Lesion/Sham, y:%, line&thickness: S/M/L
ha13 = axes;
set(ha13, 'units', 'centimeters', 'position', [xs(3) ys(1) size1], 'nextplot', 'add','tickDir', 'out',...
    'xtick',[],'xticklabel',{},'xticklabelRotation',-45,'fontsize',7,'fontname','Arial');
% 'xtick',[1.25,3.25,5.75,7.75,10.25,12.25]
plot([0.5 1.5],[TBTg.mean_Cor_S(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Cor_S(TBTg.Group==grpName(1)+"-"+"Post")],'o:','lineWidth',2,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');
plot([0.5 1.5],[TBTg.mean_Cor_M(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Cor_M(TBTg.Group==grpName(1)+"-"+"Post")],'o--','lineWidth',1.7,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');
plot([0.5 1.5],[TBTg.mean_Cor_L(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Cor_L(TBTg.Group==grpName(1)+"-"+"Post")],'o-','lineWidth',1.5,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');
plot([2.5 3.5],[TBTg.mean_Cor_S(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Cor_S(TBTg.Group==grpName(2)+"-"+"Post")],'o:','lineWidth',2,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');
plot([2.5 3.5],[TBTg.mean_Cor_M(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Cor_M(TBTg.Group==grpName(2)+"-"+"Post")],'o--','lineWidth',1.7,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');
plot([2.5 3.5],[TBTg.mean_Cor_L(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Cor_L(TBTg.Group==grpName(2)+"-"+"Post")],'o-','lineWidth',1.5,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');
plot([4.5 5.5],[TBTg.mean_Pre_S(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Pre_S(TBTg.Group==grpName(1)+"-"+"Post")],'o:','lineWidth',2,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');
plot([4.5 5.5],[TBTg.mean_Pre_M(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Pre_M(TBTg.Group==grpName(1)+"-"+"Post")],'o--','lineWidth',1.7,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');
plot([4.5 5.5],[TBTg.mean_Pre_L(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Pre_L(TBTg.Group==grpName(1)+"-"+"Post")],'o-','lineWidth',1.5,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');
plot([6.5 7.5],[TBTg.mean_Pre_S(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Pre_S(TBTg.Group==grpName(2)+"-"+"Post")],'o:','lineWidth',2,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');
plot([6.5 7.5],[TBTg.mean_Pre_M(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Pre_M(TBTg.Group==grpName(2)+"-"+"Post")],'o--','lineWidth',1.7,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');
plot([6.5 7.5],[TBTg.mean_Pre_L(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Pre_L(TBTg.Group==grpName(2)+"-"+"Post")],'o-','lineWidth',1.5,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');
plot([8.5 9.5],[TBTg.mean_Late_S(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Late_S(TBTg.Group==grpName(1)+"-"+"Post")],'o:','lineWidth',2,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');
plot([8.5 9.5],[TBTg.mean_Late_M(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Late_M(TBTg.Group==grpName(1)+"-"+"Post")],'o--','lineWidth',1.7,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');
plot([8.5 9.5],[TBTg.mean_Late_L(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Late_L(TBTg.Group==grpName(1)+"-"+"Post")],'o-','lineWidth',1.5,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');
plot([10.5 11.5],[TBTg.mean_Late_S(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Late_S(TBTg.Group==grpName(2)+"-"+"Post")],'o:','lineWidth',2,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');
plot([10.5 11.5],[TBTg.mean_Late_M(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Late_M(TBTg.Group==grpName(2)+"-"+"Post")],'o--','lineWidth',1.7,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');
plot([10.5 11.5],[TBTg.mean_Late_L(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Late_L(TBTg.Group==grpName(2)+"-"+"Post")],'o-','lineWidth',1.5,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');

text([0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5],repelem(-0.02,12),...
    {'Pr','Po','Pr','Po','Pr','Po','Pr','Po','Pr','Po','Pr','Po'},...
    'HorizontalAlignment','left','fontsize',6,'Rotation',-90);
grpAbbr = char(grpName);grpAbbr = cellstr(grpAbbr(:,1:2));
text([1,3,5,7,9,11],repelem(-0.16,6),{grpAbbr{1},grpAbbr{2},grpAbbr{1},grpAbbr{2},grpAbbr{1},grpAbbr{2}},'HorizontalAlignment','center','fontsize',6);
text([2,6,10],repelem(-0.25,3),{'Cor','Pre','Late'},'HorizontalAlignment','center','fontsize',7);

xlim([0 12]);ylim([0 1]);
ylabel('Probability','Fontsize',8,'FontName','Arial');
title('Pre vs. Post','Fontsize',9,'FontName','Arial');

% PLOT Cor pre vs. post
ha21 = axes;
set(ha21, 'units', 'centimeters', 'position', [xs2(1) ys(2) size2], 'nextplot', 'add','tickDir', 'out',...
    'fontsize',8,'fontname','Arial', 'ylim',[0.2 1]);
bh21 = bar(1:2,... % categorical({'Lesion','Sham'})
    [TBTg.mean_Cor(strcmp(TBTg.Group,'Lesion-Pre')),TBTg.mean_Cor(strcmp(TBTg.Group,'Lesion-Post'));...
    TBTg.mean_Cor(strcmp(TBTg.Group,'Sham-Pre')),TBTg.mean_Cor(strcmp(TBTg.Group,'Sham-Post'))],...
    'FaceColor','flat','EdgeColor','none');
xtips1 = bh21(1).XEndPoints;
xtips2 = bh21(2).XEndPoints;
ytips1 = bh21(1).YEndPoints;
ytips2 = bh21(2).YEndPoints;
errorbar(xtips1,ytips1,...
    [TBTg.sem_Cor(strcmp(TBTg.Group,'Lesion-Pre')),TBTg.sem_Cor(strcmp(TBTg.Group,'Sham-Pre'))],...
    '.k','capsize',3,'lineWidth',0.7);
errorbar(xtips2,ytips2,...
    [TBTg.sem_Cor(strcmp(TBTg.Group,'Lesion-Post')),TBTg.sem_Cor(strcmp(TBTg.Group,'Sham-Post'))],...
    '.k','capsize',3,'lineWidth',0.7);
ylm = ylim;ysep = (ylm(2)-ylm(1))/15;yline = ylm(2)-ysep*2; ysym = ylm(2)-ysep;
plot([xtips1(1),xtips2(1)],[yline yline],'-k','lineWidth',0.7);
plot([xtips1(2),xtips2(2)],[yline yline],'-k','lineWidth',0.7);
text(1,ysym,pValue2symbol(p.corLe),'fontsize',8,'HorizontalAlignment','center','fontname','Arial');
text(2,ysym,pValue2symbol(p.corSh),'fontsize',8,'HorizontalAlignment','center','fontname','Arial');
bh21(1).FaceColor = cGray;
bh21(2).FaceColor = cOrange;
set(gca,'xtick',[1,2],'xticklabel',{'Lesion','Sham'});
title('Correct','Fontsize',9,'FontName','Arial');

% PLOT premature pre vs. post
ha22 = axes;
set(ha22, 'units', 'centimeters', 'position', [xs2(2) ys(2) size2], 'nextplot', 'add','tickDir', 'out',...
    'fontsize',8,'fontname','Arial', 'ylim',[0 0.2]);
bh22 = bar(1:2,... % categorical({'Lesion','Sham'})
    [TBTg.mean_Pre(strcmp(TBTg.Group,'Lesion-Pre')),TBTg.mean_Pre(strcmp(TBTg.Group,'Lesion-Post'));...
    TBTg.mean_Pre(strcmp(TBTg.Group,'Sham-Pre')),TBTg.mean_Pre(strcmp(TBTg.Group,'Sham-Post'))],...
    'FaceColor','flat','EdgeColor','none');
xtips1 = bh22(1).XEndPoints;
xtips2 = bh22(2).XEndPoints;
ytips1 = bh22(1).YEndPoints;
ytips2 = bh22(2).YEndPoints;
errorbar(xtips1,ytips1,...
    [TBTg.sem_Pre(strcmp(TBTg.Group,'Lesion-Pre')),TBTg.sem_Pre(strcmp(TBTg.Group,'Sham-Pre'))],...
    '.k','capsize',3,'lineWidth',0.7);
errorbar(xtips2,ytips2,...
    [TBTg.sem_Pre(strcmp(TBTg.Group,'Lesion-Post')),TBTg.sem_Pre(strcmp(TBTg.Group,'Sham-Post'))],...
    '.k','capsize',3,'lineWidth',0.7);
ylm = ylim;ysep = (ylm(2)-ylm(1))/15;yline = ylm(2)-ysep*2; ysym = ylm(2)-ysep;
plot([xtips1(1),xtips2(1)],[yline yline],'-k','lineWidth',0.7);
plot([xtips1(2),xtips2(2)],[yline yline],'-k','lineWidth',0.7);
text(1,ysym,pValue2symbol(p.preLe),'fontsize',8,'HorizontalAlignment','center','fontname','Arial');
text(2,ysym,pValue2symbol(p.preSh),'fontsize',8,'HorizontalAlignment','center','fontname','Arial');

bh22(1).FaceColor = cGray;
bh22(2).FaceColor = cOrange;
set(gca,'xtick',[1,2],'xticklabel',{'Lesion','Sham'});
title('Premature','Fontsize',9,'FontName','Arial');

% PLOT late pre vs. post
ha23 = axes;
set(ha23, 'units', 'centimeters', 'position', [xs2(3) ys(2) size2], 'nextplot', 'add','tickDir', 'out',...
    'fontsize',8,'fontname','Arial', 'ylim',[0 0.8]);
bh23 = bar(1:2,... % categorical({'Lesion','Sham'})
    [TBTg.mean_Late(strcmp(TBTg.Group,'Lesion-Pre')),TBTg.mean_Late(strcmp(TBTg.Group,'Lesion-Post'));...
    TBTg.mean_Late(strcmp(TBTg.Group,'Sham-Pre')),TBTg.mean_Late(strcmp(TBTg.Group,'Sham-Post'))],...
    'FaceColor','flat','EdgeColor','none');
xtips1 = bh23(1).XEndPoints;
xtips2 = bh23(2).XEndPoints;
ytips1 = bh23(1).YEndPoints;
ytips2 = bh23(2).YEndPoints;
errorbar(xtips1,ytips1,...
    [TBTg.sem_Late(strcmp(TBTg.Group,'Lesion-Pre')),TBTg.sem_Late(strcmp(TBTg.Group,'Sham-Pre'))],...
    '.k','capsize',3,'lineWidth',0.7);
errorbar(xtips2,ytips2,...
    [TBTg.sem_Late(strcmp(TBTg.Group,'Lesion-Post')),TBTg.sem_Late(strcmp(TBTg.Group,'Sham-Post'))],...
    '.k','capsize',3,'lineWidth',0.7);
ylm = ylim;ysep = (ylm(2)-ylm(1))/15;yline = ylm(2)-ysep*2; ysym = ylm(2)-ysep;
plot([xtips1(1),xtips2(1)],[yline yline],'-k','lineWidth',0.7);
plot([xtips1(2),xtips2(2)],[yline yline],'-k','lineWidth',0.7);
text(1,ysym,pValue2symbol(p.lateLe),'fontsize',8,'HorizontalAlignment','center','fontname','Arial');
text(2,ysym,pValue2symbol(p.lateSh),'fontsize',8,'HorizontalAlignment','center','fontname','Arial');

bh23(1).FaceColor = cGray;
bh23(2).FaceColor = cOrange;
set(gca,'xtick',[1,2],'xticklabel',{'Lesion','Sham'});
title('Late','Fontsize',9,'FontName','Arial');

% PLOT late-3FPs
ha24 = axes;
set(ha24, 'units', 'centimeters', 'position', [xs2(4) ys(2) size2+[0.5,0]], 'nextplot', 'add','tickDir', 'out',...
    'fontsize',8,'fontname','Arial', 'ylim',[0 0.8]);
bh24 = bar(1:3,... % categorical({'S','M','L'})
    [TBTg.mean_Late_S(strcmp(TBTg.Group,'Lesion-Pre')),TBTg.mean_Late_S(strcmp(TBTg.Group,'Lesion-Post'));...
    TBTg.mean_Late_M(strcmp(TBTg.Group,'Lesion-Pre')),TBTg.mean_Late_M(strcmp(TBTg.Group,'Lesion-Post'));...
    TBTg.mean_Late_L(strcmp(TBTg.Group,'Lesion-Pre')),TBTg.mean_Late_L(strcmp(TBTg.Group,'Lesion-Post'))],...
    'FaceColor','flat','EdgeColor','none');
xtips1 = bh24(1).XEndPoints;
xtips2 = bh24(2).XEndPoints;
ytips1 = bh24(1).YEndPoints;
ytips2 = bh24(2).YEndPoints;
errorbar(xtips1,ytips1,...
    [TBTg.sem_Late_S(strcmp(TBTg.Group,'Lesion-Pre')),...
    TBTg.sem_Late_M(strcmp(TBTg.Group,'Lesion-Pre')),...
    TBTg.sem_Late_L(strcmp(TBTg.Group,'Lesion-Pre'))],...
    '.k','capsize',3,'lineWidth',0.7);
errorbar(xtips2,ytips2,...
    [TBTg.sem_Late_S(strcmp(TBTg.Group,'Lesion-Post')),...
    TBTg.sem_Late_M(strcmp(TBTg.Group,'Lesion-Post')),...
    TBTg.sem_Late_L(strcmp(TBTg.Group,'Lesion-Post'))],...
    '.k','capsize',3,'lineWidth',0.7);
ylm = ylim;ysep = (ylm(2)-ylm(1))/15;yline = ylm(2)-ysep*2; ysym = ylm(2)-ysep;
% plot([xtips1(1),xtips2(1)],[yline yline],'-k','lineWidth',0.7);
% plot([xtips1(2),xtips2(2)],[yline yline],'-k','lineWidth',0.7);
% text(1,ysym,pValue2symbol(p.lateLe),'fontsize',8,'HorizontalAlignment','center','fontname','Arial');
% text(2,ysym,pValue2symbol(p.lateSh),'fontsize',8,'HorizontalAlignment','center','fontname','Arial');

bh24(1).FaceColor = cGray;
bh24(2).FaceColor = cOrange;
set(gca,'xtick',[1,2,3],'xticklabel',{'S','M','L'});% ,'XTickLabelRotation',30
title('Late-3FPs in Lesion','Fontsize',9,'FontName','Arial');


% PLOT LesionGroup, RT distribution, legend: S/M/L * Pre/Post
ha31 = axes;
set(ha31, 'units', 'centimeters', 'position', [xs(1) ys(3) size1], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0,0.6],'xtick',[0 0.2 0.4 0.6],'xticklabel',{'0','200','400','600'},'fontsize',7,'fontname','Arial',...
    'ylim',[0 0.17]);

% plot(xedges.RT,TBTg.mean_RTdist_S(TBTg.Group==grpName(1)+"-"+"Pre",:),':','lineWidth',2,'color',cGray);
% plot(xedges.RT,TBTg.mean_RTdist_M(TBTg.Group==grpName(1)+"-"+"Pre",:),'--','lineWidth',1.7,'color',cGray);
% plot(xedges.RT,TBTg.mean_RTdist_L(TBTg.Group==grpName(1)+"-"+"Pre",:),'-','lineWidth',1.5,'color',cGray);
% plot(xedges.RT,TBTg.mean_RTdist_S(TBTg.Group==grpName(1)+"-"+"Post",:),':','lineWidth',2,'color',cOrange);
% plot(xedges.RT,TBTg.mean_RTdist_M(TBTg.Group==grpName(1)+"-"+"Post",:),'--','lineWidth',1.7,'color',cOrange);
% plot(xedges.RT,TBTg.mean_RTdist_L(TBTg.Group==grpName(1)+"-"+"Post",:),'-','lineWidth',1.5,'color',cOrange);
shadedErrorBar(xedges.RT,TBTg.mean_RTdist_S(TBTg.Group==grpName(1)+"-"+"Pre",:),...
    TBTg.sem_RTdist_S(TBTg.Group==grpName(1)+"-"+"Pre",:),...
    'lineProps',{':','lineWidth',2,'color',cGray},'patchSaturation',0.1);
shadedErrorBar(xedges.RT,TBTg.mean_RTdist_M(TBTg.Group==grpName(1)+"-"+"Pre",:),...
    TBTg.sem_RTdist_M(TBTg.Group==grpName(1)+"-"+"Pre",:),...
    'lineProps',{'--','lineWidth',1.7,'color',cGray},'patchSaturation',0.1);
shadedErrorBar(xedges.RT,TBTg.mean_RTdist_L(TBTg.Group==grpName(1)+"-"+"Pre",:),...
    TBTg.sem_RTdist_L(TBTg.Group==grpName(1)+"-"+"Pre",:),...
     'lineProps',{'-','lineWidth',1.5,'color',cGray},'patchSaturation',0.1);
shadedErrorBar(xedges.RT,TBTg.mean_RTdist_S(TBTg.Group==grpName(1)+"-"+"Post",:),...
    TBTg.sem_RTdist_S(TBTg.Group==grpName(1)+"-"+"Post",:),...
    'lineProps',{':','lineWidth',2,'color',cOrange},'patchSaturation',0.1);
shadedErrorBar(xedges.RT,TBTg.mean_RTdist_M(TBTg.Group==grpName(1)+"-"+"Post",:),...
    TBTg.sem_RTdist_M(TBTg.Group==grpName(1)+"-"+"Post",:),...
    'lineProps',{'--','lineWidth',1.7,'color',cOrange},'patchSaturation',0.1);
shadedErrorBar(xedges.RT,TBTg.mean_RTdist_L(TBTg.Group==grpName(1)+"-"+"Post",:),...
    TBTg.sem_RTdist_L(TBTg.Group==grpName(1)+"-"+"Post",:),...
    'lineProps',{'-','lineWidth',1.5,'color',cOrange},'patchSaturation',0.1);

le4 = legend({'Short','Medium','Long'},...
    'Fontsize',7,'units','centimeters','Position',[xs(3)+size2(1)+0.59,ys(3)+0.1,1,1]); % [10.7,8.9,1,1]
legend('boxoff');
le4.ItemTokenSize = [12,22];
le4.Position = le4.Position + [0.025 0.045 0 0];

xlabel('Reaction time (ms)','Fontsize',8,'FontName','Arial');
ylabel('Probability','Fontsize',8,'FontName','Arial');
title(grpName(1),'Fontsize',9,'FontName','Arial');

% PLOT ShamGroup, RT distribution, legend: S/M/L * Pre/Post
ha32 = axes;
set(ha32, 'units', 'centimeters', 'position', [xs(2) ys(3) size1], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0,0.6],'xtick',[0 0.2 0.4 0.6],'xticklabel',{'0','200','400','600'},'fontsize',7,'fontname','Arial',...
    'ylim',[0 0.17]);

% plot(xedges.RT,TBTg.mean_RTdist_S(TBTg.Group==grpName(2)+"-"+"Pre",:),':','lineWidth',2,'color',cGray);
% plot(xedges.RT,TBTg.mean_RTdist_M(TBTg.Group==grpName(2)+"-"+"Pre",:),'--','lineWidth',1.7,'color',cGray);
% plot(xedges.RT,TBTg.mean_RTdist_L(TBTg.Group==grpName(2)+"-"+"Pre",:),'-','lineWidth',1.5,'color',cGray);
% plot(xedges.RT,TBTg.mean_RTdist_S(TBTg.Group==grpName(2)+"-"+"Post",:),':','lineWidth',2,'color',cOrange);
% plot(xedges.RT,TBTg.mean_RTdist_M(TBTg.Group==grpName(2)+"-"+"Post",:),'--','lineWidth',1.7,'color',cOrange);
% plot(xedges.RT,TBTg.mean_RTdist_L(TBTg.Group==grpName(2)+"-"+"Post",:),'-','lineWidth',1.5,'color',cOrange);

shadedErrorBar(xedges.RT,TBTg.mean_RTdist_S(TBTg.Group==grpName(2)+"-"+"Pre",:),...
    TBTg.sem_RTdist_S(TBTg.Group==grpName(2)+"-"+"Pre",:),...
    'lineProps',{':','lineWidth',2,'color',cGray},'patchSaturation',0.1);
shadedErrorBar(xedges.RT,TBTg.mean_RTdist_M(TBTg.Group==grpName(2)+"-"+"Pre",:),...
    TBTg.sem_RTdist_M(TBTg.Group==grpName(2)+"-"+"Pre",:),...
    'lineProps',{'--','lineWidth',1.7,'color',cGray},'patchSaturation',0.1);
shadedErrorBar(xedges.RT,TBTg.mean_RTdist_L(TBTg.Group==grpName(2)+"-"+"Pre",:),...
    TBTg.sem_RTdist_L(TBTg.Group==grpName(2)+"-"+"Pre",:),...
     'lineProps',{'-','lineWidth',1.5,'color',cGray},'patchSaturation',0.1);
shadedErrorBar(xedges.RT,TBTg.mean_RTdist_S(TBTg.Group==grpName(2)+"-"+"Post",:),...
    TBTg.sem_RTdist_S(TBTg.Group==grpName(2)+"-"+"Post",:),...
    'lineProps',{':','lineWidth',2,'color',cOrange},'patchSaturation',0.1);
shadedErrorBar(xedges.RT,TBTg.mean_RTdist_M(TBTg.Group==grpName(2)+"-"+"Post",:),...
    TBTg.sem_RTdist_M(TBTg.Group==grpName(2)+"-"+"Post",:),...
    'lineProps',{'--','lineWidth',1.7,'color',cOrange},'patchSaturation',0.1);
shadedErrorBar(xedges.RT,TBTg.mean_RTdist_L(TBTg.Group==grpName(2)+"-"+"Post",:),...
    TBTg.sem_RTdist_L(TBTg.Group==grpName(2)+"-"+"Post",:),...
    'lineProps',{'-','lineWidth',1.5,'color',cOrange},'patchSaturation',0.1);

xlabel('Reaction time (ms)','Fontsize',8,'FontName','Arial');
ylabel('Probability','Fontsize',8,'FontName','Arial');
title(grpName(2),'Fontsize',9,'FontName','Arial');

% PLOT x:RT-Pre, y:RT-Post, marker:Lesion/Sham, each point a subjects
ha33 = axes;
set(ha33, 'units', 'centimeters', 'position', [xs(3) ys(3) size2], 'nextplot', 'add','tickDir', 'out',...
   'fontsize',7,'fontname','Arial','xlim',[0.17 0.43],'ylim',[0.17 0.43],'xtick',[0.2,0.3,0.4],'xticklabel',{'200','300','400'},...
   'ytick',[0.2,0.3,0.4],'yticklabel',{'200','300','400'}); %

plot([0,0.6],[0,0.6],':k','lineWidth',0.6)
% errorbar comes from different sessions of one subject
% h61 = errorbar(SBSsg.mean_RT(SBSsg.Group==grpName(1)+"-"+"Pre"),SBSsg.mean_RT(SBSsg.Group==grpName(1)+"-"+"Post"),...
%     SBSsg.sem_RT(SBSsg.Group==grpName(1)+"-"+"Post"),SBSsg.sem_RT(SBSsg.Group==grpName(1)+"-"+"Post"),...
%     SBSsg.sem_RT(SBSsg.Group==grpName(1)+"-"+"Pre"),SBSsg.sem_RT(SBSsg.Group==grpName(1)+"-"+"Pre"),...
%     '.','MarkerSize',12,'MarkerEdgeColor',cBlue,'color',cBlue,'lineWidth',1,'CapSize',3);
% h62 = errorbar(SBSsg.mean_RT(SBSsg.Group==grpName(2)+"-"+"Pre"),SBSsg.mean_RT(SBSsg.Group==grpName(2)+"-"+"Post"),...
%     SBSsg.sem_RT(SBSsg.Group==grpName(2)+"-"+"Post"),SBSsg.sem_RT(SBSsg.Group==grpName(2)+"-"+"Post"),...
%     SBSsg.sem_RT(SBSsg.Group==grpName(2)+"-"+"Pre"),SBSsg.sem_RT(SBSsg.Group==grpName(2)+"-"+"Pre"),...
%     '.','MarkerSize',12,'MarkerEdgeColor',cGray,'color',cGray,'lineWidth',1,'CapSize',3);
RT_CI1 = TBTsg.RT_CI(:,1); RT_CI2 = TBTsg.RT_CI(:,2); % errorbar comes from different trials of one subject
h61 = errorbar(TBTsg.RT(TBTsg.Group==grpName(1)+"-"+"Pre"),TBTsg.RT(TBTsg.Group==grpName(1)+"-"+"Post"),...
    RT_CI1(TBTsg.Group==grpName(1)+"-"+"Post"),RT_CI2(TBTsg.Group==grpName(1)+"-"+"Post"),...
    RT_CI1(TBTsg.Group==grpName(1)+"-"+"Pre"),RT_CI2(TBTsg.Group==grpName(1)+"-"+"Pre"),...
    '.','MarkerSize',12,'MarkerEdgeColor',cBlue,'color',cBlue,'lineWidth',1,'CapSize',3);
h62 = errorbar(TBTsg.RT(TBTsg.Group==grpName(2)+"-"+"Pre"),TBTsg.RT(TBTsg.Group==grpName(2)+"-"+"Post"),...
    RT_CI1(TBTsg.Group==grpName(2)+"-"+"Post"),RT_CI2(TBTsg.Group==grpName(2)+"-"+"Post"),...
    RT_CI1(TBTsg.Group==grpName(2)+"-"+"Pre"),RT_CI2(TBTsg.Group==grpName(2)+"-"+"Pre"),...
    '.','MarkerSize',12,'MarkerEdgeColor',cGray,'color',cGray,'lineWidth',1,'CapSize',3);

le6 = legend([h61,h62],{grpName(1),grpName(2)},'Fontsize',7,'units','centimeters','Position',[xs(3)+size2(1)+0.36,ys(2)+1.3,1,1]); % [14,6.4,1,1]
legend('boxoff');
le6.ItemTokenSize = [15,25];
le6.Position = le6.Position + [0.025 0.045 0 0];

xlabel('Pre RT (ms)','Fontsize',8,'FontName','Arial');
ylabel('Post RT (ms)','Fontsize',8,'FontName','Arial');

% PLOT LesionGroup, HT distribution, legend: S/M/L * Pre/Post
ha41 = axes;
set(ha41, 'units', 'centimeters', 'position', [xs(1) ys(4) size1], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0,2.5],'xtick',[0 0.5 1 1.5 2 2.5],'xticklabel',{'0','500','1000','1500','2000','2500'},...
    'fontsize',7,'fontname','Arial','ylim',[0 0.2]);
fill([0.5,1.1,1.1,0.5],[0,0,1,1],[0.8,0.8,0.8],'EdgeColor','none','FaceAlpha',0.3);
fill([1.0,1.6,1.6,1.0],[0,0,1,1],[0.8,0.8,0.8],'EdgeColor','none','FaceAlpha',0.4);
fill([1.5,2.1,2.1,1.5],[0,0,1,1],[0.8,0.8,0.8],'EdgeColor','none','FaceAlpha',0.5);
% plot(xedges.HT,TBTg.mean_HTdist_S(TBTg.Group==grpName(1)+"-"+"Pre",:),':','lineWidth',2,'color',cGray);
% plot(xedges.HT,TBTg.mean_HTdist_M(TBTg.Group==grpName(1)+"-"+"Pre",:),'--','lineWidth',1.7,'color',cGray);
% h71 = plot(xedges.HT,TBTg.mean_HTdist_L(TBTg.Group==grpName(1)+"-"+"Pre",:),'-','lineWidth',1.5,'color',cGray);
% plot(xedges.HT,TBTg.mean_HTdist_S(TBTg.Group==grpName(1)+"-"+"Post",:),':','lineWidth',2,'color',cOrange);
% plot(xedges.HT,TBTg.mean_HTdist_M(TBTg.Group==grpName(1)+"-"+"Post",:),'--','lineWidth',1.7,'color',cOrange);
% h72 = plot(xedges.HT,TBTg.mean_HTdist_L(TBTg.Group==grpName(1)+"-"+"Post",:),'-','lineWidth',1.5,'color',cOrange);
shadedErrorBar(xedges.HT,TBTg.mean_HTdist_S(TBTg.Group==grpName(1)+"-"+"Pre",:),...
    TBTg.sem_HTdist_S(TBTg.Group==grpName(1)+"-"+"Pre",:),...
    'lineProps',{':','lineWidth',2,'color',cGray},'patchSaturation',0.1);
shadedErrorBar(xedges.HT,TBTg.mean_HTdist_M(TBTg.Group==grpName(1)+"-"+"Pre",:),...
    TBTg.sem_HTdist_M(TBTg.Group==grpName(1)+"-"+"Pre",:),...
    'lineProps',{'--','lineWidth',1.7,'color',cGray},'patchSaturation',0.1);
h71 = shadedErrorBar(xedges.HT,TBTg.mean_HTdist_L(TBTg.Group==grpName(1)+"-"+"Pre",:),...
    TBTg.sem_HTdist_L(TBTg.Group==grpName(1)+"-"+"Pre",:),...
     'lineProps',{'-','lineWidth',1.5,'color',cGray},'patchSaturation',0.1);
shadedErrorBar(xedges.HT,TBTg.mean_HTdist_S(TBTg.Group==grpName(1)+"-"+"Post",:),...
    TBTg.sem_HTdist_S(TBTg.Group==grpName(1)+"-"+"Post",:),...
    'lineProps',{':','lineWidth',2,'color',cOrange},'patchSaturation',0.1);
shadedErrorBar(xedges.HT,TBTg.mean_HTdist_M(TBTg.Group==grpName(1)+"-"+"Post",:),...
    TBTg.sem_HTdist_M(TBTg.Group==grpName(1)+"-"+"Post",:),...
    'lineProps',{'--','lineWidth',1.7,'color',cOrange},'patchSaturation',0.1);
h72 = shadedErrorBar(xedges.HT,TBTg.mean_HTdist_L(TBTg.Group==grpName(1)+"-"+"Post",:),...
    TBTg.sem_HTdist_L(TBTg.Group==grpName(1)+"-"+"Post",:),...
    'lineProps',{'-','lineWidth',1.5,'color',cOrange},'patchSaturation',0.1);
h71 = h71.mainLine; h72 = h72.mainLine;

le7 = legend([h71 h72],{'Pre','Post'},'fontsize',7,'units','centimeter','Position',[xs(3)+size2(1)+0.39,ys(3)+1.4,1,1]); % [10.5,10.2,1,1]
legend('boxoff');
le7.ItemTokenSize = [12,22];
le7.Position = le7.Position + [0.025 0.045 0 0];

xlabel('Press duration (ms)','Fontsize',8,'FontName','Arial');
ylabel('Probability','Fontsize',8,'FontName','Arial');
title(grpName(1),'Fontsize',9,'FontName','Arial');

% PLOT ShamGroup, HT distribution, legend: S/M/L * Pre/Post
ha42 = axes;
set(ha42, 'units', 'centimeters', 'position', [xs(2) ys(4) size1], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0,2.5],'xtick',[0 0.5 1 1.5 2 2.5],'xticklabel',{'0','500','1000','1500','2000','2500'},...
    'fontsize',7,'fontname','Arial','ylim',[0 0.2]);
fill([0.5,1.1,1.1,0.5],[0,0,1,1],[0.8,0.8,0.8],'EdgeColor','none','FaceAlpha',0.3);
fill([1.0,1.6,1.6,1.0],[0,0,1,1],[0.8,0.8,0.8],'EdgeColor','none','FaceAlpha',0.4);
fill([1.5,2.1,2.1,1.5],[0,0,1,1],[0.8,0.8,0.8],'EdgeColor','none','FaceAlpha',0.5);
% plot(xedges.HT,TBTg.mean_HTdist_S(TBTg.Group==grpName(2)+"-"+"Pre",:),':','lineWidth',2,'color',cGray);
% plot(xedges.HT,TBTg.mean_HTdist_M(TBTg.Group==grpName(2)+"-"+"Pre",:),'--','lineWidth',1.7,'color',cGray);
% plot(xedges.HT,TBTg.mean_HTdist_L(TBTg.Group==grpName(2)+"-"+"Pre",:),'-','lineWidth',1.5,'color',cGray);
% plot(xedges.HT,TBTg.mean_HTdist_S(TBTg.Group==grpName(2)+"-"+"Post",:),':','lineWidth',2,'color',cOrange);
% plot(xedges.HT,TBTg.mean_HTdist_M(TBTg.Group==grpName(2)+"-"+"Post",:),'--','lineWidth',1.7,'color',cOrange);
% plot(xedges.HT,TBTg.mean_HTdist_L(TBTg.Group==grpName(2)+"-"+"Post",:),'-','lineWidth',1.5,'color',cOrange);
shadedErrorBar(xedges.HT,TBTg.mean_HTdist_S(TBTg.Group==grpName(2)+"-"+"Pre",:),...
    TBTg.sem_HTdist_S(TBTg.Group==grpName(2)+"-"+"Pre",:),...
    'lineProps',{':','lineWidth',2,'color',cGray},'patchSaturation',0.1);
shadedErrorBar(xedges.HT,TBTg.mean_HTdist_M(TBTg.Group==grpName(2)+"-"+"Pre",:),...
    TBTg.sem_HTdist_M(TBTg.Group==grpName(2)+"-"+"Pre",:),...
    'lineProps',{'--','lineWidth',1.7,'color',cGray},'patchSaturation',0.1);
shadedErrorBar(xedges.HT,TBTg.mean_HTdist_L(TBTg.Group==grpName(2)+"-"+"Pre",:),...
    TBTg.sem_HTdist_L(TBTg.Group==grpName(2)+"-"+"Pre",:),...
     'lineProps',{'-','lineWidth',1.5,'color',cGray},'patchSaturation',0.1);
shadedErrorBar(xedges.HT,TBTg.mean_HTdist_S(TBTg.Group==grpName(2)+"-"+"Post",:),...
    TBTg.sem_HTdist_S(TBTg.Group==grpName(2)+"-"+"Post",:),...
    'lineProps',{':','lineWidth',2,'color',cOrange},'patchSaturation',0.1);
shadedErrorBar(xedges.HT,TBTg.mean_HTdist_M(TBTg.Group==grpName(2)+"-"+"Post",:),...
    TBTg.sem_HTdist_M(TBTg.Group==grpName(2)+"-"+"Post",:),...
    'lineProps',{'--','lineWidth',1.7,'color',cOrange},'patchSaturation',0.1);
shadedErrorBar(xedges.HT,TBTg.mean_HTdist_L(TBTg.Group==grpName(2)+"-"+"Post",:),...
    TBTg.sem_HTdist_L(TBTg.Group==grpName(2)+"-"+"Post",:),...
    'lineProps',{'-','lineWidth',1.5,'color',cOrange},'patchSaturation',0.1);

xlabel('Press duration (ms)','Fontsize',8,'FontName','Arial');
ylabel('Probability','Fontsize',8,'FontName','Arial');
title(grpName(2),'Fontsize',9,'FontName','Arial');

% PLOT LesionGroup, Reaction time in 3FPs, Pre/Post
ha51 = axes;
set(ha51, 'units', 'centimeters', 'position', [xs(1) ys(5) size1], 'nextplot', 'add','tickDir', 'out',...
    'fontsize',8,'fontname','Arial', 'ylim',[0.2 0.45]);
bh51 = bar(1:3,... % categorical({'S','M','L'})
    [TBTg.mean_RT_S(strcmp(TBTg.Group,'Lesion-Pre')),TBTg.mean_RT_S(strcmp(TBTg.Group,'Lesion-Post'));...
    TBTg.mean_RT_M(strcmp(TBTg.Group,'Lesion-Pre')),TBTg.mean_RT_M(strcmp(TBTg.Group,'Lesion-Post'));...
    TBTg.mean_RT_L(strcmp(TBTg.Group,'Lesion-Pre')),TBTg.mean_RT_L(strcmp(TBTg.Group,'Lesion-Post'))],...
    'FaceColor','flat','EdgeColor','none');
xtips1 = bh51(1).XEndPoints;
xtips2 = bh51(2).XEndPoints;
ytips1 = bh51(1).YEndPoints;
ytips2 = bh51(2).YEndPoints;
errorbar(xtips1,ytips1,...
    [TBTg.sem_RT_S(strcmp(TBTg.Group,'Lesion-Pre')),...
    TBTg.sem_RT_M(strcmp(TBTg.Group,'Lesion-Pre')),...
    TBTg.sem_RT_L(strcmp(TBTg.Group,'Lesion-Pre'))],...
    '.k','capsize',3,'lineWidth',0.7);
errorbar(xtips2,ytips2,...
    [TBTg.sem_RT_S(strcmp(TBTg.Group,'Lesion-Post')),...
    TBTg.sem_RT_M(strcmp(TBTg.Group,'Lesion-Post')),...
    TBTg.sem_RT_L(strcmp(TBTg.Group,'Lesion-Post'))],...
    '.k','capsize',3,'lineWidth',0.7);
ylm = ylim;ysep = (ylm(2)-ylm(1))/15;yline = ylm(2)-ysep*2; ysym = ylm(2)-ysep;
% plot([xtips1(1),xtips2(1)],[yline yline],'-k','lineWidth',0.7);
% plot([xtips1(2),xtips2(2)],[yline yline],'-k','lineWidth',0.7);
% text(1,ysym,pValue2symbol(p.lateLe),'fontsize',8,'HorizontalAlignment','center','fontname','Arial');
% text(2,ysym,pValue2symbol(p.lateSh),'fontsize',8,'HorizontalAlignment','center','fontname','Arial');

bh51(1).FaceColor = cGray;
bh51(2).FaceColor = cOrange;
set(gca,'xtick',[1,2,3],'xticklabel',{'S','M','L'},'ytick',[0,0.2,0.3,0.4],'yticklabel',{'0','200','300','400'});% ,'XTickLabelRotation',30
ylabel('RT (ms)','Fontsize',8,'FontName','Arial');
title('RT-3FPs in Lesion','Fontsize',9,'FontName','Arial');

% PLOT ShamGroup, Reaction time in 3FPs, Pre/Post
ha52 = axes;
set(ha52, 'units', 'centimeters', 'position', [xs(2) ys(5) size1], 'nextplot', 'add','tickDir', 'out',...
    'fontsize',8,'fontname','Arial', 'ylim',[0.2 0.45]);
bh52 = bar(1:3,... % categorical({'S','M','L'})
    [TBTg.mean_RT_S(strcmp(TBTg.Group,'Sham-Pre')),TBTg.mean_RT_S(strcmp(TBTg.Group,'Sham-Post'));...
    TBTg.mean_RT_M(strcmp(TBTg.Group,'Sham-Pre')),TBTg.mean_RT_M(strcmp(TBTg.Group,'Sham-Post'));...
    TBTg.mean_RT_L(strcmp(TBTg.Group,'Sham-Pre')),TBTg.mean_RT_L(strcmp(TBTg.Group,'Sham-Post'))],...
    'FaceColor','flat','EdgeColor','none');
xtips1 = bh52(1).XEndPoints;
xtips2 = bh52(2).XEndPoints;
ytips1 = bh52(1).YEndPoints;
ytips2 = bh52(2).YEndPoints;
errorbar(xtips1,ytips1,...
    [TBTg.sem_RT_S(strcmp(TBTg.Group,'Sham-Pre')),...
    TBTg.sem_RT_M(strcmp(TBTg.Group,'Sham-Pre')),...
    TBTg.sem_RT_L(strcmp(TBTg.Group,'Sham-Pre'))],...
    '.k','capsize',3,'lineWidth',0.7);
errorbar(xtips2,ytips2,...
    [TBTg.sem_RT_S(strcmp(TBTg.Group,'Sham-Post')),...
    TBTg.sem_RT_M(strcmp(TBTg.Group,'Sham-Post')),...
    TBTg.sem_RT_L(strcmp(TBTg.Group,'Sham-Post'))],...
    '.k','capsize',3,'lineWidth',0.7);
ylm = ylim;ysep = (ylm(2)-ylm(1))/15;yline = ylm(2)-ysep*2; ysym = ylm(2)-ysep;
% plot([xtips1(1),xtips2(1)],[yline yline],'-k','lineWidth',0.7);
% plot([xtips1(2),xtips2(2)],[yline yline],'-k','lineWidth',0.7);
% text(1,ysym,pValue2symbol(p.lateLe),'fontsize',8,'HorizontalAlignment','center','fontname','Arial');
% text(2,ysym,pValue2symbol(p.lateSh),'fontsize',8,'HorizontalAlignment','center','fontname','Arial');

bh52(1).FaceColor = cGray;
bh52(2).FaceColor = cOrange;
set(gca,'xtick',[1,2,3],'xticklabel',{'S','M','L'},'ytick',[0,0.2,0.3,0.4],'yticklabel',{'0','200','300','400'});% ,'XTickLabelRotation',30
ylabel('RT (ms)','Fontsize',8,'FontName','Arial');
title('RT-3FPs in Sham','Fontsize',9,'FontName','Arial');

% PLOT LesionGroup, ReleaseTime distribution
ha61 = axes;
set(ha61, 'units', 'centimeters', 'position', [xs(1) ys(6) size1], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0,1],'xtick',[0 0.2 0.4 0.6,0.8,1.0],'xticklabel',{'0','200','400','600','800','1000'},'fontsize',7,'fontname','Arial',...
    'ylim',[0 1]);
fill([0,0.6,0.6,0],[0,0,1,1],[0.8,0.8,0.8],'EdgeColor','none','FaceAlpha',0.4);
shadedErrorBar(xedges.RelT,TBTg.mean_RelTdist_S(TBTg.Group==grpName(1)+"-"+"Pre",:),...
    TBTg.sem_RelTdist_S(TBTg.Group==grpName(1)+"-"+"Pre",:),...
    'lineProps',{':','lineWidth',2,'color',cGray},'patchSaturation',0.1);
shadedErrorBar(xedges.RelT,TBTg.mean_RelTdist_M(TBTg.Group==grpName(1)+"-"+"Pre",:),...
    TBTg.sem_RelTdist_M(TBTg.Group==grpName(1)+"-"+"Pre",:),...
    'lineProps',{'--','lineWidth',1.7,'color',cGray},'patchSaturation',0.1);
shadedErrorBar(xedges.RelT,TBTg.mean_RelTdist_L(TBTg.Group==grpName(1)+"-"+"Pre",:),...
    TBTg.sem_RelTdist_L(TBTg.Group==grpName(1)+"-"+"Pre",:),...
     'lineProps',{'-','lineWidth',1.5,'color',cGray},'patchSaturation',0.1);
shadedErrorBar(xedges.RelT,TBTg.mean_RelTdist_S(TBTg.Group==grpName(1)+"-"+"Post",:),...
    TBTg.sem_RelTdist_S(TBTg.Group==grpName(1)+"-"+"Post",:),...
    'lineProps',{':','lineWidth',2,'color',cOrange},'patchSaturation',0.1);
shadedErrorBar(xedges.RelT,TBTg.mean_RelTdist_M(TBTg.Group==grpName(1)+"-"+"Post",:),...
    TBTg.sem_RelTdist_M(TBTg.Group==grpName(1)+"-"+"Post",:),...
    'lineProps',{'--','lineWidth',1.7,'color',cOrange},'patchSaturation',0.1);
shadedErrorBar(xedges.RelT,TBTg.mean_RelTdist_L(TBTg.Group==grpName(1)+"-"+"Post",:),...
    TBTg.sem_RelTdist_L(TBTg.Group==grpName(1)+"-"+"Post",:),...
    'lineProps',{'-','lineWidth',1.5,'color',cOrange},'patchSaturation',0.1);
xlabel('Release time (ms)','Fontsize',8,'FontName','Arial');
ylabel('CDF','Fontsize',8,'FontName','Arial');
title(grpName(1),'Fontsize',9,'FontName','Arial');

% PLOT ShamGroup, ReleaseTime distribution
ha62 = axes;
set(ha62, 'units', 'centimeters', 'position', [xs(2) ys(6) size1], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0,1],'xtick',[0 0.2 0.4 0.6,0.8,1.0],'xticklabel',{'0','200','400','600','800','1000'},'fontsize',7,'fontname','Arial',...
    'ylim',[0 1]);
fill([0,0.6,0.6,0],[0,0,1,1],[0.8,0.8,0.8],'EdgeColor','none','FaceAlpha',0.4);
shadedErrorBar(xedges.RelT,TBTg.mean_RelTdist_S(TBTg.Group==grpName(2)+"-"+"Pre",:),...
    TBTg.sem_RelTdist_S(TBTg.Group==grpName(2)+"-"+"Pre",:),...
    'lineProps',{':','lineWidth',2,'color',cGray},'patchSaturation',0.1);
shadedErrorBar(xedges.RelT,TBTg.mean_RelTdist_M(TBTg.Group==grpName(2)+"-"+"Pre",:),...
    TBTg.sem_RelTdist_M(TBTg.Group==grpName(2)+"-"+"Pre",:),...
    'lineProps',{'--','lineWidth',1.7,'color',cGray},'patchSaturation',0.1);
shadedErrorBar(xedges.RelT,TBTg.mean_RelTdist_L(TBTg.Group==grpName(2)+"-"+"Pre",:),...
    TBTg.sem_RelTdist_L(TBTg.Group==grpName(2)+"-"+"Pre",:),...
     'lineProps',{'-','lineWidth',1.5,'color',cGray},'patchSaturation',0.1);
shadedErrorBar(xedges.RelT,TBTg.mean_RelTdist_S(TBTg.Group==grpName(2)+"-"+"Post",:),...
    TBTg.sem_RelTdist_S(TBTg.Group==grpName(2)+"-"+"Post",:),...
    'lineProps',{':','lineWidth',2,'color',cOrange},'patchSaturation',0.1);
shadedErrorBar(xedges.RelT,TBTg.mean_RelTdist_M(TBTg.Group==grpName(2)+"-"+"Post",:),...
    TBTg.sem_RelTdist_M(TBTg.Group==grpName(2)+"-"+"Post",:),...
    'lineProps',{'--','lineWidth',1.7,'color',cOrange},'patchSaturation',0.1);
shadedErrorBar(xedges.RelT,TBTg.mean_RelTdist_L(TBTg.Group==grpName(2)+"-"+"Post",:),...
    TBTg.sem_RelTdist_L(TBTg.Group==grpName(2)+"-"+"Post",:),...
    'lineProps',{'-','lineWidth',1.5,'color',cOrange},'patchSaturation',0.1);
xlabel('Release time (ms)','Fontsize',8,'FontName','Arial');
ylabel('CDF','Fontsize',8,'FontName','Arial');
title(grpName(2),'Fontsize',9,'FontName','Arial');

% text
grp1Sbj = unique(SBS(cellfun(@(x) ~isempty(x),strfind(SBS.Group,grpName(1))),:).Subject,'stable');
grp2Sbj = unique(SBS(cellfun(@(x) ~isempty(x),strfind(SBS.Group,grpName(2))),:).Subject,'stable');

haSbj = axes('units', 'centimeters', 'position', [15.8 ys(end) 2.642 size1(2)],'Visible','off');
text(haSbj,0,1,[upper(grpName(1)),grp1Sbj'],'fontsize',6,'VerticalAlignment','top');
text(haSbj,0.5,1,[upper(grpName(2)),grp2Sbj'],'fontsize',6,'VerticalAlignment','top');
%% Save
savename = fullfile(pwd,'LesionAfterLearning_BPOD');
saveas(hf,savename,'fig');
print(hf,'-dpng',savename);
print(hf,'-dpdf',savename,'-bestfit');
% print(hf,'-depsc2',savename);
end

%% Functions
function [SBS,TBT] = packData(btAll2d)
SBS = table; % session by session data
TBT = table; % trial by trial data
for i=1:size(btAll2d,1)
    for j=1:size(btAll2d,2)
        T = btAll2d{i,j};
        SBS = [SBS;estSBS(T,j)];
        
        nrow = size(T,1);
        if nrow>1
            tempT = addvars(T,repelem(j,nrow)','After','Date','NewVariableNames','Session');
            TBT = [TBT;tempT];
        end
    end
end
end

function outT = estSBS(data,session)

global cenMethod;

outT = table;
if isempty(data)
    return;
end

sbj = data.Subject(1);
task = data.Task(1);

typename = unique(data.TrialType);
for i=1:length(typename)
    t = struct;
    t.Subject = sbj;
    t.Group =data.Group(1);
    t.Date = data.Date(1);
    t.Session = session;
    t.Task = task;
    t.Type = typename(i);
    tdata = data(data.TrialType==t.Type,:);

    t.nBlock = length(unique(tdata.BlockNum));
    t.nTrial = length(tdata.iTrial);
    t.Dark   = sum(tdata.DarkTry)./(sum(tdata.DarkTry)+t.nTrial);
    t.Cor  = sum(tdata.Outcome=="Cor")./t.nTrial;
    t.Pre  = sum(tdata.Outcome=="Pre")./t.nTrial;
    t.Late = sum(tdata.Outcome=="Late")./t.nTrial;

    t.Cor_S = sum(tdata.Outcome=="Cor" & abs(tdata.FP-0.5)<1e-4)./sum(abs(tdata.FP-0.5)<1e-4);
    t.Cor_M = sum(tdata.Outcome=="Cor" & abs(tdata.FP-1.0)<1e-4)./sum(abs(tdata.FP-1.0)<1e-4);
    t.Cor_L = sum(tdata.Outcome=="Cor" & abs(tdata.FP-1.5)<1e-4)./sum(abs(tdata.FP-1.5)<1e-4);
    t.Pre_S = sum(tdata.Outcome=="Pre" & abs(tdata.FP-0.5)<1e-4)./sum(abs(tdata.FP-0.5)<1e-4);
    t.Pre_M = sum(tdata.Outcome=="Pre" & abs(tdata.FP-1.0)<1e-4)./sum(abs(tdata.FP-1.0)<1e-4);
    t.Pre_L = sum(tdata.Outcome=="Pre" & abs(tdata.FP-1.5)<1e-4)./sum(abs(tdata.FP-1.5)<1e-4);
    t.Late_S = sum(tdata.Outcome=="Late" & abs(tdata.FP-0.5)<1e-4)./sum(abs(tdata.FP-0.5)<1e-4);
    t.Late_M = sum(tdata.Outcome=="Late" & abs(tdata.FP-1.0)<1e-4)./sum(abs(tdata.FP-1.0)<1e-4);
    t.Late_L = sum(tdata.Outcome=="Late" & abs(tdata.FP-1.5)<1e-4)./sum(abs(tdata.FP-1.5)<1e-4);
    
    t.maxFP = max(tdata.FP);
    t.t2mFP = find(tdata.FP==t.maxFP,1,'first');
    t.minRW = min(tdata.RW);
    t.t2mRW = find(tdata.RW==t.minRW,1,'first');
    
    switch cenMethod
        case 'mean'
            t.HT = mean(rmoutliers(tdata.HT,'mean'),'omitnan');
            RT = rmoutliers(tdata(tdata.Outcome=="Cor",:).RT,'mean');
            t.RT = mean(RT(RT>=0.1),'omitnan');
            t.MT = mean(rmoutliers(tdata(tdata.Outcome=="Cor",:).MT,'mean'),'omitnan');
        case 'median'
            t.HT = median(rmoutliers(tdata.HT,'median'),'omitnan');
            RT = rmoutliers(tdata(tdata.Outcome=="Cor",:).RT,'median');
            t.RT = median(RT(RT>=0.1),'omitnan');
            t.MT = median(rmoutliers(tdata(tdata.Outcome=="Cor",:).MT,'median'),'omitnan');
        case 'geomean'
            t.HT = geomean(rmoutliers(tdata.HT,'quartiles'),'omitnan');
            RT = rmoutliers(tdata(tdata.Outcome=="Cor" & tdata.RT>0,:).RT,'quartiles');
            t.RT = geomean(RT(RT>=0.1),'omitnan');
            t.MT = geomean(rmoutliers(tdata(tdata.Outcome=="Cor",:).MT,'quartiles'),'omitnan');
    end
    outT = [outT;struct2table(t)];
end

end

function outT = estTBT_3FPs(TBT)
global cenMethod edges_RT edges_HT edges_RelT smo_win
fplist = [0.5,1.0,1.5];
nboot = 1000;

outT = table;
sbjlist = unique(TBT.Subject);

for i=1:length(sbjlist)
    data = TBT(TBT.Subject==sbjlist(i),:);
    typename = unique(data.TrialType);
    for j=1:length(typename)
        t = struct;
        t.Subject = sbjlist(i);
        t.Group = data.Group(1);
        t.Task = data.Task(1);
        t.Type = typename(j);
        tdata = data(data.TrialType==t.Type,:);

        t.nSession = length(unique(tdata.Session));
        t.nTrial = size(tdata,1);
        t.Dark = sum(tdata.DarkTry)./(sum(tdata.DarkTry)+t.nTrial);

        idxFPS = abs(tdata.FP-fplist(1))<1E-4; % small
        idxFPM = abs(tdata.FP-fplist(2))<1E-4; % medium
        idxFPL = abs(tdata.FP-fplist(3))<1E-4; % large
        idxCor = tdata.Outcome=="Cor";
        idxPre = tdata.Outcome=="Pre";
        idxLate = tdata.Outcome=="Late";

        t.Cor = sum(idxCor)./t.nTrial;
        t.Pre = sum(idxPre)./t.nTrial;
        t.Late = sum(idxLate)./t.nTrial;

        t.Cor_S = sum( idxFPS & idxCor )./sum(idxFPS);
        t.Pre_S = sum( idxFPS & idxPre )./sum(idxFPS);
        t.Late_S = sum( idxFPS & idxLate )./sum(idxFPS);
        t.Cor_M = sum( idxFPM & idxCor )./sum(idxFPM);
        t.Pre_M = sum( idxFPM & idxPre )./sum(idxFPM);
        t.Late_M = sum( idxFPM & idxLate )./sum(idxFPM);
        t.Cor_L = sum( idxFPL & idxCor )./sum(idxFPL);
        t.Pre_L = sum( idxFPL & idxPre )./sum(idxFPL);
        t.Late_L = sum( idxFPL & idxLate )./sum(idxFPL);

        switch cenMethod
            case 'mean'
                RT = rmoutliers(tdata.RT(idxCor),'mean');
                RT_S = rmoutliers(tdata.RT(idxCor&idxFPS),'mean');
                RT_M = rmoutliers(tdata.RT(idxCor&idxFPM),'mean');
                RT_L = rmoutliers(tdata.RT(idxCor&idxFPL),'mean');
                RT = RT(RT>=0.1);
                RT_S = RT_S(RT_S>=0.1);
                RT_M = RT_M(RT_M>=0.1);
                RT_L = RT_L(RT_L>=0.1);
                t.RT = mean(RT,'omitnan');
                RT_CI = bootci(nboot,{@mean,RT},'alpha',0.05)';
                t.RT_CI = [t.RT-RT_CI(1), RT_CI(2)-t.RT];
                t.RT_S = mean(RT_S,'omitnan');
                t.RT_M = mean(RT_M,'omitnan');
                t.RT_L = mean(RT_L,'omitnan');
            case 'median'
                RT = rmoutliers(tdata.RT(idxCor),'median');
                RT_S = rmoutliers(tdata.RT(idxCor&idxFPS),'median');
                RT_M = rmoutliers(tdata.RT(idxCor&idxFPM),'median');
                RT_L = rmoutliers(tdata.RT(idxCor&idxFPL),'median');
                RT = RT(RT>=0.1);
                RT_S = RT_S(RT_S>=0.1);
                RT_M = RT_M(RT_M>=0.1);
                RT_L = RT_L(RT_L>=0.1);
                t.RT = median(RT,'omitnan');
                RT_CI = bootci(nboot,{@median,RT},'alpha',0.05)';
                t.RT_CI = [t.RT-RT_CI(1), RT_CI(2)-t.RT];
                t.RT_S = median(RT_S,'omitnan');
                t.RT_M = median(RT_M,'omitnan');
                t.RT_L = median(RT_L,'omitnan');
            case 'geomean'
                RT = rmoutliers(tdata.RT(idxCor),'quartiles');
                RT_S = rmoutliers(tdata.RT(idxCor&idxFPS),'quartiles');
                RT_M = rmoutliers(tdata.RT(idxCor&idxFPM),'quartiles');
                RT_L = rmoutliers(tdata.RT(idxCor&idxFPL),'quartiles');
                RT = RT(RT>=0.1);
                RT_S = RT_S(RT_S>=0.1);
                RT_M = RT_M(RT_M>=0.1);
                RT_L = RT_L(RT_L>=0.1);
                t.RT = geomean(RT,'omitnan');
                RT_CI = bootci(nboot,{@geomean,RT},'alpha',0.05)';
                t.RT_CI = [t.RT-RT_CI(1), RT_CI(2)-t.RT];
                t.RT_S = geomean(RT_S,'omitnan');
                t.RT_M = geomean(RT_M,'omitnan');
                t.RT_L = geomean(RT_L,'omitnan');
        end
        t.RTdist = smoothdata(histcounts(tdata.RT(idxCor),...
            edges_RT,'Normalization','probability'),2,'gaussian',smo_win);
        t.RTdist_S = smoothdata(histcounts(tdata.RT(idxCor&idxFPS),...
            edges_RT,'Normalization','probability'),2,'gaussian',smo_win);
        t.RTdist_M = smoothdata(histcounts(tdata.RT(idxCor&idxFPM),...
            edges_RT,'Normalization','probability'),2,'gaussian',smo_win);
        t.RTdist_L = smoothdata(histcounts(tdata.RT(idxCor&idxFPL),...
            edges_RT,'Normalization','probability'),2,'gaussian',smo_win);

        t.HTdist = smoothdata(histcounts(tdata.HT,...
            edges_HT,'Normalization','probability'),2,'gaussian',smo_win);
        t.HTdist_S = smoothdata(histcounts(tdata.HT(idxFPS),...
            edges_HT,'Normalization','probability'),2,'gaussian',smo_win);
        t.HTdist_M = smoothdata(histcounts(tdata.HT(idxFPM),...
            edges_HT,'Normalization','probability'),2,'gaussian',smo_win);
        t.HTdist_L = smoothdata(histcounts(tdata.HT(idxFPL),...
            edges_HT,'Normalization','probability'),2,'gaussian',smo_win);
        
        t.RelTdist = smoothdata(histcounts(tdata.HT(idxCor|idxLate)-tdata.FP(idxCor|idxLate),...
            edges_RelT,'Normalization','cdf'),2,'gaussian',smo_win);
        t.RelTdist_S = smoothdata(histcounts(tdata.HT((idxCor|idxLate)&idxFPS)-tdata.FP((idxCor|idxLate)&idxFPS),...
            edges_RelT,'Normalization','cdf'),2,'gaussian',smo_win);
        t.RelTdist_M = smoothdata(histcounts(tdata.HT((idxCor|idxLate)&idxFPM)-tdata.FP((idxCor|idxLate)&idxFPM),...
            edges_RelT,'Normalization','cdf'),2,'gaussian',smo_win);
        t.RelTdist_L = smoothdata(histcounts(tdata.HT((idxCor|idxLate)&idxFPL)-tdata.FP((idxCor|idxLate)&idxFPL),...
            edges_RelT,'Normalization','cdf'),2,'gaussian',smo_win);

        outT = [outT;struct2table(t)];
    end
end

end

function symbol = pValue2symbol(p)
    if p>0.05
%         symbol = 'n.s';
        symbol = sprintf('%0.2f',p);
    elseif p>0.01 && p<=0.05
        symbol = '*';
    elseif p>0.001 && p<=0.01
        symbol = '**';
    else
        symbol = '***';
    end
end