function plotLesionAfterLearningJ(btAll2d,plotRange)
%%
set_matlab_default
altMethod = {'mean','median','geomean'};
global cenMethod edges_RT edges_HT edges_RelT smo_win;
cenMethod = altMethod{2}; % each subject's central estimate
grandCen = 'mean';
grandErr = 'sem';

edges_RT = 0:0.025:0.6; % Reaction Time (only correct)
edges_HT = 0:0.05:2.5; % Hold Time (all trials)
edges_RelT = 0:0.05:1; % Realease Time (correct + late trials)

smo_win = 8; % smoothdata('gaussian'), 
cmpWin = 3; % cmpWin sessions before/after surgery to compare
grpName = ["Lesion";"Sham"];
fplist = [0.5,1.0,1.5];
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

% Correct, pre v.s. post
corLesBefore         =      TBTsg.Cor(strcmp(TBTsg.Group,'Lesion-Pre'))';
corLesAfter            =      TBTsg.Cor(strcmp(TBTsg.Group,'Lesion-Post'))';
corShamBefore     =      TBTsg.Cor(strcmp(TBTsg.Group,'Sham-Pre'))';
corShamAfter        =      TBTsg.Cor(strcmp(TBTsg.Group,'Sham-Post'))';

p.CorrectLesionBefore = corLesBefore;
p.CorrectLesionAfter = corLesAfter;
p.CorrectShamBefore = corShamBefore;
p.CorrectShamAfter = corShamAfter;

p.CorrectLesionBeforeAfter_SignrankTest = signrank(corLesBefore,corLesAfter);
p.CorrectShamBeforeAfter_SignrankTest = signrank(corShamBefore,corShamAfter);

p.Correct_SignrankTest = sprintf('Correct Pre-Post signrank test, Lesion p=%.5f, Sham p=%.5f', ...
    p.CorrectLesionBeforeAfter_SignrankTest,p.CorrectShamBeforeAfter_SignrankTest);

[~,p.CorrectLesionBeforeAfter_ttest] = ttest(corLesBefore,corLesAfter);
[~,p.CorrectShamBeforeAfter_ttest] = ttest(corShamBefore,corShamAfter);

p.Correct_ttest = sprintf('Correct Pre-Post paired-ttest test, Lesion p=%.5f, Sham p=%.5f', ...
    p.CorrectLesionBeforeAfter_ttest,p.CorrectShamBeforeAfter_ttest);

% Premature, pre v.s. post

preLesBefore = TBTsg.Pre(strcmp(TBTsg.Group,'Lesion-Pre'))';
preLesAfter = TBTsg.Pre(strcmp(TBTsg.Group,'Lesion-Post'))';
preShamBefore = TBTsg.Pre(strcmp(TBTsg.Group,'Sham-Pre'))';
preShamAfter = TBTsg.Pre(strcmp(TBTsg.Group,'Sham-Post'))';

p.PrematureLesionBefore = preLesBefore;
p.PrematureLesionAfter = preLesAfter;
p.PrematureShamBefore = preShamBefore;
p.PrematureShamAfter = preShamAfter;

p.PrematureLesionBeforeAfter_SignrankTest = signrank(preLesBefore,preLesAfter);
p.PrematureShamBeforeAfter_SignrankTest = signrank(preShamBefore,preShamAfter);

p.Premature_SignrankTest = sprintf('Premature Pre-Post signrank test, Lesion p=%.5f, Sham p=%.5f', ...
    p.PrematureLesionBeforeAfter_SignrankTest,p.PrematureShamBeforeAfter_SignrankTest);

[~,p.PrematureBeforeAfter_ttest] = ttest(preLesBefore,preLesAfter);
[~,p.PrematureShamBeforeAfter_ttest] = ttest(preShamBefore,preShamAfter);

p.Premature_ttest = sprintf('Premature Pre-Post paired-ttest test, Lesion p=%.5f, Sham p=%.5f', ...
    p.PrematureBeforeAfter_ttest,p.PrematureShamBeforeAfter_ttest);

% Late, pre v.s. post

lateLesBefore = TBTsg.Late(strcmp(TBTsg.Group,'Lesion-Pre'))';
lateLesAfter = TBTsg.Late(strcmp(TBTsg.Group,'Lesion-Post'))';
lateShamBefore = TBTsg.Late(strcmp(TBTsg.Group,'Sham-Pre'))';
lateShamAfter = TBTsg.Late(strcmp(TBTsg.Group,'Sham-Post'))';

p.LateLesionBefore = lateLesBefore;
p.LateLesionAfter = lateLesAfter;
p.LateShamBefore = lateShamBefore;
p.LateShamAfter = lateShamAfter;

p.LateLesionBeforeAfter_SignrankTest = signrank(lateLesBefore,lateLesAfter);
p.LateShamBeforeAfter_SignrankTest = signrank(lateShamBefore,lateShamAfter);

p.Late_SignrankTest = sprintf('Late Pre-Post signrank test, Lesion p=%.5f, Sham p=%.5f', ...
    p.LateLesionBeforeAfter_SignrankTest,p.LateShamBeforeAfter_SignrankTest);

[~,p.LateBeforeAfter_ttest] = ttest(lateLesBefore,lateLesAfter);
[~,p.LateShamBeforeAfter_ttest] = ttest(lateShamBefore,lateShamAfter);

p.Late_ttest = sprintf('Late Pre-Post paired-ttest test, Lesion p=%.5f, Sham p=%.5f', ...
    p.LateBeforeAfter_ttest,p.LateShamBeforeAfter_ttest);

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
fprintf('Late 3FPs×Pre/Post friedman test in Lesion, Pre/Post effect p=%.2f',p.late3FPsLe);

DataOut.p = p

save('DataOut.mat', 'DataOut')

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
set(hf, 'name', 'Lesion effect', 'units', 'centimeters', 'position', [1 1 17 12],...
    'PaperPositionMode', 'auto','renderer','painter');  
size1 = [3.5,3.5*0.7]; 
size4 = [3.5 3.5*0.7];  
size5 = [3.5 3.5*0.7]; % compare performance
xs = [1.5 6.0 11]; % xStart
ys = [1 9  5.2 13.5 16.3,20.2]; % yStart

xs2 = [1.5 5.5 9.5 13.5];

% PLOT x:sessions, y:%, color:cor/pre/late, Lesion
ha11 = axes;
set(ha11, 'units', 'centimeters', 'position', [xs(1) ys(1) size1], 'nextplot', 'add','tickDir', 'out',...
    'fontsize',7,'ticklength', [0.02 0.025]);
lenPre = length(sess_pre)+1;
lenPost = length(sess_post)-1;

shadedErrorBar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Cor(SBSbtw.Group==grpName(1)),...
    SBSbtw.sem_Cor(SBSbtw.Group==grpName(1)),...
    'lineProps',{'o-','linewidth',1,'color',cGreen,'markerSize',4,'markerFaceColor',cGreen,'markerEdgeColor','none'});
shadedErrorBar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Pre(SBSbtw.Group==grpName(1)),...
    SBSbtw.sem_Pre(SBSbtw.Group==grpName(1)),...
    'lineProps',{'o-','linewidth',1,'color',cRed,'markerSize',4,'markerFaceColor',cRed,'markerEdgeColor','none'});
shadedErrorBar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Late(SBSbtw.Group==grpName(1)),...
    SBSbtw.sem_Late(SBSbtw.Group==grpName(1)),...
    'lineProps',{'o-','linewidth',1,'color',cGray,'markerSize',4,'markerFaceColor',cGray,'markerEdgeColor','none'});
plot([-0.5,-0.5],[0,1],'k','linewidth',0.6);

le1 = legend({'Correct','Premature','Late'}, ...
    'units','centimeters','Position',[xs(3)+size4(1)+0.5,ys(1)+0.1,1,1]); 

le1.ItemTokenSize = [12,22];
le1.Position = le1.Position + [0.025 0.045 0 0];
le1.FontSize = 7;
legend('boxoff');
% set(lemark(1:3:end),'markersize',5);
xlim([-lenPre+0.5,lenPost+0.5]);ylim([0,1]);
 
set(gca,'xtick',[-7,-5,-3,-1,0,2,4,6,8,10],'xticklabel',{'-7','-5','-3','-1','1','3','5','7','9','11'}, ...
    'ytick',[0:0.5:1], 'yticklabel',{'0', '50', '100'});
xlabel('Sessions');
ylabel('Percentage (%)');
title(grpName(1),'Fontsize',8);

% PLOT x:sessions, y:%, color:cor/pre/late, Sham
ha12 = axes;
set(ha12, 'units', 'centimeters', 'position', [xs(2) ys(1) size1], 'nextplot', 'add','tickDir', 'out','ticklength', [0.02 0.025]);
shadedErrorBar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Cor(SBSbtw.Group==grpName(2)),...
    SBSbtw.sem_Cor(SBSbtw.Group==grpName(2)),...
    'lineProps',{'o-','linewidth',1,'color',cGreen,'markerSize',4,'markerFaceColor',cGreen,'markerEdgeColor','none'});
shadedErrorBar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Pre(SBSbtw.Group==grpName(2)),...
    SBSbtw.sem_Pre(SBSbtw.Group==grpName(2)),...
    'lineProps',{'o-','linewidth',1,'color',cRed,'markerSize',4,'markerFaceColor',cRed,'markerEdgeColor','none'});
shadedErrorBar(unique(SBSbtw.Session)-lenPre,SBSbtw.mean_Late(SBSbtw.Group==grpName(2)),...
    SBSbtw.sem_Late(SBSbtw.Group==grpName(2)),...
    'lineProps',{'o-','linewidth',1,'color',cGray,'markerSize',4,'markerFaceColor',cGray,'markerEdgeColor','none'});
plot([-0.5,-0.5],[0,1],'k','linewidth',0.6);

xlim([-lenPre+0.5,lenPost+0.5]);ylim([0,1]);
set(gca,'xtick',[-7,-5,-3,-1,0,2,4,6,8,10],'xticklabel',{'-7','-5','-3','-1','1','3','5','7','9','11'}, ...
    'ytick',[0:0.5:1], 'yticklabel',{'0', '50', '100'});
xlabel('Sessions');
% ylabel('Probability');
title(grpName(2),'Fontsize',8);

%% plot performance for lesion and sham groups. 
% PLOT x:cor/pre/late * Pre/Post * Lesion/Sham, y:%, line&thickness: S/M/L

ha13 = axes;
set(ha13, 'units', 'centimeters', 'position', [xs(3) ys(1) size5], 'nextplot', 'add','tickDir', 'out',...
    'xtick',[],'xticklabel',{},'xticklabelRotation',-45, 'xlim', [0 7.5],'ylim', [0 1], 'ytick',[0:0.5:1], ...
    'yticklabel',{'0', '50', '100'},'ticklength', [0.02 0.025] );
% Lesion group
% Correct (TBTg.Group==grpName(1)) 0.75, 0.25
plot([0.5 1.25],[TBTg.mean_Cor_S(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Cor_S(TBTg.Group==grpName(1)+"-"+"Post")],'.-','lineWidth',0.5,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');
plot([0.5 1.25],[TBTg.mean_Cor_M(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Cor_M(TBTg.Group==grpName(1)+"-"+"Post")],'.-','lineWidth',1,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');
plot([0.5 1.25],[TBTg.mean_Cor_L(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Cor_L(TBTg.Group==grpName(1)+"-"+"Post")],'.-','lineWidth',1.5,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');
% Premature (TBTg.Group==grpName(1))
plot([1.5 2.25],[TBTg.mean_Pre_S(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Pre_S(TBTg.Group==grpName(1)+"-"+"Post")],'.-','lineWidth',0.5,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');
plot([1.5 2.25],[TBTg.mean_Pre_M(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Pre_M(TBTg.Group==grpName(1)+"-"+"Post")],'.-','lineWidth',1,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');
plot([1.5 2.25],[TBTg.mean_Pre_L(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Pre_L(TBTg.Group==grpName(1)+"-"+"Post")],'.-','lineWidth',1.5,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');
% Late (TBTg.Group==grpName(1))
plot([2.5 3.25],[TBTg.mean_Late_S(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Late_S(TBTg.Group==grpName(1)+"-"+"Post")],'.-','lineWidth',0.5,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');
plot([2.5 3.25],[TBTg.mean_Late_M(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Late_M(TBTg.Group==grpName(1)+"-"+"Post")],'.-','lineWidth',1,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');
plot([2.5 3.25],[TBTg.mean_Late_L(TBTg.Group==grpName(1)+"-"+"Pre"),TBTg.mean_Late_L(TBTg.Group==grpName(1)+"-"+"Post")],'.-','lineWidth',1.5,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');

% Sham group
% Correct (TBTg.Group==grpName(2))
plot([4.25 5],[TBTg.mean_Cor_S(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Cor_S(TBTg.Group==grpName(2)+"-"+"Post")],'-','lineWidth',0.5,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');
plot([4.25 5],[TBTg.mean_Cor_M(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Cor_M(TBTg.Group==grpName(2)+"-"+"Post")],'-','lineWidth',1,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');
plot([4.25 5],[TBTg.mean_Cor_L(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Cor_L(TBTg.Group==grpName(2)+"-"+"Post")],'-','lineWidth',1.5,'color',cGreen,'markersize',4,'markerFaceColor',cGreen,'markerEdgeColor','none');

% Premature (TBTg.Group==grpName(2))
plot([5.25 6],[TBTg.mean_Pre_S(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Pre_S(TBTg.Group==grpName(2)+"-"+"Post")],'-','lineWidth',0.5,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');
plot([5.25 6],[TBTg.mean_Pre_M(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Pre_M(TBTg.Group==grpName(2)+"-"+"Post")],'-','lineWidth',1,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');
plot([5.25 6],[TBTg.mean_Pre_L(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Pre_L(TBTg.Group==grpName(2)+"-"+"Post")],'-','lineWidth',1.5,'color',cRed,'markersize',4,'markerFaceColor',cRed,'markerEdgeColor','none');

% Late (TBTg.Group==grpName(2))
plot([6.25 7],[TBTg.mean_Late_S(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Late_S(TBTg.Group==grpName(2)+"-"+"Post")],'-','lineWidth',0.5,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');
plot([6.25 7],[TBTg.mean_Late_M(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Late_M(TBTg.Group==grpName(2)+"-"+"Post")],'-','lineWidth',1,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');
plot([6.25 7],[TBTg.mean_Late_L(TBTg.Group==grpName(2)+"-"+"Pre"),TBTg.mean_Late_L(TBTg.Group==grpName(2)+"-"+"Post")],'-','lineWidth',1.5,'color',cGray,'markersize',4,'markerFaceColor',cGray,'markerEdgeColor','none');

text(1, -0.125, 'Lesion', 'fontsize', 8, 'fontName', 'dejavu sans')
text(4.5, -0.125, 'Sham','fontsize', 8, 'fontName', 'dejavu sans')
ylabel('Percentage (%)');

%% PLOT LesionGroup, RT distribution, legend: S/M/L * Pre/Post
% Release time
ha31 = axes;
set(ha31, 'units', 'centimeters', 'position', [xs(1) ys(2) size4], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0 1],'xtick',[0:0.25:1],'xticklabel',{'0','250', '500','750','1000','1500', '2000'},'fontsize',7, ...
    'ylim',[0 0.17],'ytick', [0:0.05:0.15], 'yticklabel', {'0', '5', '10', '15'}, 'ticklength', [0.02 0.025]);

plotshaded([0 0.6], [0 0; max(get(ha31, 'ylim')) max(get(ha31, 'ylim'))], [0.2 0.8 0.2])

% % Before lesion
t_bins       =  xedges.RelT;
% Lesion group
% Before operative
mean_prob  =  TBTg.mean_RelTdist(TBTg.Group==grpName(1)+"-"+"Pre",:);
sem_prob    =  TBTg.sem_RelTdist(TBTg.Group==grpName(1)+"-"+"Pre",:);
plotshaded(t_bins, [mean_prob-sem_prob; mean_prob+sem_prob], cGray);
% After operative
mean_prob  =  TBTg.mean_RelTdist(TBTg.Group==grpName(1)+"-"+"Post",:);
sem_prob    =  TBTg.sem_RelTdist(TBTg.Group==grpName(1)+"-"+"Post",:);
plotshaded(t_bins, [mean_prob-sem_prob; mean_prob+sem_prob], cOrange);
 
plot(xedges.RelT, TBTg.mean_RelTdist(TBTg.Group==grpName(1)+"-"+"Pre",:), 'linewidth', 1.25, 'color', 'k')
plot(xedges.RelT, TBTg.mean_RelTdist(TBTg.Group==grpName(1)+"-"+"Post",:), 'linewidth', 1.25, 'color', cOrange)

xlabel('Release time (ms)')
ylabel('Probability') 

% Plot CDF 
ha32 = axes;
set(ha32, 'units', 'centimeters', 'position', [xs(1)+size4(1)+1.5 ys(2) size4], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0 1],'xtick',[0:0.25:1],'xticklabel',{'0','250', '500','750','1000','1500', '2000'},'fontsize',7, ...
    'ylim',[0 1], 'ytick', [0:0.25:1], 'yticklabel', {'0', '25', '50', '75', '100'},'ticklength', [0.02 0.025]);

plotshaded([0 0.6], [0 0; max(get(ha32, 'ylim')) max(get(ha32, 'ylim'))], [0.2 0.8 0.2])
% Lesion group
% Before operative
t_bins       =  xedges.RelT;

mean_cdf  =  TBTg.mean_RelTcdf(TBTg.Group==grpName(1)+"-"+"Pre",:);
sem_cdf    =  TBTg.sem_RelTcdf(TBTg.Group==grpName(1)+"-"+"Pre",:);

RelTMedBeforeL_All = TBTg.mean_RelTMed_All(TBTg.Group==grpName(1)+"-"+"Pre");
sprintf('Release time before lesion (Lesion group, all FPs): %2.2f', 1000*RelTMedBeforeL_All)

plotshaded(t_bins, [mean_cdf-sem_cdf; mean_cdf+sem_cdf], cGray)
% After operative
t_bins       =  xedges.RelT;
mean_cdf  =  TBTg.mean_RelTcdf(TBTg.Group==grpName(1)+"-"+"Post",:);
sem_cdf    =  TBTg.sem_RelTcdf(TBTg.Group==grpName(1)+"-"+"Post",:);

RelTMedAfterL_All = TBTg.mean_RelTMed_All(TBTg.Group==grpName(1)+"-"+"Post");
sprintf('Release time after lesion (Lesion group, all FPs): %2.2f', 1000*RelTMedAfterL_All)
plotshaded(t_bins, [mean_cdf-sem_cdf; mean_cdf+sem_cdf], cOrange)

plot(xedges.RelT, TBTg.mean_RelTcdf(TBTg.Group==grpName(1)+"-"+"Pre",:), 'linewidth', 1, 'color', 'k')
plot(xedges.RelT, TBTg.mean_RelTcdf(TBTg.Group==grpName(1)+"-"+"Post",:), 'linewidth', 1, 'color', cOrange)

xlabel('Release time (ms)')
ylabel('CDF')

% short
RelTMedPreLMean(1) = TBTg.mean_RelTMed_S(TBTg.Group==grpName(1)+"-"+"Pre");
RelTMedPreLSEM(1) = TBTg.sem_RelTMed_S(TBTg.Group==grpName(1)+"-"+"Pre");
RelTMedPostLMean(1) = TBTg.mean_RelTMed_S(TBTg.Group==grpName(1)+"-"+"Post");
RelTMedPostLSEM(1) = TBTg.sem_RelTMed_S(TBTg.Group==grpName(1)+"-"+"Post");
sprintf('Release time before lesion: %2.0f, after lesion: %2.0f (in ms, Lesion group, short FP)', 1000*RelTMedPreLMean(1),  1000*RelTMedPostLMean(1))

% medium
RelTMedPreLMean(2) = TBTg.mean_RelTMed_M(TBTg.Group==grpName(1)+"-"+"Pre");
RelTMedPreLSEM(2) = TBTg.sem_RelTMed_M(TBTg.Group==grpName(1)+"-"+"Pre");
RelTMedPostLMean(2) = TBTg.mean_RelTMed_M(TBTg.Group==grpName(1)+"-"+"Post");
RelTMedPostLSEM(2) = TBTg.sem_RelTMed_M(TBTg.Group==grpName(1)+"-"+"Post");
sprintf('Release time before lesion: %2.0f, after lesion: %2.0f (Lesion group, medium FP)', 1000*RelTMedPreLMean(2),  1000*RelTMedPostLMean(2))

% long
RelTMedPreLMean(3) = TBTg.mean_RelTMed_L(TBTg.Group==grpName(1)+"-"+"Pre");
RelTMedPreLSEM(3) = TBTg.sem_RelTMed_L(TBTg.Group==grpName(1)+"-"+"Pre");
RelTMedPostLMean(3) = TBTg.mean_RelTMed_L(TBTg.Group==grpName(1)+"-"+"Post");
RelTMedPostLSEM(3) = TBTg.sem_RelTMed_L(TBTg.Group==grpName(1)+"-"+"Post");
sprintf('Release time before lesion: %2.0f, after lesion: %2.0f (Lesion group, long FP)', 1000*RelTMedPreLMean(3),  1000*RelTMedPostLMean(3))

% add median release time
ha321 = axes;
set(ha321, 'units', 'centimeters', 'position', [xs(1)+(size4(1)+1.5)*2 ys(2) size4], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[250 1750],'xtick',[500 1000 1500],'xticklabel',{'500', '1000', '1500'},'fontsize',7, ...
    'ylim',[200 650],'ticklength', [0.02 0.025]);
 
plot(fplist*1000, 1000*RelTMedPreLMean, 'o-', 'linewidth', 1, 'color', 'k', 'markerfacecolor', 'k', 'markeredgecolor','w', 'markersize', 6)
line([fplist*1000; fplist*1000], 1000*[RelTMedPreLMean-RelTMedPreLSEM; RelTMedPreLMean+RelTMedPreLSEM], ...
    'color','k', 'linewidth', 1)

plot(fplist*1000, 1000*RelTMedPostLMean, 'o-', 'linewidth', 1,  'color', cOrange,'markerfacecolor', cOrange, 'markeredgecolor','w', 'markersize', 6)
line([fplist*1000; fplist*1000], 1000*[RelTMedPostLMean-RelTMedPostLSEM; RelTMedPostLMean+RelTMedPostLSEM], ...
    'color',cOrange, 'linewidth', 1)

xlabel('Foreperiod (ms)')
ylabel('Release time (ms)')
 
% add information
ha321info = axes;
set(ha321info, 'units', 'centimeters', 'position', [xs(1)+(size4(1)+1.5)*2+size4(1) ys(2) 2 size4(2)], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0 10],'xtick',[ ],'xticklabel',{},'fontsize',7, ...
    'ylim',[0 10],'ticklength', [0.02 0.025]);

plot([1 4], [4 4], '-', 'linewidth', 1, 'color', 'k', 'markerfacecolor', 'k', 'markeredgecolor','w', 'markersize', 6)
plot([2.5], [4], 'o', 'linewidth', 1, 'color', 'k', 'markerfacecolor', 'k', 'markeredgecolor','w', 'markersize', 6)

text(5, 4, 'Pre')
plot([1 4], [2 2], '-', 'linewidth', 1, 'color',cOrange, 'markerfacecolor', cOrange, 'markeredgecolor','w', 'markersize', 6)
plot([2.5], [2], 'o', 'linewidth', 1, 'color',cOrange, 'markerfacecolor', cOrange, 'markeredgecolor','w', 'markersize', 6)

text(5, 2, 'Post')
axis off
text(1, 7, 'Lesion', 'fontsize', 8, 'fontweight', 'bold')

%% PLOT Sham group, RT distribution, legend: S/M/L * Pre/Post
% Release time
ha41 = axes;
set(ha41, 'units', 'centimeters', 'position', [xs(1) ys(3) size4], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0 1],'xtick',[0:0.25:1],'xticklabel',{'0','250', '500','750','1000','1500', '2000'},'fontsize',7, ...
    'ylim',[0 0.17],'ytick', [0:0.05:0.15], 'yticklabel', {'0', '5', '10', '15'}, 'ticklength', [0.02 0.025]);

plotshaded([0 0.6], [0 0; max(get(ha41, 'ylim')) max(get(ha41, 'ylim'))], [0.2 0.8 0.2])

% % Before lesion
t_bins       =  xedges.RelT;
% Lesion group
% Before operative
mean_prob  =  TBTg.mean_RelTdist(TBTg.Group==grpName(2)+"-"+"Pre",:);
sem_prob    =  TBTg.sem_RelTdist(TBTg.Group==grpName(2)+"-"+"Pre",:);
plotshaded(t_bins, [mean_prob-sem_prob; mean_prob+sem_prob], cGray);
% After operative
mean_prob  =  TBTg.mean_RelTdist(TBTg.Group==grpName(2)+"-"+"Post",:);
sem_prob    =  TBTg.sem_RelTdist(TBTg.Group==grpName(2)+"-"+"Post",:);
plotshaded(t_bins, [mean_prob-sem_prob; mean_prob+sem_prob], cOrange);
 
plot(xedges.RelT, TBTg.mean_RelTdist(TBTg.Group==grpName(2)+"-"+"Pre",:), 'linewidth', 1.25, 'color', 'k')
plot(xedges.RelT, TBTg.mean_RelTdist(TBTg.Group==grpName(2)+"-"+"Post",:), 'linewidth', 1.25, 'color', cOrange)

xlabel('Release time (ms)')
ylabel('Probability') 

% Plot CDF 
ha42 = axes;
set(ha42, 'units', 'centimeters', 'position', [xs(1)+size4(1)+1.5 ys(3) size4], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0 1],'xtick',[0:0.25:1],'xticklabel',{'0','250', '500','750','1000','1500', '2000'},'fontsize',7, ...
    'ylim',[0 1], 'ytick', [0:0.25:1], 'yticklabel', {'0', '25', '50', '75', '100'},'ticklength', [0.02 0.025]);

plotshaded([0 0.6], [0 0; max(get(ha32, 'ylim')) max(get(ha32, 'ylim'))], [0.2 0.8 0.2])
% Lesion group
% Before operative
t_bins       =  xedges.RelT;

mean_cdf  =  TBTg.mean_RelTcdf(TBTg.Group==grpName(1)+"-"+"Pre",:);
sem_cdf    =  TBTg.sem_RelTcdf(TBTg.Group==grpName(1)+"-"+"Pre",:);

RelTMedBeforeL_All = TBTg.mean_RelTMed_All(TBTg.Group==grpName(2)+"-"+"Pre");
sprintf('Release time before lesion (Lesion group, all FPs): %2.2f', 1000*RelTMedBeforeL_All)

plotshaded(t_bins, [mean_cdf-sem_cdf; mean_cdf+sem_cdf], cGray)
% After operative
t_bins       =  xedges.RelT;
mean_cdf  =  TBTg.mean_RelTcdf(TBTg.Group==grpName(2)+"-"+"Post",:);
sem_cdf    =  TBTg.sem_RelTcdf(TBTg.Group==grpName(2)+"-"+"Post",:);

RelTMedAfterL_All = TBTg.mean_RelTMed_All(TBTg.Group==grpName(2)+"-"+"Post");
sprintf('Release time after lesion (Lesion group, all FPs): %2.2f', 1000*RelTMedAfterL_All)
plotshaded(t_bins, [mean_cdf-sem_cdf; mean_cdf+sem_cdf], cOrange)

plot(xedges.RelT, TBTg.mean_RelTcdf(TBTg.Group==grpName(2)+"-"+"Pre",:), 'linewidth', 1, 'color', 'k')
plot(xedges.RelT, TBTg.mean_RelTcdf(TBTg.Group==grpName(2)+"-"+"Post",:), 'linewidth', 1, 'color', cOrange)

xlabel('Release time (ms)')
ylabel('CDF')

% short
RelTMedPreLMean(1) = TBTg.mean_RelTMed_S(TBTg.Group==grpName(2)+"-"+"Pre");
RelTMedPreLSEM(1) = TBTg.sem_RelTMed_S(TBTg.Group==grpName(2)+"-"+"Pre");
RelTMedPostLMean(1) = TBTg.mean_RelTMed_S(TBTg.Group==grpName(2)+"-"+"Post");
RelTMedPostLSEM(1) = TBTg.sem_RelTMed_S(TBTg.Group==grpName(2)+"-"+"Post");
sprintf('Release time before lesion: %2.0f, after lesion: %2.0f (in ms, Lesion group, short FP)', ...
    1000*RelTMedPreLMean(1),  1000*RelTMedPostLMean(1))

% medium
RelTMedPreLMean(2) = TBTg.mean_RelTMed_M(TBTg.Group==grpName(2)+"-"+"Pre");
RelTMedPreLSEM(2) = TBTg.sem_RelTMed_M(TBTg.Group==grpName(2)+"-"+"Pre");
RelTMedPostLMean(2) = TBTg.mean_RelTMed_M(TBTg.Group==grpName(2)+"-"+"Post");
RelTMedPostLSEM(2) = TBTg.sem_RelTMed_M(TBTg.Group==grpName(2)+"-"+"Post");
sprintf('Release time before lesion: %2.0f, after lesion: %2.0f (Sham group, medium FP)', ...
    1000*RelTMedPreLMean(2),  1000*RelTMedPostLMean(2))

% long
RelTMedPreLMean(3) = TBTg.mean_RelTMed_L(TBTg.Group==grpName(2)+"-"+"Pre");
RelTMedPreLSEM(3) = TBTg.sem_RelTMed_L(TBTg.Group==grpName(2)+"-"+"Pre");
RelTMedPostLMean(3) = TBTg.mean_RelTMed_L(TBTg.Group==grpName(2)+"-"+"Post");
RelTMedPostLSEM(3) = TBTg.sem_RelTMed_L(TBTg.Group==grpName(2)+"-"+"Post");
sprintf('Release time before lesion: %2.0f, after lesion: %2.0f (Sham group, long FP)', 1000*RelTMedPreLMean(3),  1000*RelTMedPostLMean(3))

% add median release time
ha421 = axes;
set(ha421, 'units', 'centimeters', 'position', [xs(1)+(size4(1)+1.5)*2 ys(3) size4], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[250 1750],'xtick',[500 1000 1500],'xticklabel',{'500', '1000', '1500'},'fontsize',7, ...
    'ylim',[200 650],'ticklength', [0.02 0.025]);
 
plot(fplist*1000, 1000*RelTMedPreLMean, 'o-', 'linewidth', 1, 'color', 'k', 'markerfacecolor', 'k', 'markeredgecolor','w', 'markersize', 6)
line([fplist*1000; fplist*1000], 1000*[RelTMedPreLMean-RelTMedPreLSEM; RelTMedPreLMean+RelTMedPreLSEM], ...
    'color','k', 'linewidth', 1)

plot(fplist*1000, 1000*RelTMedPostLMean, 'o-', 'linewidth', 1,  'color', cOrange,'markerfacecolor', cOrange, 'markeredgecolor','w', 'markersize', 6)
line([fplist*1000; fplist*1000], 1000*[RelTMedPostLMean-RelTMedPostLSEM; RelTMedPostLMean+RelTMedPostLSEM], ...
    'color',cOrange, 'linewidth', 1)

xlabel('Foreperiod (ms)')
ylabel('Release time (ms)')

% add information
ha421info = axes;
set(ha421info, 'units', 'centimeters', 'position', [xs(1)+(size4(1)+1.5)*2+size4(1) ys(3) 2 size4(2)], 'nextplot', 'add','tickDir', 'out',...
    'xlim',[0 10],'xtick',[ ],'xticklabel',{},'fontsize',7, ...
    'ylim',[0 10],'ticklength', [0.02 0.025]);
  
axis off
text(1, 7, 'Sham', 'fontsize', 8, 'fontweight', 'bold')

%% Save
tosavename = fullfile(pwd,'LesionAfterLearning_J');
% saveas(hf,savename,'fig');
print(hf,'-dpng',tosavename);
print(hf,'-dpdf',tosavename,'-bestfit');
% print(hf,'-depsc2',savename);

exportgraphics(hf, [tosavename '.eps'],'ContentType','vector')
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
t_RelTbins = movmean(edges_RelT,2,'Endpoints','discard');

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
            edges_RelT,'Normalization','probability'),2,'gaussian',smo_win);
        t.RelTdist_S = smoothdata(histcounts(tdata.HT((idxCor|idxLate)&idxFPS)-tdata.FP((idxCor|idxLate)&idxFPS),...
            edges_RelT,'Normalization','probability'),2,'gaussian',smo_win);
        t.RelTdist_M = smoothdata(histcounts(tdata.HT((idxCor|idxLate)&idxFPM)-tdata.FP((idxCor|idxLate)&idxFPM),...
            edges_RelT,'Normalization','probability'),2,'gaussian',smo_win);
        t.RelTdist_L = smoothdata(histcounts(tdata.HT((idxCor|idxLate)&idxFPL)-tdata.FP((idxCor|idxLate)&idxFPL),...
            edges_RelT,'Normalization','probability'),2,'gaussian',smo_win);

        % added by JY
        t.RelTcdf = smoothdata(histcounts(tdata.HT(idxCor|idxLate)-tdata.FP(idxCor|idxLate),...
            edges_RelT,'Normalization','cdf'),2,'gaussian',smo_win);
        t.NumRelT = length(find(idxCor|idxLate)); % index possibly useful for tracking sample size
        t.RelTMed_All = med_cdf(t_RelTbins, t.RelTcdf);

        t.RelTcdf_S = smoothdata(histcounts(tdata.HT((idxCor|idxLate)&idxFPS)-tdata.FP((idxCor|idxLate)&idxFPS),...
            edges_RelT,'Normalization','cdf'),2,'gaussian',smo_win);        
        t.NumRelT_S = length(find((idxCor|idxLate)&idxFPS)); % index possibly useful for tracking sample size
        t.RelTMed_S =  med_cdf(t_RelTbins, t.RelTcdf_S)

        t.RelTcdf_M = smoothdata(histcounts(tdata.HT((idxCor|idxLate)&idxFPM)-tdata.FP((idxCor|idxLate)&idxFPM),...
            edges_RelT,'Normalization','cdf'),2,'gaussian',smo_win);
        t.NumRelT_M = length(find((idxCor|idxLate)&idxFPM)); % index possibly useful for tracking sample size
        t.RelTMed_M =  med_cdf(t_RelTbins, t.RelTcdf_M)

        t.RelTcdf_L = smoothdata(histcounts(tdata.HT((idxCor|idxLate)&idxFPL)-tdata.FP((idxCor|idxLate)&idxFPL),...
            edges_RelT,'Normalization','cdf'),2,'gaussian',smo_win);
        t.NumRelT_L = length(find((idxCor|idxLate)&idxFPL)); % index possibly useful for tracking sample size
        t.RelTMed_L =  med_cdf(t_RelTbins, t.RelTcdf_L)
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

function tmed = med_cdf(t_bins, meancdf)

t_bins2 = [t_bins(1):1/1000:t_bins(end)];
meancdf2 = interp1(t_bins, meancdf, t_bins2);
[~, indmed] = min(abs(meancdf2-0.5));
tmed = t_bins2(indmed);

end