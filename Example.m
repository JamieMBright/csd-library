% Example of all CSD methodologies
% ------------------------------------------------------------------------
%                           INPUT PARAMETERS
% ------------------------------------------------------------------------
% ASSUMPTION: all data is in 1-min resolution and all vectors perfectly 
% match each other in terms of time stamps.
%
% ghi = Global horizontal irradiance, column vector. Not necessary
%       to be continuous so long as the ghics perfectly corresponds.
% 
% ghics = Clear-sky global horizontal irradiance, column vector. 
%
% dni = Direct normal irradiance, column vector. Not necessary
%       to be continuous so long as the dnics perfectly corresponds.
% 
% dnics = Clear-sky direct normal irradiance, column vector. 
%
% dif = Diffuse horizontal irradiance, column vector. Not necessary
%       to be continuous so long as the difcs perfectly corresponds.
% 
% difcs = Clear-sky diffuse horizontal irradiance, column vector. 
%
% zen = Zenith angle in degrees. Corresponding to all inputs
%
% exth = Horizontal projection of extraterrestrial irradiance. 
%
% aod = Total aerosol extinction at 550nm.
%
% plot_figure = if this variable is defined by any class, then a figure is
%               plotted illustrating the outcome of the CSD method, e.g. 1;
% ------------------------------------------------------------------------
%                              OUTPUTS
% ------------------------------------------------------------------------
% csd = clear-sky detection. a flag of 1 means that clouds were detected
%       whereas a flag of 0 means that the hour is clear.
%% Pre-amble
clearvars
close all
addpath('models\')
%% Data
% all sample data is taken from BSRN 2013 at Adelaide Airport.
data = load('sample_data.mat');data=data.data;
% sample ghi data from Adelaide Airport, sz = [3000,1]
ghi = data(:,1);
% sample corresponding REST2 clear-sky ghi 
ghics = data(:,2);
% sample dni data from Adelaide Airport, 
dni = data(:,3);
% sample corresponding REST2 clear-sky dni
dnics = data(:,4);
% sample dif data from Adelaide Airport
dif = data(:,5);
% sample corresponding REST2 clear-sky dif
difcs = data(:,6);
% sample Eexth calculated at Adelaide Airport
exth = data(:,7);
% sample Zenith angle for Adelaide Airport
zen = data(:,8);
% sample MERRA2 aerosol extinction at 550nm
aod = data(:,9);
% sample LST for Adelaide Airport
LST = datevec(data(:,10));

%% CSD methods
% produce clear-sky detection from each model available.
perez        = Perez1990CSD(dif,dni,zen);
batlles      = Batlles1999CSD(ghi,dif,exth,zen);
longackerman = LongAckerman2000CSD(ghi,dif,zen,LST);
ineichen06   = Ineichen2006CSD(exth,dni,zen);
ineichen09   = Ineichen2009CSD(ghi,exth,zen);
polo         = Polo2009CSD(ghi,ghics,zen,LST);
xie          = Xie2013CSD(ghi,ghics,dni,dnics,zen,exth);
gueymard     = Gueymard2013CSD(dni,dnics,zen);
lefevre      = Lefevre2013CSD(ghi,dif,exth,zen);
garcia       = Garcia2014CSD(ghi,dif,zen,aod,LST);
quesadaruiz  = Quesadaruiz2015CSD(ghi,ghics);
reno         = Reno2016CSD(ghi,ghics);
inman        = Inman2015CSD(ghi,ghics,dni,dnics);
ineichen16   = Ineichen2016CSD(ghi,dni,dif,zen,exth,aod);
aliamartinez = AliaMartinez2016CSD(ghi,ghics,zen,LST);
larraneta    = Larraneta2017CSD(dni,dnics,zen);
shen         = Shen2018CSD(dni,dnics,aod);
ellis        = Ellis2018CSD(ghi,ghics);
zhang        = Zhang2018CSD(ghi,ghics);
ruizarias    = RuizArias2018CSD(ghi,dif,zen,LST);

%% Figures
% produce a large figure showing the data from all models
methods = {'Perez','Batlles','LongAckerman','Ineichen06','Ineichen09','Polo','Garcia','Xie','Gueymard','Lefevre','QuesadaRuiz','Reno','Inman','Ineichen16','AliaMartinez','Larraneta','Shen','Ellis','Zhang','RuizArias'};
ymax = 1300;
f=figure('name','CSD methods','color','w');
f.Position = [250,250,1200,500];
for i =1:length(methods)
    subplot(ceil(length(methods)/2),2,i)
    hold on
    plot(ghi)
    plot(ghics,'k:')
    CSD = ghi;
    csd = eval(lower(methods{i}));
    CSD(csd==1)=NaN;
    plot(CSD,'linewidth',2,'color','k')
    hold off
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    ylim([0 ymax])
    ylabel('GHI','FontSize',5)
    xlim([1,length(ghi)])
    text(length(ghi)/2,ymax.*1.2,methods{i},'HorizontalAlignment','Center')
    if i>length(methods)-2
        xlabel('Time','FontSize',8)
    end
end
