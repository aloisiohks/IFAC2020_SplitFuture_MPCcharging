function figFormat(width,height,xlabels,subplab,plotyyFlag)
% FIGFORMAT Set graphics to (width x height) axes frame; set line thicknesses
%
% Call this routine after plotting, in order to make it suitable for
% presentation
%
%	width = width of axis frame, default: golden section to height
%	height = height of axis frame, default: golden sect of width, or 3.5 inches
%	xlabels = control over space between frame and x-axis
%	    'xtick': allow xticklabels
%	    'nolabels': no labels at all
%	    'plot' (default): allow xticklabels and xlabel (for x, etc)
%	    'subplot': allow xticklabels and xlabel (for x and (a) etc)
% subplab = '(a)' etc, only for subplot
%	The default frame size is 3.5 x 5.66 inches

% GLP Global style variables.
global THINLINEWIDTH MEDLINEWIDTH THICKLINEWIDTH GRAY
global SMALLMARKERSIZE LARGEMARKERSIZE
THINLINEWIDTH = 0.8;
MEDLINEWIDTH = 1.0;
THICKLINEWIDTH = 2.5;
SMALLMARKERSIZE = 10;
LARGEMARKERSIZE = 14;
GRAY = 'g';

% First, account for the x-axis labels
if nargin<3, 
  xlabels=''; 
end
if isempty(xlabels), 
  xlabels='plot'; 
end
if strcmp(xlabels,'nolabels'), 
  fcdist=0.1;
elseif strcmp(xlabels,'xtick'), 
  fcdist=0.1;
elseif strcmp(xlabels,'plot'), 
  fcdist=0.6;
elseif strcmp(xlabels,'subplot'), 
  fcdist=0.9;
else
  error(['xlabels = ''',xlabels,''' is not valid']);
end

% Now, account for the subplot title
if nargin<4, 
  subplab=''; 
end
if nargin<5,
  plotyyFlag = 0;
end

% Take care of the width and height attributes
if nargin<1, 
  width = NaN; 
end
if isempty(width),
  width = NaN;
end
if nargin<2, 
  height = NaN; 
end
if isempty(height),
  height = NaN;
end

% Compute the golden ratio and see whether user-specified dimensions agree
% with this "perfect" shape
gold=(sqrt(5)-1)/2; %0.6180
if ~isnan(height)&~isnan(width)
  if (height~=width)&(abs(height/width-gold)>1e-2)
    disp('Warning! Ratio differs from golden section in figSetup')
  end
end
if isnan(height)&isnan(width), 
  height=3.5; 
end
if isnan(height), 
  height=gold*width; 
end
if isnan(width), 
  width=height/gold; 
end

TopMarg = 0.5;	BotMarg = fcdist;
LefMarg = 1;	
if plotyyFlag,
  RiMarg = 0.75;
else
  RiMarg = 0.25;
end

pwidth = width+RiMarg+LefMarg;
pheight = height+TopMarg+BotMarg;

% Set default drawing properties
set(gcf,'DefaultAxesLineWidth',THINLINEWIDTH);
set(gcf,'DefaultLineLineWidth',THICKLINEWIDTH);
set(gca,'LineWidth',THINLINEWIDTH);
% for presentations, axis font size 14; legend 16; x/ylabel 18; title 20
set(gca,'FontSize',14);
set(gca,'XMinorTick','On');
set(gca,'YMinorTick','On');
set(get(gca,'XLabel'),'FontSize',18);
set(get(gca,'YLabel'),'FontSize',18);
set(get(gca,'Title'),'FontSize',20);

%  set axes position and frame size
funits=get(gcf,'units');
aunits=get(gca,'units');
set( gca, 'units', 'inches')
set( gca, 'position', [ LefMarg fcdist width height ])
set( gcf, 'units','inches');
a = get(gcf,'position');
set(gcf,'position',[a(1) a(2) pwidth pheight])
set(gcf,'units',funits);
set( gca, 'units',aunits)
set(findobj(gca,'Type','Line'),'LineWidth',THICKLINEWIDTH);

figure(gcf)

if strcmp(xlabels,'subplot')
  if ~isempty(subplab)
    ax=axis;
    if strcmp(get(gca,'XScale'),'linear'), xc=mean(ax(1:2));
    elseif strcmp(get(gca,'XScale'),'log'), xc=sqrt(prod(ax(1:2)));
    end
    if strcmp(get(gca,'YScale'),'linear'), yc=ax(3)-(fcdist-0.1)*(ax(4)-ax(3))/height;
    elseif strcmp(get(gca,'YScale'),'log'), yc=ax(3)/(ax(2)/ax(1))^0.12; % broken?
    end
    delete(findobj(gca,'Type','text','string',subplab))
    text(xc,yc,subplab,'fontsize',14,...
    	'Horizontalalignment','center','Verticalalignment','bottom');
  end
end

drawnow; pause(0.01);