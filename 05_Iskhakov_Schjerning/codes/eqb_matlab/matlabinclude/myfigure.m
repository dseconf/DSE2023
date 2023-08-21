function h=myfigure(arg)
  % This function creates the figure with white background
  % that can be fullscreen, and creates axes inside of it
  % arg is a string that can contain:
  %     fullscreen, and/or
  %     figure, and/or
  %     name for the figure
  error('Should not be used: migrate to graphs.myfigure'); 
  if nargin>0
    args=strsplit(arg);
    name=setdiff(args,{'fullscreen','figure'}); %remove preset args
    if numel(name)==1
      name=name{1}; %convert to string
    elseif numel(name)>0
      name=join(name);
      name=name{1};
    else
      name='';
    end
  else
    args={};
    name='';
  end
  if ismember('fullscreen',args)
    %default position is full screen
    position=get(0,'ScreenSize');
    position([1 2])=0;
    h=figure('Color',[1 1 1],'Position',position,'NextPlot','new','Name',name);
  else
    h=figure('Color',[1 1 1],'NextPlot','new','Name',name);
  end
  if ~ismember('figure',args)
    h=axes('Parent',h);
    hold(h,'all');
  end
end