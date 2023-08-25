% Matlab class implementing a set of tools to work with linearly interpolated functions
% Written by Fedor Iskhakov, Australian National University, 2016

classdef polyline
% This class defines a function linearly interpolated on a grid
% and implements the tools for working with arrays of these lines
%
% Objects of this class can be put in arrays and manipulated in bulk
% Methods: len() to return the number of points
%					 sort() to sort the points
%					 interpolate() to compute the function 
%					 upper_envelope() to find upper envelope of a set of lines
%					 inset() to insert a new point
%					 chop() to break one line into two
% 				 grow() to join multiple lines together
%					 thinout() to delete points by index
%					 diff() to compare two polylines
%					 plot() to make a plot of the function
%
% Fedor Iskhakov, 2016
% Australian National University
%
properties (Access=public)
	x=[]; 		%x-grid points
	y=[]; 		%functional values
	label='';	%label of this line
	% correct_for_incomplete_loopbacks_for_dcegm=true;
end
methods
	%----------------------------------------------
	function obj = polyline(varargin)
		%create interpolated function
		%optional inputs: x, y, label
		if nargin>0
			obj.x=reshape(varargin{1},1,[]);
		end
		if nargin>1
			obj.y=reshape(varargin{2},1,[]);
		end
		if nargin>2 && ischar(varargin{3})
			obj.label=varargin{3};
		end
	end
	%----------------------------------------------
	% function res=get.len(obj)
	function res=len(obj)
		%number of gridpoints in this line/function
		%if obj is not scalar, multiple outputs are produced
		for k=1:numel(obj)
			res(k)=numel(obj(k).x);
			if res(k)~=numel(obj(k).y);
				error (sprintf('Polyline with mismatch of the number of points: x has %d while y has %d elements',numel(obj(k).x),numel(obj(k).y)))
			end
		end
		res=reshape(res,size(obj)); %return the same shape when numel(obj)>1
	end
	%----------------------------------------------
	function res=sort(obj)
		%this function sorts the polylines on x grid
		%returns the sorted polylines, original obj is not modified
		for i=1:numel(obj)
			[~,i1]=sort(obj(i).x);
			res(i)=polyline(obj(i).x(i1),obj(i).y(i1),obj(i).label);
		end
		res=reshape(res,size(obj));
	end
	%----------------------------------------------
	function [res extrapflag]=interpolate(obj,x)
		%compute the func value in given point(s)
		%linear extrapolation allowed, indicator is created if nargsout==2
		if numel(obj)==1
			%when polyline is scalar (only one)
			try
				res=interp1(obj.x,obj.y,x,'linear','extrap');
			catch ME
				steps=reshape(obj.x(2:end)-obj.x(1:end-1),[],1);
				fprintf('\n\nInterpolation error!\npolyline dims are (%d,%d)\nRepeated points: %d\nInverted steps: %d\n', ...
					size(obj,1),size(obj,2),sum(abs(steps)<eps),sum(steps<0));
				fprintf('Identifier: %s\nCall stack: ',ME.identifier);
				for i=0:numel(ME.stack)-1
					if i>0
						fprintf(' >>> ');
					end
					fprintf('%s (line %d)',ME.stack(end-i).name,ME.stack(end-i).line)
				end
				fprintf('\n')
				fprintf('INTERACTIVE KEYBOARD MODE\nVariable steps contains the steps of the grid\n');
				keyboard
			end
			if nargout>1
				extrapflag=~(x>=min(obj.x) & x<=max(obj.x));
			end
		else
			%interpolate each polyline if there are many
			%output matrix with rows for each polyline, colums of interpolated values
			xx=reshape(x,1,[]);%but make x row vector to output one matrix
			for k=1:numel(obj)
				try
					res(k,1:numel(x))=interp1(obj(k).x,obj(k).y,xx,'linear','extrap');
				catch ME
					steps=@(nn) reshape(obj(nn).x(2:end)-obj(nn).x(1:end-1),[],1);
					fprintf('\n\nInterpolation error!\nk=%d, polyline dims are (%d,%d)\nRepeated points: %d\nInverted steps: %d\n', ...
						k,size(obj,1),size(obj,2),sum(abs(steps(k))<eps),sum(steps(k)<0));
					fprintf('Identifier: %s\nCall stack: ',ME.identifier);
					for i=0:numel(ME.stack)-1
						if i>0
							fprintf(' >>> ');
						end
						fprintf('%s (line %d)',ME.stack(end-i).name,ME.stack(end-i).line)
					end
					fprintf('\n')
					fprintf('INTERACTIVE KEYBOARD MODE\nUse steps(i) function to analyze the grids\n');
					keyboard
				end
				if nargout>1
					extrapflag(k,1:numel(x))=~(xx>=min(obj(k).x) & xx<=max(obj(k).x));
				end
			end
		end
	end
	%----------------------------------------------
	function res=grow(obj,in2,front)
		%add polyline2 to the end (or the front) of polyline1
		%many-to-many extending is allowed
		%return is matrix with obj in rows, in2 in columns
		if nargin<3
			front=false;%to the end by default
		end
		for k1=1:numel(obj)
			for k2=1:numel(in2)
				if front
					res(k1,k2)=polyline( ...
						[reshape(in2(k2).x,1,[]) reshape(obj(k1).x,1,[])], ...
						[reshape(in2(k2).y,1,[]) reshape(obj(k1).y,1,[])], ...
						sprintf('%s + %s',in2(k2).label,obj(k1).label));
				else
					res(k1,k2)=polyline( ...
						[reshape(obj(k1).x,1,[]) reshape(in2(k2).x,1,[])], ...
						[reshape(obj(k1).y,1,[]) reshape(in2(k2).y,1,[])], ...
						sprintf('%s + %s',obj(k1).label,in2(k2).label));
				end
			end
		end
	end
	%----------------------------------------------
	function res=inset(obj,x,y,j)
		%this function inserts a point (x,y) after position j
		%in the current grid (so j=0 grows first point)
		%if x and y are vectors, multiple points are inserted
		%if j is missing, the point is added in the end
		%inset into each polyline if there are many
		%Cases in the function to minimize overal runtime
		res=obj;
		if numel(x)>1
			mx=numel(x);
			my=numel(y);
			for k=1:numel(obj)
				if nargin<=3 || j>=numel(obj(k).x);
					res(k).x(end+1:end+mx)=x;
					res(k).y(end+1:end+my)=y;
				else
					res(k).x(j+mx+1:end+mx)=obj(k).x(j+1:end);
					res(k).x(j+1:j+mx)=x(:);
					res(k).y(j+my+1:end+my)=obj(k).y(j+1:end);
					res(k).y(j+1:j+my)=y(:);
				end
			end
		elseif numel(obj)>1
			for k=1:numel(obj)
				if nargin<=3 || j>=numel(obj(k).x);
					res(k).x(end+1)=x;
					res(k).y(end+1)=y;
				else
					res(k).x(j+2:end+1)=obj(k).x(j+1:end);
					res(k).x(j+1)=x;
					res(k).y(j+2:end+1)=obj(k).y(j+1:end);
					res(k).y(j+1)=y;
				end
			end
		else
			if nargin<=3 || j>=numel(obj.x);
				res.x(end+1)=x;
				res.y(end+1)=y;
			else
				res.x(j+2:end+1)=obj.x(j+1:end);
				res.x(j+1)=x;
				res.y(j+2:end+1)=obj.y(j+1:end);
				res.y(j+1)=y;
			end
		end
	end
	%----------------------------------------------
	function res=thinout(obj,indx)
		%this removes the indexed points from polylines
		res=obj;%copy input to output
		for k=1:numel(res)
			ii=intersect(1:numel(res(k).x),indx);
			res(k).x(ii)=[];
			res(k).y(ii)=[];
		end
	end
	%----------------------------------------------
	function indx=diff(obj,pl2,significance)
		%this function returns the indexes of points in obj that are not in pl2
		%obj can have multiple elements, pl2 is treated as scalar polyline
		if nargin<3
			significance=5; %equality is measured up to 10^-singif
		end
		x1=round(pl2(1).x*(10^significance)) * 10^(-significance);
		y1=round(pl2(1).y*(10^significance)) * 10^(-significance);
		for k=1:numel(obj)
			x=round(obj(k).x*(10^significance)) * 10^(-significance);
			y=round(obj(k).y*(10^significance)) * 10^(-significance);
			indx{k}=find(~ismember(x,x1) | ~ismember(y,y1));
		end
		if numel(indx)==1
			indx=indx{1}; %return vector if scalar obj
		else
			indx=reshape(indx,size(obj));
		end
	end
	%----------------------------------------------
	function [res1 res2]=chop(obj,j,repeat)
		%separate the grid into 1,..,j and j+1,.. parts
		%if repeat=true, the boundary points are repeated in both resulting polylines
		if nargin<2
			error 'Have to have one input for .chop'
		end
		%chop each polyline if there are many
		for k=1:numel(obj)
			if j>obj(k).len
				warning 'Producing empty polyline by chopping at index j>len'
				j=obj(k).len;
			end
			res1(k)=polyline(obj(k).x(1:j),obj(k).y(1:j), ...
										sprintf('%s (1:%d)',obj(k).label,j));
			if nargout>1
				if exist('repeat') && repeat
					res2(k)=polyline(obj(k).x(j:end),obj(k).y(j:end), ...
												sprintf('%s (%d:%d)',obj(k).label,j,obj(k).len));
				else
					res2(k)=polyline(obj(k).x(j+1:end),obj(k).y(j+1:end), ...
												sprintf('%s (%d:%d)',obj(k).label,j+1,obj(k).len));
				end
			end
		end
		res1=reshape(res1,size(obj)); %return same shape
		if nargout>1
			res2=reshape(res2,size(obj)); %return same shape
		end
	end
	%----------------------------------------------	
	function [res intersections]=upper_envelope(obj,fullinterval)
		%This function computes the upper envelop over the array of polylines
		%It assumes that all grids are sorted or should be sorted, and treats them as sorted
		%By default only the overlapping segments are used for upper envelope calculation
		%When fullinterval=true all polylines are extended to union of intervals
		if numel(obj)==1
			warning 'Upper envelope is meant for an array of polylines'
			res=obj;
			if nargout>1
				intersections=polyline([],[],'intersection points');
			end
			return;
		end
		l=obj.len;%check that x and y are of the same size
		fullinterval=exist('fullinterval');%when second arg given, full interval
		obj=obj(l>0); %disregard all polylines of zero length
		pt=sort(unique([obj.x]));%collect all the x points in sorted row vector
		%interpolate all lines on all points recording the extrapolation cases
		[intr extr]=obj.interpolate(pt);
		if ~fullinterval		
			%disregard points where at least one line is extrapolated
			mask=sum(extr,1)>0;
			intr(:,mask)=[];
			pt(mask)=[];
			n=sum(~mask);%number of point in the overlap region
		else
			%disregard only points where particular lines are extrapolated (full interval!)
			intr(extr)=-Inf;
			n=numel(pt);
		end
		%find lines on the top
		maxintr=repmat(max(intr),size(intr,1),1);
		top=intr==maxintr;
		%build up the upper envelope
		res=polyline(pt(1),maxintr(1,1),'upper envelope');
		if nargout>1
			intersections=polyline([],[],'intersection points');
		end
		k0=find(top(:,1),1,'first');%index of top line
		%loop through all points
		for j=2:n
			k1=find(top(:,j),1,'first');%index of next top line
			if k1~=k0
				%check if there is an intersection point
				%between the lines:
				ln1=k0;ln2=k1; %intersections between these lines
				pt1=pt(j-1);pt2=pt(j); %which lies between these points
				[y1 extr1]=obj(ln1).interpolate([pt1 pt2]);
				[y2 extr2]=obj(ln2).interpolate([pt1 pt2]); %and these function values (maybe extrapolated)
				%check that neither is extrapolated in both points,
				%and that intersection point is inside the interval <= func values are different at the borders
				if all(~[extr1 extr2]) & all(abs(y1-y2)>0)
					%find the intersection point or points 
					while true
						pt3=fzero(@(x) obj(ln2).interpolate(x)-obj(ln1).interpolate(x),[pt1 pt2]);
						pt3f=obj(ln1).interpolate(pt3);
						%check if there are lines above the found intersection
						[intr2 exrt2]=obj.interpolate(pt3);%interpolate all lines in the new point
						intr2(exrt2)=-Inf; %disregard the extrapolated points
						maxintr2=repmat(max(intr2),size(intr2,1),1);
						ln3=find(intr2==maxintr2,1,'first');
						if ln3==ln1 | ln3==ln2
							%there are no other functions above!
							%add the intersection point
							res=res.inset(pt3,pt3f,res.len);
							if nargout>1
								intersections=intersections.inset(pt3,pt3f); %inset in the end
							end
							%maybe there are some more intersections before next point?
							if ln2==k1
								%no actually, because the left line is the one we started with
								break;
							else
								%indeed, so update the interval of new search
								ln1=ln2;
								pt1=pt3;
								ln2=k1;
								pt2=pt(j);
							end
						else
							%there is line ln3 above the found intersection point
							%so, it is not on the upper envelope
							%need to search again
							ln2=ln3;%new candidate
							pt2=pt3;%new border
						end
					end
				end
			end
			if any(abs(obj(k1).x-pt(j))<eps) || j==n
			% if ismember(pt(j),obj(k1).x) || j==n
				%add new grid point to the end
				res=res.inset(pt(j),maxintr(1,j)); %inset in the end
			end
			k0=k1;%next step
		end
	end
	%----------------------------------------------	
	function [res indxremoved newdots]=secondary_envelope(obj)
		%this function computes the secondary envelope of the polyline
		%returns cleaned polylines, indexes of removed points and new points as a polyline
		%in case of many polylines, clean up each
		%indxremoved is cell array when obj has many elements
		for k=1:numel(obj)
			cur=obj(k);%current line
			%identify the loop-back regions
			ii=cur.x(2:end)>cur.x(1:end-1); %zeros for the loopback regions
			sect=polyline;%sections
			i=1;
			while true
				j=find(ii~=ii(1),1,'first');
				if isempty(j)
					%exit the loop if there are no more loop-backs
					if i>1
						%if sections already started, add the last one
						sect(i)=cur;
					end
					break;
				end
				[sect(i) cur]=cur.chop(j,true); %true to repeat the boundary points
				ii(1:j-1)=[];%chop ii array accordingly
				i=i+1;
			end
			% perform secondary envelope if sections created
			if numel(sect)>1
				sect=sect.sort; %sort all sections since half of them are in opposite direction
				[res(k) newdots(k)]=sect.upper_envelope(true); %true for full interval
				%removed points indexes
				indxremoved{k}=obj(k).diff(res(k),10); %index of dots in obj(k) but not in res(k)
			else
				%without loopbacks -- just copy the input
				res(k)=obj(k);
				indxremoved{k}=[];
				newdots(k)=polyline;
			end
		end %next k 
		%return in same dimensions
		res=reshape(res,size(obj));
		indxremoved=reshape(indxremoved,size(obj));
		if numel(obj)==1
			%for a single polyline return indx as a vector
			indxremoved=indxremoved{1};
		end
	end
	%----------------------------------------------
	function res=plot(obj,ax,varargin)
		%plot the function on the given axes with optional line properties
		%returns line handle(s)
		%additional input arguments - line properties
		%also additional arguments may contain 'sort' to sort the data
		l=obj.len;%check that x and y are of the same size
		if all(l==0)
			warning 'Nothing to plot'
			return
		end
		if nargin>1
			if isempty(ax)
				fig1=figure('Color','white');
				ax=axes('Parent',fig1);
			elseif ~ishandle(ax)
				error 'First argument should be axes handle or []'
			end
		end
		if nargin>2
			mask=cellfun(@(x) ischar(x)&strcmp(x,'sort'),varargin);%find 'sort'
			if sum(mask)>0
				dosort=true;
			else
				dosort=false;
			end
			opts=varargin(~mask);
		else
			opts={'Marker','o','MarkerSize',3,'MarkerFaceColor','auto','MarkerEdgeColor','auto'};
			dosort=false;%no sorting by default
		end
		for k=1:numel(obj) %draw all polylines if there are many
			if l(k)>0
				if dosort %sorting the points of the grid
					[xx j]=sort(obj(k).x);
					xx=reshape(xx,[],1);
					yy=reshape(obj(k).y(j),[],1);
				else
					xx=reshape(obj(k).x,[],1);
					yy=reshape(obj(k).y,[],1);
				end
				if exist('ax')
					hold(ax,'all');
					res(k)=plot(ax,xx,yy,opts{:},'DisplayName',obj(k).label);
				else
					fig1 = figure('Color','white');
					ax = axes('Parent',fig1);
					res(k)=plot(ax,xx,yy,opts{:},'DisplayName',obj(k).label);
				end
			end
		end
		res(l==0)=[];%delete null handles
	end
	%----------------------------------------------
	function demo(obj,n,k)
		%This is just a demo for the upper envelope function:
		%computes the upper envelope of a set of randomized polylines
		if ~exist('n')
			n=10; %lines
		end
		if ~exist('k')
			k=10; %points per line
		end
		for i=1:2:n
			a(i)=polyline(sort([0 5+5*rand(1,k) 10]),[0 rand(1,k+1).*linspace(.5,10,k+1)],sprintf('%d',i));
			a(i+1)=polyline([0 10],[10-i -(10-i)*10/i+10-i],sprintf('%d',i));
		end
		h=a.plot;
		ax=get(h(1),'Parent');
		[ue ix]=a.upper_envelope;
		ue.plot(ax,'LineWidth',1,'Color','red','LineWidth',2);
		ix.plot(ax,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],...
		    				'MarkerSize',4,'Marker','o','LineStyle','none');
		set(ax,'Ylim',[0 10]);
	end
end %methods
end %of classdef

