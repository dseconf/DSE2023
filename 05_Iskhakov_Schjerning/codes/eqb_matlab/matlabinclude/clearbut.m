function clearbut(varargin)
	%This function clears all variables in the workspace except those listed in the input arguments
	vars={};
	% evalin('base','assignin(''caller'',''vars'',who)');
  vars=evalin('base','who');
	for v=reshape(vars,1,[])
		if ~ismember(v,varargin);
			evalin('base',['clear(''' v{1} ''',''global'')']);
		end
	end
end %function