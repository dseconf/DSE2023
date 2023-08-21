function varargout=combinevecs(varargin)
    % build matrix of grids from several grids given in order
    % OUTPUT: combined vector
    %         matrix of masks where each row is mask for the 
    %         last grid by all values of the rest of the grids
    for i=nargin:-1:1
        grid=reshape(varargin{i},[],1);
        if i==nargin
            varargout{1}=grid; %grid in column
            % if nargout>1
            %     varargout{2}=(ones(size(varargout{1}))==1);
            % end
        else
            varargout{1}=[kron(grid,ones(mm,1)) repmat(varargout{1},numel(grid),1)];
            % if i==1 && nargout>1
            %     maskgen=varargout{1}(:,1:end-1)
            %     maskgenu=unique(maskgen,'rows')
            %     varargout{2}=[];
            %     for i=1:size(maskgenu,1)
            %         varargout{2}=[sum(maskgen==ones(size(maskgen,1),1)*maskgenu(i,:),2)==size(maskgen,2) varargout{2}];
            %     end
            % end
        end
        mm=size(varargout{1},1);
    end
    % mask for the last grid totals
    if nargout>1
        mm=mm/numel(varargin{nargin});
        varargout{2}=kron(eye(mm),ones(numel(varargin{nargin}),1));
    end
end