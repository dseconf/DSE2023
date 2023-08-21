function p1 = combinestructs(varargin)
    %STRUCT2VEC: Procedure that extract vector from structure with matrix fields. 
    %
    % SYNTAX:
    %   pvec = struct2vec(struct1, struct1, struct3)
    % 
    % INPUTS
    %   pstruct:   Structure with fields taking scalar values (at least one field for each value element in pnames) 
    %
    %   pnames:    k dimensional cell array, that specify names of the fields in pstruct that should be written to pvec
    %
    % See also: vec2struct

    p1=varargin{1};
    for j=2:nargin; 
        p2=varargin{j};
        pnames=fieldnames(p2);       
        k=numel(pnames);
        for i=1:k;
            p1.(char(pnames(i)))=p2.(char(pnames(i)));
        end
    end
end % end of struct2vec
