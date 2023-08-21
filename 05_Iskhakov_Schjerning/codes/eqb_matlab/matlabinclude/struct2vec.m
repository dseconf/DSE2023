function pvec = struct2vec(pstruct,pnames,pnames_map)
    %STRUCT2VEC: Procedure that extract vector from structure with matrix fields. 
    %
    % SYNTAX:
    %   pvec = struct2vec(pstruct,pnames)
    % 
    % OUTPUTS
    %   pvec:      k dimentional vector of scalars
    %
    % INPUTS
    %   pstruct:   Structure with fields taking scalar values (at least one field for each value element in pnames) 
    %
    %   pnames:    k dimensional cell array, that specify names of the fields in pstruct that should be written to pvec
    %
    %   pnames_map:    result of 'cellfun(@(x) numel(pstruct.(x)),pnames)'
    %
    % See also: vec2struct
    if nargin < 2;
        pnames=fieldnames(pstruct);
        m=cellfun(@(x) isnumeric(pstruct.(x)),pnames);
        pnames(~m)=[]; %drop non-numeric fields
    end

    version=2;

    if version==1

        % old version
        k=numel(pnames);
        j = 1;
        for i=1:k;
            r_i = size(pstruct.(char(pnames(i))),1);
            c_i = size(pstruct.(char(pnames(i))),2);
            for c=1:c_i;
                for r=1:r_i;
                    pvec(j,1)=pstruct.(char(pnames(i)))(r, c);
                    j = j + 1; 
                end
            end    
        end

    else

        %new version
        if nargin < 3;
            pnames_map=cellfun(@(x) numel(pstruct.(x)),pnames); %number of parameters under each pname
        end
        pvec=zeros(sum(pnames_map),1);  %preallocate memory, column vector
        j=0;
        for i=1:numel(pnames_map)
            if iscell(pstruct.(char(pnames(i))))
                pvec(j+1:j+pnames_map(i))=[pstruct.(char(pnames(i))){:}]; %by column fit into the vector
            else
                pvec(j+1:j+pnames_map(i))=pstruct.(char(pnames(i)))(:); %by column fit into the vector
            end
            j=j+pnames_map(i);
        end
    
    end
end % end of struct2vec
