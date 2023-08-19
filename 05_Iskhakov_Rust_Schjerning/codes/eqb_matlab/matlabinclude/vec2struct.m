function pstruct = vec2struct(pvec,pnames, pstruct,pnames_map);
    %VEC2STRUCT: Procedure to convert a vector to a structure with scalar fields. 
    %
    %SYNTAX:      
    %   pstruct = vec2struct(pvec,pnames, pstruct)
    %
    %OUTPUTS:
    %   pstruct:   Structure with fields taking scalar values (at least one for each value element in pnames) 
    %
    %INPUTS
    %   pvec:      k dimentional vector of scalars
    %
    %   pnames:    k dimensional cell array, that holds names of the fields corresponding to each row in pvec
    %
    %   pstruct:   
    %              If this argument is passed, the fields pstruct corresponding to pnames will be be updated 
    %              with the values in pvec (if they exist already) or added if they do not already exists. 
    %
    %   pnames_map:    result of 'cellfun(@(x) numel(pstruct.(x)),pnames)'
    %
    % See also: struct2vec

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
                    pstruct.(char(pnames(i)))(r, c)=pvec(j,1);
                    j = j + 1; 
                end
            end    
        end

    else

        %new version
        if nargin < 4;
            pnames_map=cellfun(@(x) numel(pstruct.(x)),pnames); %number of parameters under each pname
        end
        j=0;
        for i=1:numel(pnames_map)
            sz=size(pstruct.(char(pnames(i))));
            if iscell(pstruct.(char(pnames(i))))
                pstruct.(char(pnames(i)))=num2cell(reshape(pvec(j+1:j+pnames_map(i)),sz)); %by column read from the vector
            else
                pstruct.(char(pnames(i)))=reshape(pvec(j+1:j+pnames_map(i)),sz); %by column read from the vector
            end
            j=j+pnames_map(i);
        end
    end
        
end % end of vec2struct
