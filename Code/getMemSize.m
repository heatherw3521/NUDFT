function [ bytes ] = getMemSize( variable )
    
    props = properties(variable); 
    if size(props, 1) < 1
        
        bytes = whos(varname(variable)); 
        bytes = bytes.bytes;
    else
        
        bytes = 0;
        for ii=1:length(props)
            currentProperty = getfield(variable, char(props(ii)));
            pp = props(ii);
            bytes = bytes + getMemSize(currentProperty);
        end    
    
    end   
        
end

function [ name ] = varname( ~ )
name = inputname(1);
end