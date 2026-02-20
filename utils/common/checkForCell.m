function out = checkForCell(in,variable_type)
if nargin<2
    variable_type = 'numeric'; 
else
    if islogical(variable_type) %backward compatibility
        if variable_type==true
            variable_type = 'logical'; 
        else
            variable_type = 'numeric'; 
        end 
    end
end

if iscell(in)
       out = [];
       for it = 1:length(in)
           switch variable_type
               case 'logical'
                currval = ismember(in{it},{'True' 'TRUE'});
               case 'numeric'
                currval = str2num(in{it}); 
               case 'string'
                currval = in{it}; 
            end
           if isempty(currval), currval = NaN; end
           out(it,1:length(currval)) = currval;
       end
else
    out = in;
end