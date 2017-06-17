function [ j ] = movecell( i, dir )
% compute next cell j given initial cell i and direction dir
% if movement is invalid, return j = -1
j = -1;
switch dir
    case 0
        j = i;
    case 1
        if mod(i,3) == 0
            j = -1;
        else
            j = i+1;
        end
    case 2
        j = i-3;
        if j<= 0
            j = -1;
        end
    case 3
        if mod(i-1,3) == 0
            j = -1;
        else
            j = i-1;
        end
    case 4
        j = i+3;
        if j> 9
            j = -1;
        end
end

return

end

