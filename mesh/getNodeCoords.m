function coords = getNodeCoords(n,k)

    coords = zeros(1,3);
    
    %get X
    %find column
    column = rem(k,n+1);
    if column == 0
        column = n+1;
    end
    coords(1,1) = column - 1 - n/2;
    %get Y
    %find row
    row = ceil(k/(n+1));
    if row > n+1
        %cracked node
        coords(1,2) = 0;
    else
        coords(1,2) = -row + 1 + n/2;
    end
    
end
