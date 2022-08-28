function elem = getElementNodes(n,k,cracked)

    row = ceil(k/n) - 1;

    %elem = [k+row, k+row+1, k+row+n+2, k+row+n+1];
    elem = [k+row, k+row+n+1, k+row+n+2, k+row+1];
    
    if cracked
    
        %crackedElems = (n^2)/2+1:1:(n^2)/2+n/2;

        %if ismember(k,crackedElems)
        if k >= (n^2)/2+1

            elem(1) = elem(1) + (n+1)*(n/2 + 1);
            elem(4) = elem(4) + (n+1)*(n/2 + 1);

        end

        if k == (n^2)/2+n/2
            elem(4) = elem(4) - (n+1)*(n/2 + 1);
        end
    
    end
     
end