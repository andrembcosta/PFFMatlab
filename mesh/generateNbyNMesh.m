function mesh = generateNbyNMesh(n, size, cracked)

    mesh.numele = n*n;
    if cracked
        mesh.numnod = (n+1)*(n+1) + floor((n+1)/2);
    else
        mesh.numnod = (n+1)*(n+1);
    end
    mesh.xnod = zeros(mesh.numnod, 4);
    mesh.ien = zeros(mesh.numele, 5);
    mesh.ienb = [];
    mesh.physid = ones(mesh.numele,1);
    
    for i = 1:mesh.numnod
        
        mesh.xnod(i,:) = [i getNodeCoords(n,i)];
        
    end
    
    for j = 1:mesh.numele
        
        mesh.ien(j,:)= [j getElementNodes(n,j,cracked)];
        
    end
    
    mesh.xnod(:,2:end) = mesh.xnod(:,2:end)/(n/2);
    
    mesh.xnod(:,2:end) = (size/2)*mesh.xnod(:,2:end);
    
end