function definePars(deltaT, alpha)
    global Z phi psiZ psiW trans 
    
    % phi
    phiTilde = [1, deltaT, deltaT^2/2; 0, 1, deltaT; 0, 0, alpha];
    phi = [phiTilde, zeros(3); zeros(3), phiTilde];
    
    % psiZ
    psiTildeZ = [deltaT^2/2; deltaT; 0];
    psiZ = [psiTildeZ, zeros(3,1); zeros(3,1), psiTildeZ];
    
    % psiW
    psiTildeW = [deltaT^2/2; deltaT; 1];
    psiW = [psiTildeW, zeros(3,1); zeros(3,1), psiTildeW];
    
    % transition prob matrix
    m = 5;
    trans = ones(m)/20;
    trans(1:m+1:m*m) = 16/20;
    
    % Z vector
    Z = [0, 3.5, 0, 0, -3.5; 0, 0, 3.5, -3.5, 0];
