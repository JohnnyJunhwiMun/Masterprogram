function [F, G] = matrices_q4_ag(A, B, h, tau)
    if tau < h
        Fx = expm(A*h);
        Fu = (expm(A*h) - expm(A*(h-tau)))/A * B;
        G1 = (expm(A*(h-tau)) - eye(size(A)))/A * B;
        F = [Fx, Fu, zeros(2,1); zeros(1,4); zeros(1,2), 1, 0];
        G = [G1; 1; 0];
    elseif tau <= 2*h
        Fxh = expm(A*h);
        Fuh = (expm(A*h) - expm(A*(2*h-tau)))/A * B;
        Gh = (expm(A*(2*h-tau)) - eye(size(A)))/A * B;
        F = [Fxh, Fuh, Gh; zeros(size(Fxh)),[0;1],zeros(2,1)]; 
        G = [zeros(size(Fxh(:,1))); 1; 0];
    else
        F = [];
        G = [];
    end
end