function [u,iout]= fcn(x,sigma)
    persistent uk xk i;
    k1=zeros(1,2);
    A=zeros(2);
    B=zeros(2,1);
    P=zeros(2);
    Q=zeros(2);
    coder.extrinsic('evalin');
    A = evalin('base', 'A'); 
    B = evalin('base', 'B'); 
    P = evalin('base', 'P'); 
    Q = evalin('base', 'Q'); 
    k1 = evalin('base', 'k1'); 
    
    
    if isempty(uk)||isempty(xk)
        i=0;
        xk = x;
        uk=-k1*x;
    end

    if [x' xk']*[A'*P+P*A+sigma*Q -P*B*k1;-(B*k1)'*P zeros(2)]*[x;xk]>=0
        i = i+1;
        uk = -k1*x;
        xk = x;
    end

    u=uk;
    iout=i;
end
