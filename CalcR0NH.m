function betaE=CalcR0NH(R0E,N,sigma,gamma,M)
    betav=linspace(0.01,0.1,1001); % Span to search for beta
    R0=zeros(1001,1); % Record R0

    for rr=1:1001
        beta=betav(rr); % set the beta value
        F=[];
        for ii=1:4
            F=[F; 0 0 0 0 beta*M(ii,1)*N(ii)/N(1)  beta*M(ii,2)*N(ii)/N(2) beta*M(ii,3)*N(ii)/N(3) beta*M(ii,4)*N(ii)/N(4)]; 
        end

        F=[F; zeros(4,length(F(1,:)))];

        V=zeros(size(F));
        for ii=1:4
            V(ii,ii)=sigma;
            V(4+ii,ii)=-sigma;
        end

        for ii=5:8
            V(ii,ii)=gamma;
        end

        test=abs(eig(F*inv(V))); % Calculate based on the next generation matrix
        R0(rr)=max(test);
    end
    betaE=pchip(R0,betav,R0E);
end

  