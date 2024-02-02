function results = getResults (y, d, time, nclass)
%%it yields the main output parameters of the calculation, namely mMH/ms,
%%H2g, Temperature, Pressure, Equilibrium Pressure, reaction rate, by
%%evaluating input data.

results = struct(...
    'MassesRatio', y(:,1:d.MH_nodes), ...
    'H2_g', y(:,d.MH_nodes+1), ...
    'Temperature', y(:,d.MH_nodes+2:end), ...
    'P', y(:,d.MH_nodes+1)*d.R.*mean(y(:,d.MH_nodes+2:2*d.MH_nodes+1), 2)/((1+d.fV)*sum(d.V)*d.MW_H2) ...
    );

if strcmp(nclass,'abs')
    
    results.Peq = ...
        d.P0*exp(d.deltaH_a./(d.R*y(:,d.MH_nodes+2:end)) - ...
        d.deltaS_a/d.R + d.sl*(y(:,1:d.MH_nodes) - 0.5));
    
    for i = 1:d.MH_nodes
       results.r(:,i) = ...
          d.Ca*exp(-d.Ea./(d.R*y(:,d.MH_nodes+2+i-1))).*...
          log(results.P./results.Peq(:,i)).*(1-y(:,i));
    end

    results.Ra = 9.81*d.beta*...
        abs(y(:,end)-d.Tpcm)...
        *d.L^3*d.Pr./d.ni^2;
    St = d.cpPCM*(y(:,end) - d.Tpcm)/d.lambdaPCM;
    Fo = d.alfaL*time/d.L^2;
    theta = St.*Fo;
    results.Nu = (2*theta).^-0.5 + (0.35*results.Ra.^0.25 - (2*theta).^-0.5).*(1+(0.0175*results.Ra.^.75.*theta.^(3/2)).^-2).^(1/-2);
%     results.Nu = (0.35*results.Ra.^.25)./(1+(.143/d.Pr)^(9/16))^(4/9);
    results.h = results.Nu*d.kL/d.L;
    results.Bi = (results.h*d.D/4)/d.kEff;

    results.Qreact = -d.deltaH_a*results.r*d.m_s*d.SC/d.MW_MH;
    results.Qexch = pi*d.D*d.L*results.h.*...
        (d.Tpcm-results.Temperature(:,end));
  

elseif strcmp(nclass, 'drm')
    results.Peq = nan(size(y,1), 1);
    results.r   = nan(size(y,1), 1);
    
elseif strcmp(nclass, 'des')
    
    results.Peq = ...
        d.P0*exp( -d.deltaH_d./(d.R*y(:,d.MH_nodes+2:2*d.MH_nodes+1)) + ...
        d.deltaS_d/d.R + 0.13*(y(:,1:d.MH_nodes)-0.5) );
    
    pressuresRatio = zeros(size(y,1), d.MH_nodes);
    for i = 0:d.MH_nodes-1
        activeReaction = results.P < results.Peq(i+1);
        cutIdx = length(activeReaction)-sum(activeReaction);
            pressuresRatio(:,i+1) = ((y(:,d.MH_nodes+1)*d.R.*mean(y(:,d.MH_nodes+2:2*d.MH_nodes+1), 2)/((1+d.fV)*sum(d.V)*d.MW_H2)) - ...
                    (d.P0*exp( -d.deltaH_d./(d.R*y(:,d.MH_nodes+2+i)) + d.deltaS_d/d.R + 0.13*(y(:,i+1)-0.5) )))...
                    ./(d.P0*exp( -d.deltaH_d./(d.R*y(:,d.MH_nodes+2+1)) + d.deltaS_d/d.R + 0.13*(y(:,i+1)-0.5) ));
            results.r(:,i+1) = ...
                d.Cd*exp( -d.Ed./(d.R*y(:,d.MH_nodes+2+i)) ).*...
                ( pressuresRatio(:,i+1) ).*y(:,i+1);
    end
    results.r(1:cutIdx,:) = 0;
    if size(y,2)>2*d.MH_nodes+1
        results.solidFront = y(:,2*d.MH_nodes+2);
        results.deltaS = ( y(:,2*d.MH_nodes+2) - d.r0 )/(d.r-1);
        results.Qexch = d.kS*pi*(d.D+2*results.deltaS)*d.L.*(y(:,2*d.MH_nodes+4)+d.Tpcm-y(:,2*d.MH_nodes+1))./(2*results.deltaS);
        results.h_star = results.Qexch./(d.A*( y(:,2*d.MH_nodes+4)+d.Tpcm-y(:,2*d.MH_nodes+1) ));
        results.Bi_star = results.h_star*d.D/(4*d.kEff);
        results.Nu_star = results.h_star*d.L/d.kL;
    end
    results.Qreact = d.deltaH_d*results.r*d.m_s*d.SC/d.MW_MH;

end
end