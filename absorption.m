function [t, y, absorptionEnd] = absorption(d, iniC, options)

    function dy = absorptionODE (t, y)
        dy = zeros(2*d.MH_nodes+1,1);

        P = y(d.MH_nodes+1)*d.R*mean(y(d.MH_nodes+2:end))/(d.MW_H2*sum(d.V)*(1+d.fV)); %25 % free space
        if P > d.Pin
            P = d.Pin;
        end
%         d.rhoH_foam = d.porosity*d.foam_porosity*d.MW_H2*P/(d.R*mean(y(d.MH_nodes+2:2*d.MH_nodes+1)));
%         d.capEff   = d.rho_foam*d.cp_foam + (d.rhoM + d.rhoH)*d.cpMH + d.rhoH_foam*d.cpH;
%         d.alfaEff  = d.kEff/d.capEff;

        Peq = d.P0*exp(d.deltaH_a./(d.R*y(d.MH_nodes+2:end)) - d.deltaS_a/d.R + d.sl*(y(1:d.MH_nodes) - 0.5));
        
        criticPressure = d.Pin*(2/2.4)^(1.4/0.4);
        if criticPressure > P
            criticT = d.Tin*(2/2.4);
            criticDensity = criticPressure/(4125*criticT);
            fH_in = criticDensity*d.Svalve*sqrt(1.4*criticT*4125);
        else
            fH_in = d.Svalve/(4125*d.Tin)*d.Pin^(0.4/1.4)*(P)^(1/1.4)*...
                        sqrt( 2.8/0.4*4125*d.Tin*(1 - (P/d.Pin)^(0.4/1.4)) );
        end
        
        Ra = 9.81*d.beta*abs(d.Tpcm-y(end))*d.L^3*d.Pr/d.ni^2;
%         Nu = (0.35*Ra^.25)/(1+(.143/d.Pr)^(9/16))^(4/9);
        St = d.cpPCM*( abs(y(end) - d.Tpcm) )/d.lambdaPCM;
        Fo = d.alfaL*t/d.L^2;
        theta = St*Fo;
        Nu = (2*theta)^-0.5 + (0.35*Ra^0.25 - (2*theta)^-0.5)*(1+(0.0175*Ra^.75*theta^(3/2))^-2)^(1/-2);
        h = Nu*d.kL/d.L;
%         h = 2500;

        dy(1:d.MH_nodes) = d.Ca*exp(-d.Ea./(d.R*y(d.MH_nodes+2:end))).*log(P./Peq).*(1-y(1:d.MH_nodes));

        dy(d.MH_nodes+1) = fH_in - ...
            sum( dy(1:d.MH_nodes).*d.m_s*d.wt );

%         Qsource = (-d.deltaH_a*dy(1:d.MH_nodes)*d.rhoMH*d.SC/d.MW_MH)/d.capEff;
        Qsource = (-d.deltaH_a*dy(1:d.MH_nodes).*d.m_s*d.SC/d.MW_MH)./(d.capEff*d.V);

        dy(d.MH_nodes+2) = (2*d.alfaEff/d.deltaR^2)*y(d.MH_nodes+2+1) - ...
                (2*d.alfaEff/d.deltaR^2)*y(d.MH_nodes+2) + ...
                Qsource(1);

        for i = d.MH_nodes+2+1:2*d.MH_nodes+1-1
            ri = (i-(d.MH_nodes+2))*d.deltaR;
            dy(i) = (d.alfaEff/(2*d.deltaR*ri) + ...
                d.alfaEff/d.deltaR^2)*y(i+1) - ...
                (2*d.alfaEff/d.deltaR^2)*y(i) + ...
                (d.alfaEff/d.deltaR^2 - d.alfaEff/(2*d.deltaR*ri))*y(i-1) + ...
                Qsource(i-d.MH_nodes-1);
        end

        ri = ri+d.deltaR;
        dy(2*d.MH_nodes+1) = (2*d.alfaEff/d.deltaR^2)*y(2*d.MH_nodes+1-1) - ...
                (d.alfaEff*h/(ri*(d.kEff)) + 2*h*d.alfaEff/((d.kEff)*d.deltaR) + 2*d.alfaEff/d.deltaR^2)*y(2*d.MH_nodes+1) + ...
                (d.alfaEff*h/(ri*(d.kEff)) + 2*h*d.alfaEff/((d.kEff)*d.deltaR))*d.Tpcm + ...
                Qsource(d.MH_nodes);

            
    end

    [t,y,absorptionEnd] = ode15s(@absorptionODE, [1e-3 3600*10], iniC, options);
    
end