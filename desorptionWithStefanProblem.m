function [t, y] = desorptionWithStefanProblem(d, iniC, options, time)

    function dy = desorptionODE(t,y)
        
        dy = zeros(2*d.MH_nodes+1+d.N,1);

        %ODE system
        P = y(d.MH_nodes+1)*d.R*mean([y(d.MH_nodes+2:2*d.MH_nodes); y(2*d.MH_nodes+2)])/((1+d.fV)*sum(d.V)*d.MW_H2);
        if P < d.Pout
            P = d.Pout;
        end

        Peq = ones(d.MH_nodes,1);
        for i = 1:d.MH_nodes-1
            Peq(i) = d.P0*exp( -d.deltaH_d./(d.R*y(i+d.MH_nodes+1)) + d.deltaS_d/d.R + 0.13*(y(i)-0.5) );
        end
           Peq(d.MH_nodes) = d.P0*exp( -d.deltaH_d./(d.R*y(2*d.MH_nodes+2)) + d.deltaS_d/d.R + 0.13*(y(d.MH_nodes)-0.5) );
        
        criticPressure = P*(2/2.4)^(1.4/0.4); 
       
        if criticPressure > d.Pout
            criticT = mean([y(d.MH_nodes+2:2*d.MH_nodes); y(2*d.MH_nodes+2)])*(2/2.4);
            criticDensity = criticPressure/(4125*criticT);
            fH_out = criticDensity*d.Svalve*sqrt(1.4*criticT*4125);
        else
            fH_out = d.Svalve/(4125*mean([y(d.MH_nodes+2:2*d.MH_nodes); y(2*d.MH_nodes+2)]))*P^(0.4/1.4)*(d.Pout)^(1/1.4)*...
                    sqrt( 2.8/0.4*4125*mean([y(d.MH_nodes+2:2*d.MH_nodes); y(2*d.MH_nodes+2)])*(1 - (d.Pout/P)^(0.4/1.4)) );
        end
            
        for i = 1:d.MH_nodes-1
            if P < Peq(i)
                dy(i) = d.Cd*exp( -d.Ed/(d.R*y(i+d.MH_nodes+1)) )*( (P-Peq(i))/Peq(i) )*y(i);
            else
                dy(i) = 0;
            end
        end
        if P < Peq(d.MH_nodes)
            dy(d.MH_nodes) = d.Cd*exp( -d.Ed/(d.R*y(2*d.MH_nodes+2)) )*( (P-Peq(d.MH_nodes))/Peq(d.MH_nodes) )*y(d.MH_nodes);
        else
            dy(d.MH_nodes) = 0;
        end

        dy(d.MH_nodes+1) = -fH_out - sum( dy(1:d.MH_nodes).*d.m_s*d.SC*d.MW_H2/d.MW_MH );

        Qsource = (d.deltaH_d*dy(1:d.MH_nodes).*d.porosity*d.rhoMH*d.SC/d.MW_MH);
        Qsource = Qsource/( d.capEff );


        dy(d.MH_nodes+2) = (2*d.alfaEff/d.deltaR^2)*y(d.MH_nodes+2+1) - ...
                (2*d.alfaEff/d.deltaR^2)*y(d.MH_nodes+2) + ...
                Qsource(1);

        for i = d.MH_nodes+2+1:2*d.MH_nodes-1
            ri = (i-(d.MH_nodes+2))*d.deltaR;
            dy(i) = (d.alfaEff/(2*d.deltaR*ri) + d.alfaEff/d.deltaR^2)*y(i+1) - ...
                (2*d.alfaEff/d.deltaR^2)*y(i) + ...
                (d.alfaEff/d.deltaR^2 - d.alfaEff/(2*d.deltaR*ri))*y(i-1) + ...
                Qsource(i-(d.MH_nodes+1));
        end
        ri = ri+d.deltaR;
        dy(2*d.MH_nodes) = (d.alfaEff/(2*d.deltaR*ri) + d.alfaEff/d.deltaR^2)*y(2*d.MH_nodes+2) - ...
                (2*d.alfaEff/d.deltaR^2)*y(2*d.MH_nodes) + ...
                (d.alfaEff/d.deltaR^2 - d.alfaEff/(2*d.deltaR*ri))*y(2*d.MH_nodes-1) + ...
                Qsource(2*d.MH_nodes-(d.MH_nodes+1));
        
        %y(101) == epsilon
        deltarS = ( y(2*d.MH_nodes+1) - d.r0 )/(d.r-1);
        deltarL = (d.E-y(2*d.MH_nodes+1))/(d.N-d.r);
        dy(2*d.MH_nodes+1) = 1/(d.rhoPCM*d.lambdaPCM)*...
            ( d.kS*(y(d.rreal-2)-4*y(d.rreal-1))/(2*deltarS) + ...
            d.kL*(y(d.rreal+2)-4*y(d.rreal+1))/(2*deltarL) );

        %y(102) == temperatura interfaccia (salto in i)
        ri = ri+d.deltaR;
        C = (ri+deltarS)/ri*d.kS/d.kMH*(y(2*d.MH_nodes+4)+d.Tpcm-y(2*d.MH_nodes+2))/(2*deltarS);
        TNplus = y(2*d.MH_nodes)+(d.deltaR+deltarS)*C;
        dy(2*d.MH_nodes+2) = d.alfaEff*...
            ( (TNplus-y(2*d.MH_nodes))/(2*ri*d.deltaR) + ...
            (y(2*d.MH_nodes)-2*y(2*d.MH_nodes+2)+TNplus)/d.deltaR^2 ) +...
            Qsource(end);
        T0 = y(2*d.MH_nodes+2) - d.Tpcm;

        Ql = 0;
               
        dy(2*d.MH_nodes+3) = (d.r0 + deltarS)/y(2*d.MH_nodes+1)*(y(2*d.MH_nodes+4)-T0)/(2*deltarS)*dy(2*d.MH_nodes+1) + ...
                d.alfaS*( (y(2*d.MH_nodes+4)-T0)/(2*deltarS*(d.r0+deltarS)) + ...
                (T0-2*y(2*d.MH_nodes+3)+y(2*d.MH_nodes+4))/deltarS^2 );
        
        for i = 2*d.MH_nodes+4:d.rreal-1
            dy(i) = (d.r0 + (i-(2*d.MH_nodes+2))*deltarS)/y(2*d.MH_nodes+1)*(y(i+1)-y(i-1))/(2*deltarS)*dy(2*d.MH_nodes+1) + ...
                d.alfaS*( (y(i+1)-y(i-1))/(2*deltarS*(d.r0+(i-(2*d.MH_nodes+2))*deltarS)) + ...
                (y(i-1)-2*y(i)+y(i+1))/deltarS^2 );
        end

        dy(d.rreal) = 0;

        for i = d.rreal+1:d.Nreal-1
            dy(i) = (d.E-(d.r0+(i-(2*d.MH_nodes+2))*deltarL))/(d.E-y(2*d.MH_nodes+1))*(y(i+1)-y(i-1))/(2*deltarL)*dy(2*d.MH_nodes+1) + ...
                d.alfaL*( (y(i+1)-y(i-1))/(2*deltarL*(d.r0+(i-(2*d.MH_nodes+2))*deltarL)) + ...
                (y(i-1)-2*y(i)+y(i+1))/deltarL^2 );
        end

        dy(d.Nreal) = 2*d.alfaL/deltarL*...
            ( (y(d.Nreal-1)-y(d.Nreal))/deltarL + Ql/d.kL);
        
    end

    [t, y] = ode15s(@desorptionODE, [1e-3 3600e10], iniC, options);

end
