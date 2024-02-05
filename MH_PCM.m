clear, clc

%% Data and Initial conditions
d = inputData();

Peq = ( d.P0*exp(d.deltaH_a/(d.R*d.Tpcm) -...
    d.deltaS_a/d.R + d.sl*(0 - 0.5)) )*...
    ones(1,d.MH_nodes);
iniC(1:d.MH_nodes) = 0;
iniC(d.MH_nodes+1) = (1+d.fV)*Peq(1).*sum(d.V)*d.MW_H2/(d.R*d.Tpcm);
iniC(d.MH_nodes+2:2*d.MH_nodes+1) = d.Tpcm+0.01;

options = odeset('RelTol', 1e-10, 'AbsTol', 1e-13,...
    'Events', @absorptionCondition);

%% Absorption
fprintf('Modeling absorption...')

[t_abs, y] = absorption(d, iniC, options);

results.abs = getResults(y, d, t_abs, 'abs');
results.abs.liqFront = sqrt( d.D^2/4 + ...
        (cumtrapz(t_abs, abs(results.abs.Qexch)/d.lambdaPCM)/d.rhoPCM)/(pi*d.L) );
results.abs.fH_in = gradient(results.abs.H2_g)+sum( results.abs.r*d.m_s*d.wt ,2 );
results.abs.H2_abs = results.abs.MassesRatio*d.m_s*d.wt;
results.abs.phi = (results.abs.H2_abs)/(d.wt*sum(d.m_s));

time = t_abs;

fprintf('\nAbsorption complete. (%.0f min) \n\n', t_abs(end)/60)

input('Do you want to model desorption?\n[y] \n[n]\n', 's');
if ans == 'y'
    
elseif ans == 'n'
    return
else
    fprintf('Wrong input. Continuing with desorption\n')
end

%% Absorption post-processing

 figure; semilogy(time*d.alfaEff/d.L^2, results.abs.fH_in*120e6*d.D/(d.alfaEff^3*d.rhoMH*1e18))
 ylabel('P^*'); xlabel('Fo')
 grid
 
 figure; semilogy(results.abs.phi, results.abs.fH_in*120e6*d.D/(d.alfaEff^3*d.rhoMH*1e18))
 ylabel('P^*'); xlabel('\phi^*')
 grid

 figure; plot(time/60, results.abs.H2_abs*1e3, 'k', 'LineWidth', 1.0)
 xlabel('time [min]')
 ylabel('H_2 absorbed [g]')
 grid on
 
%% Desorption with Stefan Problem
iniC = [...
    ones(1,d.MH_nodes) ...
    results.abs.H2_g(end) ...
    d.Tpcm*ones(1,d.MH_nodes-1) ...
    d.r0+0.001 ...
    d.Tpcm ...
    zeros(1,d.N-1)]; %results.drm.MassesRatio(end,:) & results.drm.Temperature(end,1:end-1) & results.drm.Temperature(end,end)
options = odeset('NonNegative', 1:d.MH_nodes,...
    'RelTol', 1e-12, 'AbsTol', 1e-15, 'Events', @desorptionStefanCondition);

fprintf('Modeling desorption (Stefan)...')

[t_des, y] = desorptionWithStefanProblem(d, iniC, options, time, cycle);
t_des = t_des+time(end);

% swap y(101)<->y(102)
temp = y(:,2*d.MH_nodes+1);
y(:,2*d.MH_nodes+1) = y(:,2*d.MH_nodes+2);
y(:,2*d.MH_nodes+2) = temp;

results.des = getResults(y,d,t_des,'des');

results.des.H2_des = (1-y(:,1:d.MH_nodes))*d.m_s*d.wt;
results.des.phi = 1 - (results.des.H2_des)/(d.wt*sum(d.m_s));

fprintf('\nDesorption complete. (%.0f min)\n\n', (t_des(end)-t_des(1))/60)

for i = 1:length(results.des.P)
    if results.des.P(i) < d.Pout
        results.des.P(i) = d.Pout;
    end
end
critic.P = results.des.P*(2/2.4).^(1.4/0.4);
results.des.fH_out = zeros(length(t_des),1);

for i = 1:length(t_des)
    if critic.P(i) > d.Pout
        critic.T = mean(y(i,d.MH_nodes+2:2*d.MH_nodes+1),2)*(2/2.4);
        critic.D = critic.P(i)/(4125*critic.T);
        results.des.fH_out(i) = critic.D*d.Svalve.*sqrt(1.4*critic.T*4125);
    else
        results.des.fH_out(i) = d.Svalve./( 4125*mean(y(i,d.MH_nodes+2:2*d.MH_nodes+1), 2) ).*...
            results.des.P(i)^(0.4/1.4)*(d.Pout)^(1/1.4).*...
                    sqrt( 2.8/0.4*4125*mean(y(i,d.MH_nodes+2:2*d.MH_nodes+1), 2).*...
                    (1 - (d.Pout./results.des.P(i)).^(0.4/1.4)) );
    end
end

%% global post-processing

PCMmass = (d.wt*sum(d.m_s)*d.deltaH_d/d.MW_H2)/d.lambdaPCM;
results.abs.phi = (results.abs.H2_abs)/(d.wt*sum(d.m_s));
results.des.phi = 1-(results.des.H2_des)/(d.wt*sum(d.m_s));
results.des.Power = results.des.fH_out*120e6;

fprintf(['The storage system has:\n' ...
    'Total mass: %.2f kg\n'...
    'PCM mass: %.2f kg\n'...
    'Total volume: %.3f m3\n'...
    'Global gravimetric density: %.3f %%\n'...
    'Global volume density: %.2f kWh/m3\n'...
    'Total amount of H2 stored: %.2f g\n'...
    ],...
    sum(d.m_s)+PCMmass,...
    PCMmass,...
    pi*results.des.solidFront(end)^2*d.L,...
    results.abs.H2_abs(end)/(sum(d.m_s)+PCMmass)*100,...
    results.abs.H2_abs(end)*33.3333/(pi*results.des.solidFront(end)^2*d.L),...
    results.abs.H2_abs(end)*1e3)

%non-dimensional desorption
figure
plot(results.des.phi, results.des.Power/1e3, 'LineWidth',1.0); grid
set(gca, 'yscale', 'log', 'FontSize',18)
% des_avgPower = trapz(outMatrix{1,1}.des.t, outMatrix{1,1}.des.Power)...
%     /outMatrix{1,1}.des.t(end);
set(gca, 'xdir', 'reverse')

title('Desorption')
xlabel('\phi^*')
ylabel('P [kW]')

% dimensional desorption
figure
plot((t_des-t_des(1))/60,results.des.Power/1e3, 'LineWidth',1.0)
grid
xlabel('Time [min]'); ylabel('Desorption power [kW]')

figure; plot((t_des-t_des(1))/60, results.des.H2_des*1e3, 'k', 'LineWidth', 1.0)
xlabel('time [min]')
ylabel('H_2 desorbed [g]')
grid on


%% clearance
clear Peq iniC options temp critic
