function [values, terminal, direction] = desorptionStefanCondition (t,y)

data = inputData();

values1 = max( round(y(1:data.MH_nodes), 2) );
values2 = -round( y(2*data.MH_nodes+1) - data.E, 3);
values3 = round( t - data.timeLimit, 0);

values = [values1; values2; values3];
% values = values1;
terminal = [1; 1; 1];
direction = [-1; -1; 1];

end