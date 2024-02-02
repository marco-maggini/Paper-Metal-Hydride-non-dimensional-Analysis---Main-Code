function [value, isterminal, direction] = absorptionCondition(t,y)

data=inputData();
H2_abs=sum(y(1:data.MH_nodes).*data.m_s);

% value = round(H2_abs/sum(data.m_s)-1.05,4);
value = round( mean(y(1:data.MH_nodes)) - 0.995, 4);
% value = round( y(20) - 1.0 ,4);
isterminal = 1;
direction = 1;

end