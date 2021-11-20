function f=Modular_Period_ObsvPlot(t,V,errors,P,loop)
% P; selected period; 
% loop: 
t_modu=rem(t,P);
figure,
hold on;
for i=1:loop
errorbar(t_modu+(i-1)*P,V,errors);
end
hold off;