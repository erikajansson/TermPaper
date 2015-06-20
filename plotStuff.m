% plot option price
function plotStuff(i, S1, S2, I1,I2, BRCnohitnocall, BRCnohit, BRChitnocall, BRChit, Diffnohit, Diffhit)
figure(i)
subplot(3,1,1)
surf(S1,S2,BRCnohitnocall(I1,I2))
hold on
mesh(S1,S2,BRCnohit(I1,I2))
hold off
box on
grid on
legend('Never called','Called optimally')
xlabel('s_1')
ylabel('s_2')
zlabel('Price','FontSize',14)
title('Barrier not yet hit')

%figure(i*3 + 5)
subplot(3,1,2)
surf(S1,S2,BRChitnocall(I1,I2))
hold on
mesh(S1,S2,BRChit(I1,I2))
hold off
legend('Never called','Called optimally')
title('Hit already, called never vs called optimally')

%figure(i*3 + 4)
subplot(3,1,3)
surf(S1,S2,Diffnohit(I1,I2))
hold on
mesh(S1,S2,Diffhit(I1,I2))
hold off
legend('No hit so far','Barrier already hit')
%set(h,'FontSize',14)
axis on
xlabel('s_1')
ylabel('s_2')
zlabel('Price Difference','FontSize',14)
title('Optimal Call Strategy')
end  
