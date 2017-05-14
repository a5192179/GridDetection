close all
a=1:4;
b=2:5;
figure1=figure();
plot(a)
legend(['a=',num2str(a(1),'%d')])
hold on
a=3:6;
figure2=figure();
figure(figure2)
plot(b)
figure(figure1)
plot(a)
legend(['a=',num2str(a(1),'%d')],['a=',num2str(a(2),'%d')])