P = [
    100
    200
    300
    400
    500
    600
    700
    800
    900
    1000
    ];

h = P*-0.010197/100

Q1 = [
    -3.14099e-8
    -6.28199e-8
    -9.42298e-8
    -1.2564e-7
    -1.5705e-7
    -1.8846e-7
    -2.19869e-7
    -2.51279e-7
    -2.82689e-7
    -3.14099e-7
    ];
Q2 = [
    -1.5705e-8
    -3.14099e-8
    -4.71149e-8
    -6.28199e-8
    -7.85248e-8
    -9.42298e-8
    -1.09935e-7
    -1.2564e-7
    -1.41345e-7
    -1.5705e-7
    ];

A = 0.03141;

q1 = Q1/A*100;
q2 = Q2/A*100;

p1 = polyfit(h,q1,1)
y1 = p1(1)*h + p1(2)
p2 = polyfit(h,q2,1)
y2 = p2(1)*h + p2(2)

figure
plot(h,q1,'bo',h,y1,'b',h,q2,'ko',h,y2,'k')
xlabel('Head gradient (cm H_2O/cm)')
ylabel('Darcy flux q (cm/s)')
legend('Darcy1','Darcy1 fit','Darcy2','Darcy2 fit')
set (gca,'Ydir','reverse')
set (gca,'Xdir','reverse')
