clear;
clc;
close all;

iterfem = [ ...
  35, 33, 28, 25, 24, 22;
  44, 39, 33, 31, 31, 28;
  55, 46, 42, 42, 41, 37;
  64, 56, 55, 55, 54, 48;
  78, 76, 76, 75, 71, 59;
];

iternu = [ ...
  23, 41, 188, 300, 300, 300;
  38, 46, 223, 300, 300, 300;
  53, 49, 240, 300, 300, 300;
  79, 64, 288, 300, 300, 300;
  108, 87, 300, 300, 300, 300;
];

iterfem = log10(iterfem);
iternu = log10(iternu);
nen = [32,64,128,256,512,1024];
hh = nen./nen./nen;
loghh = log10(hh);
figure(1)
hold on
plot(loghh,iternu(1,:),'-o','LineWidth',2);
plot(loghh,iternu(2,:),'-o','LineWidth',2);
plot(loghh,iternu(3,:),'-o','LineWidth',2);
plot(loghh,iternu(4,:),'-o','LineWidth',2);
plot(loghh,iternu(5,:),'-o','LineWidth',2);
xlabel("log10(h)","FontSize",14,"FontWeight","bold");
ylabel("log10(iternum)","FontSize",14,"FontWeight","bold");
legend("p=3","p=4","p=5","p=6","p=7","FontSize",14,"Location","best");
hold off

figure(2)
hold on
plot(loghh,iterfem(1,:),'-o','LineWidth',2);
plot(loghh,iterfem(2,:),'-o','LineWidth',2);
plot(loghh,iterfem(3,:),'-o','LineWidth',2);
plot(loghh,iterfem(4,:),'-o','LineWidth',2);
plot(loghh,iterfem(5,:),'-o','LineWidth',2);
xlabel("log10(h)","FontSize",14,"FontWeight","bold");
ylabel("log10(iternum)","FontSize",14,"FontWeight","bold");
legend("p=3","p=4","p=5","p=6","p=7","FontSize",14,"Location","best");
hold off