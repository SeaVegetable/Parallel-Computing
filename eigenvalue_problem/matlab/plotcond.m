clear;
clc;
close all;

CONDg = [128.772376688697,160.501378779429,210.156024285443,284.115519555739,391.68469533028,546.396439924173,767.249132633978; ...
         609.137128871121,697.062394589085,855.574932666225,1099.35830957425,1462.8924726307,1996.95272402776,2770.20713695242; ...
         3518.413515647,3735.96214884153,4355.49030144981,5335.08101818439,6821.66977271751,9052.02341729086,12338.0596023783; ...
         22975.4971178248,22754.0368485412,25577.9160331892,30213.2175051883,37299.2349380685,48088.0311205785,64238.5432623905; ...
         162178.659973624,148919.688175903,162741.637312242,187040.049723924,224209.359138189,281287.86443683,367839.323346322];
COND  = [28.5604130405391,82.0998641532794,327.208502603839,1307.67226784349,5229.54136662966,20917.0246489625,83666.961192343; ...
         222.54992821228,215.004496441155,381.733159913437,1525.39862784631,6100.11680094581,24399.0175960881,97594.634821238; ...
         1812.5795760605,1700.6315544374,1688.11008936008,1781.51080582718,7124.05845026019,28494.3251064842,113975.429851442; ...
         15132.6441280637,13634.8151004491,13493.1503553076,13505.2718704483,13508.6624082014,32799.6162618335,131196.031260472; ...
         129178.94005403,110155.733696557,108472.867401802,108589.13288144,108621.638748566,108629.550059182,148805.537199058];

NEN = [16,32,64,128,256,512,1024];
PP = [3,4,5,6,7];
HH = NEN./NEN./NEN;

log_HH = log(HH);
log10_HH = log10(HH);
log_COND = log10(COND);
log_COND_FULL_g = log10(CONDg);

figure(1)
hold on
for i = 1 : length(PP)
    plot(log10_HH, log_COND(i,:), '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
    for j = 1 : length(HH)-1
        slope = (log_COND(i,j+1) - log_COND(i,j)) / (log10_HH(j+1)-log10_HH(j));
        x_mid = (log10_HH(j+1) + log10_HH(j)) / 2;
        y_mid = (log_COND(i,j+1) + log_COND(i,j)) / 2;
        text(x_mid, y_mid, sprintf('%.2f', slope), 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    end
end
hold off

legendLabels = arrayfun(@(x) sprintf('p = %.0f', x), PP, 'UniformOutput', false);
legend(legendLabels, 'FontSize', 12,'FontWeight', 'bold', 'Location', 'best');
xlabel("log10(h)","FontSize",12,"FontWeight","bold");
ylabel("log10(kappa iga)","FontSize",12,"FontWeight","bold");

figure(2)
hold on
for i = 1 : length(PP)
    plot(log10_HH, log_COND_FULL_g(i,:), '-o', 'LineWidth', 1.5, 'MarkerSize', 8);
    for j = 1 : length(HH)-1
        slope = (log_COND_FULL_g(i,j+1) - log_COND_FULL_g(i,j)) / (log10_HH(j+1)-log10_HH(j));
        x_mid = (log10_HH(j+1) + log10_HH(j)) / 2;
        y_mid = (log_COND_FULL_g(i,j+1) + log_COND_FULL_g(i,j)) / 2;
        text(x_mid, y_mid, sprintf('%.2f', slope), 'HorizontalAlignment', 'center', 'FontSize', 12, 'FontWeight', 'bold');
    end
end
hold off

legendLabels = arrayfun(@(x) sprintf('p = %.0f', x), PP, 'UniformOutput', false);
legend(legendLabels, 'FontSize', 12,'FontWeight', 'bold', 'Location', 'best');
xlabel("log10(h)","FontSize",12,"FontWeight","bold");
ylabel("log10(kappa full greville)","FontSize",12,"FontWeight","bold");
