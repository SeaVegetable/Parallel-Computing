clear;
clc;
close all;

p_vals     = [3 4 5 6 7];
nElemX_vec = [16 32 64 128 256 512 1024];
hh = nElemX_vec./nElemX_vec./nElemX_vec;
loghh = log10(hh);

kappa_full_mat = [
     60.5003622,    60.3974355,    60.5468206,    60.5845103,    60.5937839,    60.5960838,    60.5966564;
    413.0758790,   404.6447060,   405.9656530,   406.3525870,   406.4468510,   406.4701100,   406.4758860;
   3000.6016160,  2853.8775180,  2863.1638960,  2866.8753850,  2867.7713000,  2867.9912340,  2868.0457190;
  22373.9648380, 20370.0343240, 20408.2261290, 20441.9063720, 20449.9762340, 20451.9469210, 20452.4338530;
 171498.5060970,146940.7072690,146722.0057250,147016.8908010,147087.2842810,147104.3833670,147108.5974050;
];
log_kappa_full_mat = log10(kappa_full_mat);

kappa_iga_mat = [
   425.061724,   420.378988,   421.754162,   422.270214,   422.398558,   422.430508,   422.438478;
  3361.016769,  3287.130826,  3302.881076,  3308.049253,  3309.319084,  3309.633581,  3309.711834;
 26288.872424, 25116.732923, 25247.211060, 25295.640045, 25307.418756, 25310.320376, 25311.040427;
 205876.611184,189409.780430,190251.946257,190686.715992,190791.471796,190817.134788,190823.485571;
1633835.840000,1423354.800000,1426265.110000,1430063.690000,1430971.940000,1431193.170000,1431247.760000;
];
log_kappa_iga_mat = log10(kappa_iga_mat);

figure(1)
hold on;
for pp = 1:length(p_vals)
    plot(loghh, log_kappa_full_mat(pp,:), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    text(loghh(end), log_kappa_full_mat(pp, end), ['p = ' num2str(p_vals(pp))], ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 20, 'FontWeight', 'bold');
end
xlabel('log10(h)');
ylabel('log10(\kappa)');

figure(2)
hold on;
for pp = 1:length(p_vals)
    plot(loghh, log_kappa_iga_mat(pp,:), 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    text(loghh(end), log_kappa_iga_mat(pp, end), ['p = ' num2str(p_vals(pp))], ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 20, 'FontWeight', 'bold');
end
xlabel('log10(h)');
ylabel('log10(\kappa)');