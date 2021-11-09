%% set default figure settings and get the color/line/marker order for plotting

[figSettings, color, rgbcmy, mList, lineList, dashList] = setFigureDefaults();

%%

T = readtable('mam4_plotting.csv', 'NumHeaderLines', 0);
plot_vars = T.Properties.VariableNames(5 : 14);

num = length(plot_vars);
t = T.time;

figure(1)
clf
hold on
for i = 1 : num
    
    name = plot_vars{i};
    field = T.(name);
    subplot(2, 5, i);
    plot(t, field)
    name = erase(name, 'CONC_');
    title(name);
    
end
