%% 3D Laffer curve
% 
% This script runs a sequence of counterfactuals varying both the fuel tax
% and the registration tax. 
% 
% 2021-07-14
%

close all; 
clear all; 
addpath('matlabinclude');
addpath('autotrade');

%% SETUP 

DOPLOTS = 1; % Show plots along the way 
DOSAVE  = 1; % Save the final graphs to disc 
colormap(summer); 

% load parameters
loaded = load('results/estimation/mle_converged.mat');
mp0 = loaded.mp_mle;
sol_mle = loaded.sol_mle; 
s = trmodel.index(mp0); 

% List of taxes 
num_fuel = 22; 
num_car = 22; 
fueltax_fraction_list = linspace(0., 3., num_fuel)';  
cartax_fraction_list = linspace(0., 1.5, num_car)'; % the relative change in the two car tax rates (baseline = 1.0). 

nIt = numel(fueltax_fraction_list) * numel(cartax_fraction_list); 

figOutPath = './results/laffer'; 

% we want (1.0,1.0) to be on grid for the plots 
assert(ismember(1., fueltax_fraction_list), 'Please ensure that the baseline value, 1.0, is in fueltax_fraction_list')
assert(ismember(1., cartax_fraction_list),  'Please ensure that the baseline value, 1.0, is in cartax_fraction_list')


%% STEP 0: baseline 

mp0.fixprices=0;
mp0.modeltype = 'structuralform';

p0=sol_mle.p; % price vector at mle estimates
sol0=equilibrium.solve(mp0, s, p0); % solve model in baseline
outcomes0 = stats.compute_outcomes(mp0, s, sol0); % comptue market outcomes 

%% STEP 1: Laffer iteration

h = graphs.waitbar('Computing Laffer curve'); 
it = 0; 

outcomeList = cell(num_fuel, num_car); 
mpList = cell(num_fuel, num_car); 
wait_times = nan(num_fuel*num_car, 1); 

for ifueltax = 1:numel(fueltax_fraction_list)
    for icartax=1:numel(cartax_fraction_list)
        it = it+1; 
        t_ = tic; 
        if it>1
        graphs.waitbar('Computing Laffer curve',it,nIt,h,wait_times); 
        end
        
        mp1=mp0;
        fueltax_fraction = fueltax_fraction_list(ifueltax); 
        cartax_fraction = cartax_fraction_list(icartax); 

        % set counterfactual taxes and compute implied prices
        mp1.tax_fuel=fueltax_fraction * mp0.tax_fuel;
        mp1.cartax_hi=mp0.cartax_hi*cartax_fraction;
        mp1.cartax_lo=mp0.cartax_lo*cartax_fraction;
        mp1=trmodel.update_mp(mp1); 

        sol1=equilibrium.solve(mp1, s, p0);
        p0=sol1.p; % update starting values for next iteration
        outcomes1 = stats.compute_outcomes(mp1, s, sol1);

        outcomeList{ifueltax, icartax} = outcomes1; 
        mpList{ifueltax, icartax} = mp1; 
        
        wait_times(it) = toc(t_); tic; 
    end
end

close(h) 


%% 3D Laffer curve
graphs.myfigure();  
xx = fueltax_fraction_list; 
yy = cartax_fraction_list; 

% loop and fill out tax revenues for each fuel price that we have solved
% above.
zz_tax = nan(num_fuel, num_car); 
zz_co2 = nan(num_fuel, num_car); 
zz_wel = nan(num_fuel, num_car); 
zz_wel_co2_low  = nan(num_fuel, num_car); 
zz_wel_co2_high = nan(num_fuel, num_car); 
zz_wel_co2_new = nan(num_fuel, num_car); 

for i=1:num_fuel
    for j=1:num_car
        zz_tax(i,j) = outcomeList{i,j}.total_revenue; 
        zz_co2(i,j) = outcomeList{i,j}.total_co2; 
        zz_wel(i,j) = outcomeList{i,j}.social_surplus_ex_co2; 
        zz_wel_co2_low(i,j)  = outcomeList{i,j}.social_surplus_total; 
        zz_wel_co2_high(i,j) = outcomeList{i,j}.social_surplus_total_highco2;
    end
end

% add dot at baseline (if the values are on grid)
i_x = find(fueltax_fraction_list == 1);
i_y = find(cartax_fraction_list == 1); 
if not(isempty(i_x) || isempty(i_y))
    h=scatter3(xx(i_x), yy(i_y), zz_tax(i_x, i_y)*1.02, 'r', 'filled'); 
    h.SizeData = 250; 
else 
    warning('Unable to find the point (1,1) on the grid!!! Cannot show the baseline points on the 3d graphs.'); 
    pause(2); 
    % it is in principle possible to interpolate the nearby points on the
    % zz-matrix to find the corresponding value off grid 
end



hold on; 

surf(xx,yy,zz_tax'); 
view([-140, 30]); 

xlabel('Fuel tax'); 
ylabel('Car tax'); 
zlabel('Tax revenue'); 
graphs.set_fig_layout_post(gcf, true, graphs.fontsize_3_wide);
hold off; 

if DOSAVE 
    name_ = sprintf('%s/laffer_3d.eps', figOutPath);
    saveas(gcf, name_, 'epsc'); 
    fprintf('Figure saved as <a href="%s">%s</a>\n', figOutPath, name_);
end

%% CO2 3d graph 
graphs.myfigure(); 
hold on; 

surf(xx,yy,zz_co2'); 
view([140, 30]); 

if not(isempty(i_x) || isempty(i_y))
    h=scatter3(xx(i_x), yy(i_y), zz_co2(i_x, i_y)*1.02, 'r', 'filled'); 
    h.SizeData = 250; 
end

xlabel('Fuel tax'); 
ylabel('Car tax'); 
zlabel('CO2'); 
graphs.set_fig_layout_post(gcf, true, graphs.fontsize_3_wide);
hold off; 

if DOSAVE 
    name_ = sprintf('%s/laffer_3d_co2.eps', figOutPath);
    saveas(gcf, name_, 'epsc'); 
    fprintf('Figure saved as <a href="%s">%s</a>\n', figOutPath, name_);
end


%% Welfare 3d graph 
graphs.myfigure(); 
hold on; 

surf(xx,yy,zz_wel'); 
view([140, 30]); 

if not(isempty(i_x) || isempty(i_y))
    h=scatter3(xx(i_x), yy(i_y), zz_wel(i_x, i_y)*1.02, 'r', 'filled'); 
    h.SizeData = 250; 
end

xlabel('Fuel tax'); 
ylabel('Car tax'); 
zlabel('Welfare (excl. CO2)'); 
graphs.set_fig_layout_post(gcf, true, graphs.fontsize_3_wide);
hold off; 

if DOSAVE 
    name_ = sprintf('%s/laffer_3d_welfare_ex_co2.eps', figOutPath);
    saveas(gcf, name_, 'epsc'); 
    fprintf('Figure saved as <a href="%s">%s</a>\n', figOutPath, name_);
end


%% Cross-sections of the 3d graph 

selected_outcomes = {'total_revenue', 'total_co2', 'social_surplus_total', 'consumer_surplus'}; 

% Line Style Order (LSO) 
myLSO = {'-o', '-*', '-x', '-d', '-v', '-+', '-<', '->', '-s', '-^'}; 
myLSO = repmat(myLSO, 1, 100); % just need plenty of line styles 

% cross-sections
for ivar=1:numel(selected_outcomes)

    var=selected_outcomes{ivar};
    figure(100+ivar); 
    hold on; 

    xx = cartax_fraction_list; 
    for i=1:num_fuel
        yy = nan(size(cartax_fraction_list)); 
        for j=1:numel(yy)
            yy(j) = outcomeList{i,j}.(var);
        end
        plot(xx, yy, myLSO{i}); 
    end

    hold off; 
    leg = legend(sprintfc('%5.0f%%',  fueltax_fraction_list * 100), ...
        'location','eastoutside'); 
    title(leg, 'Fuel tax');
    xlabel('Car tax (relative to baseline)'); 
    ylabel(stats.get_outcome_name(var)); 
    ax = gca; 
    ax.LineStyleOrder = {'-','--','-.',':'};
    graphs.set_fig_layout_post(gcf);

    if DOSAVE 
        name_ = sprintf('%s/laffer_3d_cross_%s.eps', figOutPath, var);
        saveas(gcf, name_, 'epsc'); 
        fprintf('Figure saved as <a href="%s">%s</a>\n', figOutPath, name_);
    end
end

%% ------------------------------------------------------------------------
%% Level curves 
%% ------------------------------------------------------------------------
% reshape data from array of structs to matrix form 

% variables that we want to explore 
vv = {'consumer_surplus', 'social_surplus_ex_co2', 'total_co2', 'total_revenue'}; 

% set up the grid
xx = fueltax_fraction_list; 
yy = cartax_fraction_list; 
[X,Y] = meshgrid(xx,yy); 

% fill out outcomes 
Z = nan(numel(xx),numel(yy),numel(vv)); 
for iv=1:numel(vv)
    v = vv{iv}; 
    it = 1; % resets at each iteration
    for ix=1:numel(xx)
        for iy=1:numel(yy)
            % verify that dimensions have not been misunderstood 
            assert(X(iy,ix) == xx(ix));
            assert(Y(iy,ix) == yy(iy));
            Z(iy,ix,iv) = outcomeList{it}.(v);            
            it = it+1;
        end
    end
end

% Plot separate level curve graphs 
for j=1:numel(vv)
    f=graphs.myfigure(); 
    v = vv{j}; 
    contour(X,Y,Z(:,:,j),'linewidth',2); 
    title(v, 'interpreter', 'none'); 
    graphs.set_fig_layout_post(f); 
    xlabel('Fuel tax'); ylabel('Registration tax (relative to baseline)'); 
end

%% Plot with all the level curves together  

close all

vv_subset = {'social_surplus_ex_co2', 'total_co2', 'total_revenue'}; 
cc = colororder; 
f=graphs.myfigure();
cons_ = {}; lines_ = {}; 
hold on; 
ic = 1; 
for j=1:numel(vv)
    v = vv{j}; 
    if ismember(v, vv_subset)
        level0 = outcomes0.(v); 
        [~, cons_{end+1}] = contour(X',Y',Z(:,:,j),[level0,level0],'color',cc(ic,:),'linewidth',2); 
        ic = ic+1; 
    end
end

% --- extract curves --- 

i1_ = 2; 
i2_ = 3; 
cm1_ = cons_{i1_}.ContourMatrix(:,2:end); % social surplus 
n_ = cons_{i2_}.ContourMatrix(2,1); % # of points in the contour
cm2_ = cons_{i2_}.ContourMatrix(:,2:n_+1); % co2

% --- add lower "triangle" area --- 

% find edges of the graphing area  
xlim = [min(fueltax_fraction_list), max(fueltax_fraction_list)]; 
ylim = [min(cartax_fraction_list), max(cartax_fraction_list)]; 

% % add lower point: duplicate for CO2 but edge for social surplus
% cm1_ = [[xlim(2);ylim(1)], cm1_];
% cm2_ = [cm2_(:,1), cm2_];
% 
% % add second point 
% cm1_ = [cm1_, cm1_(:,end)];
% cm2_ = [cm2_, [xlim(1);ylim(2)]];

cm1_ = [ cm1_, [xlim(2);ylim(1)]];
cm2_ = [ cm2_, cm2_(:,end)];

% add second point 
cm1_ = [cm1_];
cm2_ = [cm2_];

% area between curves below baseline (worse CO2 or welfare)
% cm1 = cm1_(:, cm1_(1,:) <= 1.0);
% cm2 = cm2_(:, cm2_(1,:) <= 1.0);
% xx_ = [cm1(1,:), fliplr(cm2(1,:))]; 
% yy_ = [cm1(2,:), fliplr(cm2(2,:))]; 
% fill(xx_,yy_,0.6 * [1,.1,.1]); alpha(0.2);

% area between curves above baseline: higher CO2 and/or welfare 
cm1 = [[1;1], cm1_(:, cm1_(1,:) >= 1.0)];
cm2 = [[1;1], cm2_(:, cm2_(1,:) >= 1.0)];
xx_ = [cm1(1,:), fliplr(cm2(1,:))]; 
yy_ = [cm1(2,:), fliplr(cm2(2,:))]; 

% --- do the actual plot --- 
mycolor = 0.6 * [.1,1,.1]; % green-ish hue 
fill(xx_,yy_,mycolor); 
alpha(0.2); 

% --- add "dots" --- 
offset = 1.03; 

% baseline dot 
plot(1.,1.,'or','markersize',10., 'markerfacecolor', 'r'); 
t_ = text(1.*.7, 1.*.95, 'Baseline'); 
t_.FontName = graphs.fontname; 
set(gcf, 'units', 'normalized', 'outerposition', [0 0 .55 .5]);

% tax maximized 
[i1,i2] = ind2sub(size(zz_tax), find(zz_tax(:) == max(zz_tax(:)))); 
plot(fueltax_fraction_list(i1), cartax_fraction_list(i2), 'or', 'markersize', 10., 'markerfacecolor', 'r'); 
t_ = text(fueltax_fraction_list(i1)*.9, cartax_fraction_list(i2)*1.36, ... 
    'Revenue max'); 
t_.FontName = graphs.fontname; 

% welfare maximized 
[i1,i2] = ind2sub(size(zz_wel), find(zz_wel_co2_low(:) == max(zz_wel_co2_low(:)))); 
plot(fueltax_fraction_list(i1), cartax_fraction_list(i2), 'or', 'markersize', 10., 'markerfacecolor', 'r'); 
t_ = text(fueltax_fraction_list(i1)*.68, cartax_fraction_list(i2)+.06, ...
    'Welf. max at 50$/ton'); 
t_.FontName = graphs.fontname; 

% welfare maximized 
[i1,i2] = ind2sub(size(zz_wel), find(zz_wel_co2_high(:) == max(zz_wel_co2_high(:)))); 
plot(fueltax_fraction_list(i1), cartax_fraction_list(i2), 'or', 'markersize', 10., 'markerfacecolor', 'r'); 
t_ = text(fueltax_fraction_list(i1)*.88, cartax_fraction_list(i2)+.06, ...
    'Welf. max at 250$/ton'); 
t_.FontName = graphs.fontname; 

% --- text --- 
legend('Social welfare (excl. CO2) constant', 'CO2 emissions constant', 'Tax revenue constant', ...
    'Reduced CO2, higher welfare and tax revenue', 'interpreter', 'none'); 
xlabel('Fuel tax'); ylabel('Registration tax');
hold off; 

graphs.set_fig_layout_post(f, false, 21); 


% --- save --- 
if DOSAVE 
    name_ = sprintf('%s/laffer_contours.eps', figOutPath);
    saveas(gcf, name_, 'epsc');
    fprintf('Figure saved as <a href="%s">%s</a>\n', figOutPath, name_);
end
