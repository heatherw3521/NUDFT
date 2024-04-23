% imaging examples
%%
clear all
close all
load('image_test_noise.mat')  %The "image" signal
%load('decay_test_noise.mat') %smooth signal assoc w/decaying coeffs.

%%
% 
% outputs: 
%
% each cell corresponds to a noise level. See "noises" var. 
%
% within each cell: another cell array: 2 by 5
%    two rows. Row 1 = jittered grid nodes. Row 2 = idd rand + gaps nodes. 
%    five cols: 

% col 1 = backslash 
% col 2 = CG normal eqns
% col 3 = CG adjoint normal eqns
% col 4 = ADI-based LS (INUDFT)
% col 5 = ADI-based Toeplitz on normal eqns

% within the inner cells there are structs: 
% within each struct: 
% yp = reconstructed signal on equi grid: yp = V*cfs, where cfs were computed by
% solving V*(unknown) = rhs. 
% y = values at sample locations
% xp = locations on [0, 2*pi] where we reconstruct
% x = locations where signal was sampled for reconstruction
% err = rel err on coeff recovery 
% resid = rel resid 
% sig = true signal on grid pts
% cfs = true coeffs
% pts = name of point set used for sampling
% solver = name of solver used to find cfs. 

%%
% Noise-free case:

out_noisefree = outputs{1}; 
%%
% bad node case: 
direct = out_noisefree{2,1}; 
CGnor = out_noisefree{2,2}; 
CGadj = out_noisefree{2,3}; 
INUDFT = out_noisefree{2,4};
FDToep = out_noisefree{2,5}; 

x = direct.x/2/pi; % where the signal was sampled, rescaled to [0, 1).
xp = (direct.xp)/2/pi; % equi grid rescaled to [0, 1). 
sig = direct.sig; %true signal 
 
%%
% plot the signal and a the set of sampled points. 
colorOrder = get(gca, 'ColorOrder');
clf
plot(xp, real(sig), '-k', 'linewidth', 1); 
hold on
plot(x, real(direct.y),'.', 'markersize', 6, 'color', colorOrder(1,:)); 
set(gca, 'fontsize', 14)
set(gcf,'color','w')
%export_fig image_true_w_samples.pdf
%export_fig image_true_w_no_samples.pdf

%%
% plot reconstruction errors on semilog scale: 
clf
semilogy(xp, abs(sig -CGnor.yp), 'color', colorOrder(2,:))
hold on
semilogy(xp, abs(sig - CGadj.yp), 'color', colorOrder(3,:))
semilogy(xp, abs(sig - INUDFT.yp), 'color', colorOrder(4,:))
%semilogy(xp, abs(sig - FDToep.yp), 'color', colorOrder(3,:))
legend('CGnor', 'CGadj', 'INDUFT', 'FDToep')
set(gca, 'fontsize', 14)
set(gcf,'color','w')
xlim([.4, 1]); %for SHORT version
export_fig errors_image_nonoise_gaps_error_short.pdf

%%
% a zoomed-in plot of the reconstructed signal vs. true signal. 
clf
plot(xp, real(sig), '-k', 'linewidth', 1.2); %true signal
hold on
plot(x, real(direct.y),'.', 'markersize', 12, 'color', colorOrder(1,:)); % sampled values
plot(xp, real(CGnor.yp), 'color', colorOrder(2,:), 'linewidth', 1.2 ) %CGnorsig reconstructed
plot(xp, real(CGadj.yp), 'color', colorOrder(3,:), 'linewidth', 1.2 )  
plot(xp, real(FDToep.yp), 'color', colorOrder(4,:), 'linewidth', 1.2 )  
%plot(xp, real(FDToep.yp), 'color', colorOrder(3,:) ) 
xl = [.4671, .4835]; 
yl = [-.6476, 1.5];
xlim(xl); 
ylim(yl); 
set(gca, 'fontsize', 14)
set(gcf,'color','w')
legend('', '', 'CGnor', 'CGadj', 'ADI-based LS solver')
%export_fig imagezoomin_nonoise_gaps.pdf

%%
% same plot as above but not zoomed in: 
clf
plot(x, real(direct.y),'.', 'markersize', 8, 'color', colorOrder(1,:)); % sampled values
hold on
plot(xp, real(sig), '--k', 'linewidth', 1); %true signal
plot(xp, real(CGnor.yp), 'color', colorOrder(2,:), 'linewidth', 1 ) %CGnorsig reconstructed
plot(xp, real(CGadj.yp), 'color', colorOrder(3,:), 'linewidth', 1)  
plot(xp, real(INUDFT.yp), 'color', colorOrder(4,:), 'linewidth', 1 )  
%plot(xp, real(FDToep.yp), 'color', colorOrder(3,:) ) 
xl = [.4671, .4835]; 
yl = [-.6476, 1.5];
%xlim(xl); 
%ylim(yl); 
set(gca, 'fontsize', 14)
set(gcf,'color','w')
%add a little box around the zoomed in part: 
gry = [1, 1, 1]/2; 
%p = patch([xl(1), xl(1), xl(2), xl(2), xl(1)],[yl(1), yl(2), yl(2), yl(1), yl(1)], gry, 'FaceAlpha',0.2);
p = rectangle('Position', [xl(1), yl(1), xl(2)-xl(1), yl(2)-yl(1)], 'linewidth', 2, 'Edgecolor', gry);
xlim([.4, 1]); %for SHORT version
legend('data', 'true signal', 'CGnor', 'CGadj', 'ADI-based LS solver', 'location', 'southeast')
ylim([-1.5, 1.5])
set(gca, 'XTick', [])
%export_fig imagefull_nonoise_gaps_rect_short.pdf

%%
% a plot of noisy data and the signals recovered:
out_noisemed = outputs{4};
% good node case: 
directn = out_noisemed{1,1}; 
CGnorn = out_noisemed{1,2}; 
CGadjn = out_noisemed{1,3}; 
INUDFTn = out_noisemed{1,4};
FDToepn = out_noisemed{1,5}; 

xp = (directn.xp)/2/pi; %equi grid rescaled;  
sig = directn.sig; % true signal; 
x = (directn.x)/2/pi; %sample locations rescaled
data = directn.rhs; % data (sb noisy)
%%
% plot noisy data with the recovered signals no zoom: 
figure(2)
plot(x, real(data), '.k', 'markersize', 5); 
hold on
plot(xp, real(sig), '--', 'color', [1, 1, 1]/2)
plot(xp, real(CGnorn.yp), 'color', colorOrder(2,:), 'linewidth', 1 ) %CGnorsig reconstructed
plot(xp, real(CGadjn.yp), 'color', colorOrder(3,:), 'linewidth', 1)  
plot(xp, real(INUDFTn.yp), 'color', colorOrder(4,:), 'linewidth', 1 ) 
legend('data', 'true signal', 'CGnor', 'CGadj', 'ADI-based LS solver', 'location', 'southeast')
set(gca, 'fontsize', 14)
set(gcf,'color','w')
xlim([.4, 1]); %for SHORT version
ylim([-1.5, 1.5])
xl = [.7030, .7183];
yl = [.8972, 1.0594];
p = rectangle('Position', [xl(1), yl(1), xl(2)-xl(1), yl(2)-yl(1)], 'linewidth', 2, 'Edgecolor', gry);
%export_fig imagefull_noise_equi.pdf

%%
% plot noisy data with the recovered signals, zoomed in: 
clf
plot(x, real(data), '.', 'markersize', 12, 'color', colorOrder(1,:)); 
hold on
plot(xp, real(sig), '--k')
plot(xp, real(CGnorn.yp), 'color', colorOrder(2,:), 'linewidth', 1 ) %CGnorsig reconstructed
plot(xp, real(CGadjn.yp), 'color', colorOrder(3,:), 'linewidth', 1)  
plot(xp, real(INUDFTn.yp), 'color', colorOrder(4,:), 'linewidth', 1 ) 
xl = [.7030, .7183];
yl = [.8972, 1.0594];
xlim(xl); 
ylim(yl); 
legend('data', 'true signal', 'CGnor', 'CGadj', 'ADI-based LS solver', 'location', 'southeast')
set(gca, 'fontsize', 14)
set(gcf,'color','w')
%export_fig imagezoomed_noise_equi.pdf
%p = rectangle('Position', [xl(1), yl(1), xl(2)-xl(1), yl(2)-yl(1)], 'linewidth', 2, 'Edgecolor', gry);




