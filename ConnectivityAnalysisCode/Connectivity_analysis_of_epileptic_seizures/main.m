% MATLAB code implementation of the connectivity analysis, the
% classification of short spikes segments vs seizures and the complexity
% analysis (eigenvalue decomposition) of one experiment example in:
% "Layer and cell specific recruitment dynamics during epileptic seizures
% in-vivo” by F. Aeed, T. Shnitzer, R. Talmon and Y. Schiller.
% ***************************************************************@
%
% Author: Tal Shnitzer.
% Created:  6/17/19.
% ***************************************************************@

function main
%MAIN creates Figure 3e from the paper and Figures 3b and 3d corresponding
% to one experiment example (different than the ones in the paper).

load('data.mat');

ntr  = size(data,3);
klag = 2; % number of time-lags to use in the construction of the affinity kernels (per ROI).

%% Preprocessing
% Removing the mean and dividing by the standard deviation of each ROI to 
% avoid classification based on intensity:

norm_data = bsxfun(@minus,data,mean(data,2));
norm_data = bsxfun(@rdivide,data,std(norm_data,[],2));

% toim = [];
% for ii = 1:size(data,3)
%     toim = [toim; norm_data(:,:,ii)];
% end
% figure
% imagesc(toim);

%% Calculating the connectivity and correlations for each trial:
corr_all = cell(1,ntr);
conn_all = cell(1,ntr);

h = waitbar(0,'Calculating connectivity and correlation matrices for all trials');
for tr = 1:ntr
    waitbar(tr/ntr);
    
    %% Correlation (for refernce):
    corr_all{tr} = corrcoef(norm_data(:,:,tr).');
    
    % Adding a small value to the diagonal of the correlation matrix to 
    % avoid zero eigenvalues:
    corr_all{tr} = corr_all{tr} + 1e-3*eye(size(corr_all{tr}));
    
    %% Connectivity:
    conn_all{tr} = calcconn(norm_data(:,:,tr).', klag);

end
close(h)

%% Riemannian metric for connectivity and correlation analysis
% Analyzing the connectivity matrices and the correlation matrices using
% the Riemannian manifold of symmetric positive-definite (SPD) matrices:

% Thresholds for calculating the Riemannian mean of the matrices:
spd_eps     = 1e-6;
spd_maxiter = 200;

connfeat = calc_rieman_feat( cat(3,conn_all{:}), spd_eps, spd_maxiter );
corrfeat = calc_rieman_feat( cat(3,corr_all{:}), spd_eps, spd_maxiter );

% Dimensionality reduction using diffusion maps (dm):
ncomps  = 3;
ks_conn = 1;    % kernel scale (for the affinity kernel of dm)
ks_corr = 10;   % kernel scale (for the affinity kernel of dm) 
                % the correlation features require a larger kernel scale to avoid degeneracy
psi_conn = dmaps( connfeat, ncomps, ks_conn );
psi_corr = dmaps( corrfeat, ncomps, ks_corr ); 

%% Plotting Figure 3e
% Scatter plot of the connectivity features and the corresponding 
% correlation features figure for comparison:

cmap          = [1,0,0; 0,0,1];
segtype_color = cmap(labels+1,:);

figure('Name','Connectivity scatter plot (Figure 3c)')
scatter3(psi_conn(:,1),psi_conn(:,2),psi_conn(:,3),[],segtype_color,'fill');
hold on
h(1) = scatter3(inf,inf,inf,[],[1,0,0],'fill');
h(2) = scatter3(inf,inf,inf,[],[0,0,1],'fill');
hold off
legend(h,'Prolonged seizure onset','Short spikes onset','Location','northeast');

figure('Name','Correlation scatter plot')
scatter3(psi_corr(:,1),psi_corr(:,2),psi_corr(:,3),[],segtype_color,'fill');
hold on
h(1) = scatter3(inf,inf,inf,[],[1,0,0],'fill');
h(2) = scatter3(inf,inf,inf,[],[0,0,1],'fill');
hold off
legend(h,'Prolonged seizure onset','Short spikes onset','Location','northeast');

%% Leave-one-out cross-validation SVM classification based on the diffusion maps features

[ ~, pred_loss_conn ] = loocv_svm( psi_conn, labels );
[ ~, pred_loss_corr ] = loocv_svm( psi_corr, labels );

disp('------------------------------------------------------------------------------------------------------')
disp('Classification accuracy based on a leave-one-out cross-validation with an SVM classifier (RBF kernel):')
disp(['Connectivity features: ',num2str(100*(1-pred_loss_conn)),'%']);
disp(['Correlation  features: ',num2str(100*(1-pred_loss_corr)),'%']);
disp('------------------------------------------------------------------------------------------------------')

%% Plotting the corresponding Figure 3b (not on the same data as in the paper)
% The Riemannian mean of the connectivity matrices of short spikes onset vs
% prolonged seizure onset:

seiz_mats = cat(3,conn_all{labels==0});
shrt_mats = cat(3,conn_all{labels==1});
seiz_mean = calc_rieman_mean( seiz_mats, spd_eps, spd_maxiter );
shrt_mean = calc_rieman_mean( shrt_mats, spd_eps, spd_maxiter );

% Removing the diagonal values which are significantly higher than the rest
% and replacing them with the mean value, to obtain a clearer presentation 
% of the matrices:
seiz_mean_d = seiz_mean - diag(diag(seiz_mean)) + mean(seiz_mean(:))*eye(size(seiz_mean));
shrt_mean_d = shrt_mean - diag(diag(shrt_mean)) + mean(shrt_mean(:))*eye(size(shrt_mean));

[X,Y] = meshgrid(0:size(seiz_mean_d)-1,0:size(seiz_mean_d)-1);
figure('Name','Riemannian mean of connectivity matrices')
ax(1) = subplot(1,2,1); mesh(X,Y,seiz_mean_d); title('Prolonged seizure onset'); zlabel('Average connectivity matrix'); xlabel('Cells'); ylabel('Cells');
colormap(jet); zlim([0.5,0.8]); cax = caxis;
hold on
[~,hcl] = contour(X,Y,seiz_mean_d);
hcl.Fill = 'on';
hcl.ContourZLevel = 0.5;
hold off

ax(2) = subplot(1,2,2); mesh(X,Y,shrt_mean_d); title('Short spikes onset'); zlabel('Average connectivity matrix'); xlabel('Cells'); ylabel('Cells');
colormap(jet); zlim([0.5,0.8]); caxis(cax);
hold on
[~,hcs] = contour(X,Y,shrt_mean_d);
drawnow;
hcs.Fill = 'on';
hcs.ContourZLevel = 0.5;
hold off

Link = linkprop(ax,'CameraPosition');
setappdata(gcf, 'StoreTheLink', Link);
set(gcf,'position',[680,560,695,420]);

%% Plotting the corresponding Figure 3d (for only one experiment)

% Normalizing the average connectivity matrix to be row-stochastic:
seiz_mean_norm = diag(1./sum(seiz_mean)) * seiz_mean;
shrt_mean_norm = diag(1./sum(shrt_mean)) * shrt_mean;

seiz_eigs = eig(seiz_mean_norm);
shrt_eigs = eig(shrt_mean_norm);

erng = 2:20;
figure
plot(erng,seiz_eigs(erng),'r',erng,shrt_eigs(erng),'b');
xlabel('$$\ell$$','Interpreter','Latex','FontSize',16);
ylabel('Average $$\lambda_\ell$$','Interpreter','Latex','FontSize',16);
xlim([erng(1)-1,erng(end)+1]);
legend({'Prolonged seizure onset','Short spikes onset'},'FontSize',12);

end