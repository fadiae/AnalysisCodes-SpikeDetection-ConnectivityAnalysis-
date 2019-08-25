function [ pred_labels, pred_loss ] = loocv_svm( data, labels )
%LOOCV_SVM computes the leave-one-out (loo) cross-validation (cv) 
% classification loss and predicted labels for ``data'' and the ground 
% truth labels: ``labels''. The classification is performed using an SVM 
% classifier with a radial basis function (RBF) kernel.
%
% Input:    data        -   data matrix of size (samples X features)
%           labels      -   true data labels of size (features)
% Output:   pred_labels -   predicted labels after the loo-cv classification
%           pred_loss   -   the loss of the loo-cv classification

svm_model   = fitcsvm(data,labels,'Standardize','on','KernelFunction','RBF','KernelScale','auto');
loocv_model = crossval(svm_model,'Leaveout','on');

pred_loss   = kfoldLoss(loocv_model);
pred_labels = kfoldPredict(loocv_model);

end

