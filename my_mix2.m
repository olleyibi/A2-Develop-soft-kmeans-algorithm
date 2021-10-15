clc
clear all
close all


data = dlmread('old_faithful.dat','\t',26,0);
data = data(:,2:3);
x2 = linspace(min(data(:,2)),max(data(:,2)));
x1 = linspace(min(data(:,1)),max(data(:,1)));
[X1 X2] = meshgrid(x1, x2);


sigma1=[.2 0;0 46];
sigma2=[.4 0;0 93];

centroid = [randperm(floor(max(data(:,1))),2)' randperm(floor(max(data(:,2))),2)'];
mu1 = centroid(1, :);
mu2 = centroid(2, :);

% Mixing proportion
prop = [0.4 0.6];
mix = @(x) prop(1)*mvnpdf(x,mu1,sigma1) + prop(2)*mvnpdf(x,mu2,sigma2);



h1 = scatter(data(:, 1), data(:, 2),'cyan');axis square; box on;
h2 =  contour(X1, X2, reshape(mix([X1(:) X2(:)]), 100, 100),5);
h3 = plot(centroid(:, 1), centroid(:, 2), 'ksq','markersize',10);
xlabel('eruption duration (min)'); ylabel('time to next eruption (min)')
hold on


for i=1:20
    pause(1)
    clf
    c_respp = centroid_c_respponsibility(data, centroid);
    mu1 = centroid(1, :);
    mu2 = centroid(2, :);
    mix = @(x) prop(1)*mvnpdf(x,mu1,sigma1) + prop(2)*mvnpdf(x,mu2,sigma2);
    figure(1)
    scatter(data(:, 1), data(:, 2), 'cyan');axis square; box on;
    hold on; contour(X1, X2, reshape(mix([X1(:) X2(:)]), 100, 100),5);
    hold on; plot(centroid(:, 1), centroid(:, 2), 'ksq','markersize',10);
    %uistack(h1, 'bottom');
    centroid = data_assign(data, centroid, c_respp);
end
centroid



function c_resp=centroid_c_respponsibility(data, centroid)
[clust, ~] = size(centroid);
[num_points, ~] = size(data);
b = -1; % The beta value
c_resp = zeros(clust, num_points);

for i=1:clust
    c_resp(i, :) = exp(b*vecnorm((data-centroid(i, :))'));
end
c_resp = c_resp./sum(c_resp); % Column-wise normalisation
end


function centroidter=data_assign(data, centroid, c_respponsibility)
[clust, dim] = size(centroid);
[num_points, ~] = size(data);
centroidter = zeros(clust, dim);
total_c_resp = sum(c_respponsibility');
for i=1:clust
    centroid_new = zeros(1, dim);
    for j=1:num_points
        centroid_new = centroid_new + data(j, :).*c_respponsibility(i, j);
    end
    centroidter(i, :) = centroid_new./total_c_resp(1, i);
end
end