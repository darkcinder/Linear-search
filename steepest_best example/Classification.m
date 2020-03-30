% This file is a template to help you complete HW 2, problem 4.

global DATA
global LABELS

% Loads data as the variable DATA from folder CURRENT_DIR/dataset
% DATA: 100x2, each x_i is DATA(i, :)' 
load('dataset/DATA.mat');

% Loads labels as the variable LABELS from folder CURRENT_DIR/dataset
% LABELS: 100x1 vector of +1 and -1
load('dataset/LABELS.mat');

% Data is layed out horizontally: [x1, y1; x2, y2; ...]
[N_pts, N_dim] = size(DATA);

% Solve for unknown parameter to fit the distribution to data
%theta_opt = steepest_descent([1;1])
theta_opt = newton([1;1])

% Loop over each point, and label according to prob(label, x, theta*)
label_pred = zeros(size(LABELS));

% Classify the points
for i = 1:N_pts
    % If prob(label = +1| x, theta*) > 0.5, assign to class +1
    if round(prob_fun(1, DATA(i, :)', theta_opt)) == 1
        label_pred(i) = 1;
    else
        label_pred(i) = -1;
    end
end

% Plot to see if our distribution correctly classifies the points
hold on
plot(LABELS, '*')
plot(label_pred, '-')
title('Prediction vs Actual', 'fontsize', 16)
legend('Labels', 'Predicted Labels')
xlabel('Data Point', 'fontsize', 14)
ylabel('Class', 'fontsize', 14)
hold off
