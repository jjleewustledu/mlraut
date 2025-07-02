%% Debug Ward Construction
% Simple test to verify the implementation works

clear; close all; clc;

fprintf('Debug Ward Construction Test\n');
fprintf('===========================\n\n');

%% Step 1: Create and verify object
ward = WardConstructionFull(4, 2.0);  % Very small grid for debugging

fprintf('Object created successfully\n');
fprintf('Contour points structure:\n');
fprintf('  n_total = %d\n', ward.contour_points.n_total);
fprintf('  length(zeta) = %d\n', length(ward.contour_points.zeta));
fprintf('  length(weights) = %d\n', length(ward.contour_points.weights));

%% Step 2: Test single point computation
x_test = [0, 0, 0, 0];
twistor_test = @(Z) 1 / (Z(1) + 0.1);

fprintf('\nTesting single point computation...\n');
try
    A_test = ward.computeGaugeAtPoint(x_test, twistor_test, 'U1');
    fprintf('Success! A = [%.3e, %.3e, %.3e, %.3e]\n', ...
            A_test(1), A_test(2), A_test(3), A_test(4));
catch ME
    fprintf('Error: %s\n', ME.message);
    fprintf('Stack trace:\n');
    for i = 1:length(ME.stack)
        fprintf('  %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
end

%% Step 3: Test derivative computation directly
if exist('ward', 'var') && isvalid(ward)
    fprintf('\nTesting derivative computation...\n');
    
    % Create test twistor
    Z_test = [1; 0; 1; 0];
    
    % Test at different contour points
    test_indices = [1, ward.contour_points.n_total/2, ward.contour_points.n_total];
    
    for idx = test_indices
        if idx <= ward.contour_points.n_total && idx >= 1
            try
                dg = ward.computeDerivative(twistor_test, Z_test, round(idx), 1);
                fprintf('  Index %d: dg = [%.3e, %.3e, %.3e, %.3e]\n', ...
                        round(idx), dg(1), dg(2), dg(3), dg(4));
            catch ME
                fprintf('  Index %d: Error - %s\n', round(idx), ME.message);
            end
        end
    end
end

%% Step 4: Minimal field computation
fprintf('\nAttempting minimal field computation (2x2x2x2 grid)...\n');
ward_mini = WardConstructionFull(2, 1.0);

simple_twistor = @(Z) 1 / (Z(3) + 0.1);  % Simple pole

tic;
try
    [A_mu, success] = ward_mini.computeGaugeField(simple_twistor, 'U1');
    comp_time = toc;
    
    if success
        fprintf('Success! Computation time: %.2f seconds\n', comp_time);
        
        % Check field magnitude
        A1_max = max(abs(A_mu{1}(:)));
        A2_max = max(abs(A_mu{2}(:)));
        fprintf('Max |A_t| = %.3e\n', A1_max);
        fprintf('Max |A_x| = %.3e\n', A2_max);
    end
catch ME
    fprintf('Error in field computation: %s\n', ME.message);
end

fprintf('\nDebug test complete.\n');