classdef Bootstrap
    %% Implements block bootstrap, stratified subsampling, especially useful for FultzMulti.
    %
    %  https://www.mathworks.com/matlabcentral/answers/660618-how-to-perform-block-boostrap-resampling-with-matlab
    %  https://www.mathworks.com/help/finance/example-stratified-sampling.html
    %  https://www.mathworks.com/matlabcentral/fileexchange/52741-statistics-resampling?status=SUCCESS
    %
    % JJL:
    % The underlying time-series are obtained from fMRI, which is known to produce time-series that are non-Gaussian,
    % possibly related to EPI pulse sequences used.   Since the report of spurious statistics in fMRI made by Eklund,
    % Nichols, and Knutsson in PNAS 13(28) 7900 2016, there's strong preference for using permutation resampling with
    % fMRI time-series, so I do not wish to use theoretical distributions.   For N=4000 recordings of the f and g
    % time-series, how should I proceed to produce confidence levels that are interpretable and accurate?   I would like
    % to produce figures that show coherence and phase spectra with trend curves with respect to frequency.  I would
    % also like to have colored shadows indicating the confidence intervals for the trends.EditGiven the fMRI context
    % and the Eklund et al. findings about spurious statistics from parametric assumptions, you're absolutely right to
    % avoid theoretical distributions. Here's how to proceed with N=4000: Recommended Approach: 1. Block Bootstrap for
    % Temporal Dependencies
    %
    % Claude Sonnet 4:
    % Given the fMRI context and the Eklund et al. findings about spurious statistics from parametric assumptions,
    % you're absolutely right to avoid theoretical distributions. Here's how to proceed with N=4000: Recommended
    % Approach: 
    % 
    % 1. Block Bootstrap for Temporal Dependencies
    %
    % Use moving block bootstrap rather than simple bootstrap, since fMRI time series have strong temporal
    % autocorrelation Block length should be chosen based on your TR and the autocorrelation structure (typically 10-30
    % TRs) This preserves the temporal dependencies that are crucial in fMRI data
    %
    % 2. Stratified Subsampling for Computational Efficiency
    %
    % With N=4000, full bootstrap is computationally intensive Use stratified random subsampling: repeatedly sample
    % (without replacement) subsets of ~500-1000 trials Perform 1000-5000 such subsamples to build your confidence
    % distributions. This maintains the non-parametric approach while being computationally feasible
    %
    % 3. Frequency-Specific Confidence Intervals
    %
    % Compute coherence and phase at each frequency bin across all bootstrap/subsample iterations Use percentile-based
    % confidence intervals (e.g., 2.5th and 97.5th percentiles for 95% CI) Consider bias-corrected percentile intervals
    % if you observe systematic bias
    %
    % 4. Visualization Strategy
    %
    % Plot mean coherence/phase curves across frequency Add colored confidence bands using fill_between() or similar Use
    % semi-transparent shading (alpha ~0.3-0.5) Consider different colors for different confidence levels (e.g., darker
    % for 68%, lighter for 95%)
    %
    % 5. Multiple Comparisons Consideration
    %
    % With many frequency bins, consider cluster-based permutation testing (popular in neuroimaging) Or use false
    % discovery rate (FDR) correction across frequencies
    %
    % This approach respects the non-Gaussian nature of fMRI data while being computationally tractable and producing
    % publication-quality figures with interpretable confidence bands.
    %  
    %  Created 06-Jun-2025 16:53:07 by jjlee in repository /Users/jjlee/MATLAB-Drive/mlraut/src/+mlraut.
    %  Developed on Matlab 24.2.0.2923080 (R2024b) Update 6 for MACA64.  Copyright 2025 John J. Lee.
    
    methods
        function this = Bootstrap(varargin)
        end
    end
    
    %  Created with mlsystem.Newcl, inspired by Frank Gonzalez-Morphy's newfcn.
end
