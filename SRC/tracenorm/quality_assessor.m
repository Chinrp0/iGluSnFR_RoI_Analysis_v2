function quality_metrics = quality_assessor(raw_data, baseline_stats, dfof_stats, event_stats, config)
    % QUALITY_ASSESSOR - Assess ROI quality based on multiple criteria
    % Separate module for quality evaluation using baseline, dF/F, and event data
    %
    % Inputs:
    %   raw_data       - [frames x ROIs] original fluorescence data
    %   baseline_stats - Statistics from baseline_detector
    %   dfof_stats     - Statistics from dfof_calculator 
    %   event_stats    - Statistics from schmitt_event_detector
    %   config         - Configuration struct
    %
    % Outputs:
    %   quality_metrics - Comprehensive quality assessment for each ROI
    
    if nargin < 5
        config = tracenorm_config();
    end
    
    [numFrames, numROIs] = size(raw_data);
    
    if config.verbose
        fprintf('Assessing quality for %d ROIs...\n', numROIs);
    end
    
    %% === Extract Key Metrics ===
    % Baseline quality
    baseline_cv = baseline_stats.baseline_cv;
    valid_fraction = baseline_stats.valid_fraction;
    transport_flags = baseline_stats.transport_rois;
    
    % dF/F quality  
    snr = dfof_stats.snr;
    dynamic_range = dfof_stats.dynamic_range;
    fraction_valid_per_roi = dfof_stats.fraction_valid_per_roi;
    
    % Event characteristics
    events_per_roi = event_stats.events_per_roi;
    event_frames_per_roi = event_stats.event_frames_per_roi;
    
    %% === Define Quality Criteria ===
    % Use config parameters if available, otherwise use defaults
    if isfield(config, 'quality_criteria')
        criteria = config.quality_criteria;
    else
        % Default quality criteria
        criteria = struct();
        criteria.min_snr = 2.0;                    % Minimum signal-to-noise ratio
        criteria.max_baseline_cv = 0.2;            % Maximum baseline coefficient of variation (20%)
        criteria.min_valid_fraction = 0.8;         % Minimum fraction of valid baseline data (80%)
        criteria.min_dynamic_range = 0.01;         % Minimum dF/F dynamic range (1%)
        criteria.max_transport_slope = 0.1;        % Maximum slope for transport detection
        criteria.min_dfof_valid = 0.8;             % Minimum fraction of valid dF/F data
        criteria.min_events_for_active = 2;        % Minimum events to consider ROI "active"
        criteria.max_event_fraction = 0.5;         % Maximum fraction of frames that can be events
    end
    
    %% === Apply Quality Flags ===
    quality_flags = struct();
    
    % Baseline quality flags
    quality_flags.poor_baseline_stability = baseline_cv > criteria.max_baseline_cv;
    quality_flags.insufficient_baseline_data = valid_fraction < criteria.min_valid_fraction;
    quality_flags.transport_artifact = transport_flags;
    
    % Signal quality flags
    quality_flags.low_snr = snr < criteria.min_snr;
    quality_flags.low_dynamic_range = dynamic_range < criteria.min_dynamic_range;
    quality_flags.insufficient_dfof_data = fraction_valid_per_roi < criteria.min_dfof_valid;
    
    % Event-related flags
    quality_flags.no_events = events_per_roi == 0;
    quality_flags.few_events = events_per_roi < criteria.min_events_for_active & events_per_roi > 0;
    quality_flags.excessive_events = (event_frames_per_roi / numFrames) > criteria.max_event_fraction;
    
    %% === Overall Quality Classification ===
    % High quality ROIs: pass all major criteria
    high_quality = ~quality_flags.poor_baseline_stability & ...
                   ~quality_flags.insufficient_baseline_data & ...
                   ~quality_flags.transport_artifact & ...
                   ~quality_flags.low_snr & ...
                   ~quality_flags.low_dynamic_range & ...
                   ~quality_flags.insufficient_dfof_data;
    
    % Active ROIs: high quality + sufficient events
    active_rois = high_quality & events_per_roi >= criteria.min_events_for_active;
    
    % Transport ROIs: specifically flagged for transport
    transport_only = quality_flags.transport_artifact & quality_flags.no_events;
    
    % Poor quality: fail critical criteria
    poor_quality = quality_flags.low_snr | ...
                   quality_flags.poor_baseline_stability | ...
                   quality_flags.insufficient_baseline_data | ...
                   quality_flags.insufficient_dfof_data;
    
    %% === Quality Categories ===
    quality_categories = cell(1, numROIs);
    for roi = 1:numROIs
        if active_rois(roi)
            quality_categories{roi} = 'active';
        elseif transport_only(roi)
            quality_categories{roi} = 'transport_only';
        elseif high_quality(roi)
            quality_categories{roi} = 'high_quality_inactive';
        elseif poor_quality(roi)
            quality_categories{roi} = 'poor_quality';
        else
            quality_categories{roi} = 'moderate_quality';
        end
    end
    
    %% === Calculate Quality Scores ===
    % Composite quality score (0-1, higher is better)
    quality_scores = zeros(1, numROIs);
    
    for roi = 1:numROIs
        score = 0;
        total_weight = 0;
        
        % SNR contribution (weight: 3)
        if ~isnan(snr(roi)) && snr(roi) > 0
            score = score + 3 * min(1, snr(roi) / (2 * criteria.min_snr));
            total_weight = total_weight + 3;
        end
        
        % Baseline stability contribution (weight: 2)
        if ~isnan(baseline_cv(roi)) && baseline_cv(roi) > 0
            stability_score = max(0, 1 - baseline_cv(roi) / criteria.max_baseline_cv);
            score = score + 2 * stability_score;
            total_weight = total_weight + 2;
        end
        
        % Dynamic range contribution (weight: 2)
        if ~isnan(dynamic_range(roi)) && dynamic_range(roi) > 0
            range_score = min(1, dynamic_range(roi) / (2 * criteria.min_dynamic_range));
            score = score + 2 * range_score;
            total_weight = total_weight + 2;
        end
        
        % Event activity contribution (weight: 1)
        if events_per_roi(roi) > 0
            event_score = min(1, events_per_roi(roi) / (2 * criteria.min_events_for_active));
            score = score + 1 * event_score;
        end
        total_weight = total_weight + 1;
        
        % Data validity contribution (weight: 1)
        validity_score = min(valid_fraction(roi), fraction_valid_per_roi(roi));
        score = score + 1 * validity_score;
        total_weight = total_weight + 1;
        
        % Normalize by total weight
        if total_weight > 0
            quality_scores(roi) = score / total_weight;
        end
        
        % Penalize transport artifacts
        if transport_flags(roi)
            quality_scores(roi) = quality_scores(roi) * 0.5;
        end
    end
    
    %% === Compile Quality Metrics ===
    quality_metrics = struct();
    
    % Individual ROI metrics
    quality_metrics.flags = quality_flags;
    quality_metrics.categories = quality_categories;
    quality_metrics.scores = quality_scores;
    quality_metrics.criteria_used = criteria;
    
    % Summary statistics
    quality_metrics.summary = struct();
    quality_metrics.summary.total_rois = numROIs;
    quality_metrics.summary.high_quality_rois = sum(high_quality);
    quality_metrics.summary.active_rois = sum(active_rois);
    quality_metrics.summary.transport_only_rois = sum(transport_only);
    quality_metrics.summary.poor_quality_rois = sum(poor_quality);
    
    % Fractions
    quality_metrics.summary.fraction_high_quality = sum(high_quality) / numROIs;
    quality_metrics.summary.fraction_active = sum(active_rois) / numROIs;
    quality_metrics.summary.fraction_transport = sum(transport_only) / numROIs;
    quality_metrics.summary.fraction_poor = sum(poor_quality) / numROIs;
    
    % Quality score statistics
    quality_metrics.summary.mean_quality_score = mean(quality_scores, 'omitnan');
    quality_metrics.summary.median_quality_score = median(quality_scores, 'omitnan');
    quality_metrics.summary.quality_score_range = [min(quality_scores), max(quality_scores)];
    
    % Integration with existing pipeline (for backward compatibility)
    quality_metrics.fraction_good_rois = quality_metrics.summary.fraction_high_quality;
    quality_metrics.num_low_quality = quality_metrics.summary.poor_quality_rois;
    
    if config.verbose
        fprintf('  Quality assessment complete:\n');
        fprintf('    High quality: %d (%.1f%%)\n', sum(high_quality), 100*sum(high_quality)/numROIs);
        fprintf('    Active ROIs: %d (%.1f%%)\n', sum(active_rois), 100*sum(active_rois)/numROIs);
        fprintf('    Transport only: %d (%.1f%%)\n', sum(transport_only), 100*sum(transport_only)/numROIs);
        fprintf('    Poor quality: %d (%.1f%%)\n', sum(poor_quality), 100*sum(poor_quality)/numROIs);
    end
end