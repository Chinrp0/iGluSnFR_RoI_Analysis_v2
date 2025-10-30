%% Step 1: Run setup (verify everything is loaded)
setup_pipeline

%% Step 2: Set your data folder path
folder_path = 'D:\Data\GluSnFR\Ms\2025-10-25_CP_Ms-hipp_Snfr4f-NGR_Doc2b_1mMCa_resave\spont\GPU_SNR_Processed\5_raw_mean';

%% Step 3: Run the integrated batch analysis (WT vs R213W comparison)
results = integrated_batch_analysis_final(folder_path);

%% Step 4: Optional - Create detailed trace comparison plots
trace_fig = plot_condition_roi_traces(results);

% ```
% 
% ## What This Will Do:
% 
% 1. **Automatically classify files** by searching for "WT" or "R213W" in filenames
% 2. **Process all 17 files** using the corrected Schmitt trigger detector
% 3. **Generate 4 comparison plots:**
%    - Event frequency comparison (with 0.01 Hz bins)
%    - Event amplitude comparison
%    - Active ROI comparison
%    - File-level summary
% 4. **Print statistical comparison** (Mann-Whitney U tests for frequency and amplitude)
% 
% ## Expected Output:
% 
% You should see console output like:
% ```
% === Complete Integrated Batch Analysis: WT vs R213W ===
% Folder: D:\Data\...
%   Found 10 WT files, 7 R213W files (17 total)
%   Processing files...
%   Statistical difference: p < 0.001 ***