%% ========================================================================
%% LOAD_DATASET - Load a data folder + Recordings.csv into the workspace
%% ========================================================================
% Asks for a data folder (dialog opens at tracenorm_config data_root), reads
% every _mean.csv in it once together with the folder's Recordings.csv, and
% leaves the result in the workspace as `dataset`. run_pipeline uses
% `dataset` whenever it is present, so the CSVs are not re-read on every run.
%
% Usage:
%   load_dataset              % pick a folder in the dialog
%   run_pipeline              % analyses the loaded dataset
%   run_pipeline              % ...again, still without reloading
%
% To skip the dialog (or when running without a display), set `data_folder`
% before calling:
%   data_folder = 'E:\...\E_RoI_Raw-Mean'; load_dataset
%% ========================================================================

setup_pipeline

evalc('cfg_ = tracenorm_config();');
if ~exist('data_folder', 'var') || isempty(data_folder)
    data_folder = select_data_folder(cfg_);
end
dataset = load_dataset_folder(data_folder, cfg_);
clear cfg_ data_folder
