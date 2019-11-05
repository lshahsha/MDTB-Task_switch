function [ varargout ] = sc1_sc2_mdtb_task_switch( what, varargin )
%[VARARGOUT] = SC1_SC2(WHAT, VARARGIN) This functions does every analysis
%that can be done on the mdtb dataset.
%   Detailed explanation goes here

numDummys = 3;   % number of dummy scans per run
numTRs    = 601; % number of scans per run

%%% setting path for the working directories
baseDir = '/Users/ladan/Documents/Project-Cerebellum/Cerebellum_Data';
% baseDir = '/home/ladan/Documents/Data/Cerebellum-MDTB';

%%% setting directory names
suitDir      = 'suit';                  %% directory where the anatomicals used in creating flatmaps are stored.
regDir       = 'RegionOfInterest';      %% The ROI directory 

suitToolDir  = '/Users/ladan/Documents/MATLAB/suit';

runLst  = 1:16;    %% run numbers to use in GLM

% cd(baseDir)

% Hard-coding some variables
%%% subjects
subj_name = {'s01','s02','s03','s04','s05','s06','s07','s08','s09','s10','s11',...
    's12','s13','s14','s15','s16','s17','s18','s19','s20','s21','s22','s23','s24',...
    's25','s26','s27','s28','s29','s30','s31'};
returnSubjs = [2, 3, 4, 6, 8, 9, 10, 12, 14, 15, 17, 18, 19, 20, 21, 22, 24, 25, 26, 27, 28, 29, 30, 31]; % "good" subjects


hemI    = {'L', 'R'}; % left and right hemisphere
hemName = {'CortexLeft', 'CortexRight'};

warning('off')
switch what
    case 'GLM:mdtb:contrast'
        %%% Calculating contrast images.
        % 'SPM_light' is created in this step (xVi is removed as it slows
        % down code for FAST GLM).
        % This case is written so that it works with both GLM 7 and GLM 8.
        % Reminder: GLM 7 was written with each condition as a separate
        % regressor and a regressor for each of the instructions. GLM 8 was
        % written with each task modeled as a 30 sec block and instructions
        % modeled as a separate regressor.
        % Example1: sc1_sc2_mdtb_task_switch('GLM:mdtb:contrast', 'sn', [17, 18], 'glm', 8, 'which', 'task')
        % Example2: sc1_sc2_mdtb_task_switch('GLM:mdtb:contrast', 'sn', [3], 'glm', 72, 'which', 'cond')
        
        sn             = returnSubjs;        %% list of subjects
        glm            = 8;              %% The glm number :)
        experiment_num = 1;
        con_vs         = 'average_1'; %% set it to 'rest' or 'average' (depending on the contrast you want)
        which          = 'task';      %% it can be set to either cond or task. set it to 'task for GLM_8 and 'cond' for GLM_7
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num', 'con_vs', 'which'})
        
        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'
        
        %%% setting directory paths I need
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        
        for s = sn
            fprintf('******************** calculating contrasts for %s ********************\n', subj_name{s});
            load(fullfile(glmDir, subj_name{s}, 'SPM.mat'))
            
            SPM  = rmfield(SPM,'xCon');
            T    = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            
            % t contrast for tasks
            ucondition = unique(T.(which));
            idx = 1;
            for tt = 1:length(ucondition) % 0 is "instruct" regressor
                switch con_vs
                    case 'rest_task' % contrast against rest
                        con                                  = zeros(1,size(SPM.xX.X,2));
                        con(:,logical(T.(which) == ucondition(tt))) = 1;
                        tmp = zeros(1, size(SPM.xX.X, 2));
                        tmp(:, strcmp(T.TN, 'rest')) = 1;
                        sum_rest = sum(tmp);
                        tmp = tmp./sum_rest;
                        con                                  = con/abs(sum(con));
                        con(:, tmp~=0) = -1/sum_rest;
                    case 'average_1' % contrast against the average of the other tasks including the instructions
                        % New: Eva, Oct 2nd
                        con                     = zeros(1,size(SPM.xX.X, 2));
                        con(1,logical(T.(which) == ucondition(tt))) = 1./sum(logical(T.(which) == ucondition(tt)));
                        con(1,logical(T.(which) ~= ucondition(tt))) = -1./sum(logical(T.(which) ~= ucondition(tt)));
                    case 'average_2' % contrast against the average of the other tasks not including the instructions
                        con        = zeros(1,size(SPM.xX.X, 2));
                        % TO TRY below - no instructions as contrasts
                        con(1,logical(T.(which) == ucondition(tt))) = 1./sum(logical(T.(which) == ucondition(tt)));
                        con(1,logical((T.(which) ~= ucondition(tt)) .* (T.inst == 0))) = -1./sum(logical((T.(which) ~= ucondition(tt)) .* (T.inst == 0)));
                    case 'rest'
                        con                     = zeros(1,size(SPM.xX.X,2));
                        con(:,logical(T.(which) == ucondition(tt))) = 1;
                        con                     = con/abs(sum(con));
                end
                
                name = sprintf('%s-%s',char(unique(T.TN(T.(which) == ucondition(tt)))), con_vs);
                
                SPM.xCon(idx) = spm_FcUtil('Set',name, 'T', 'c',con',SPM.xX.xKXs);
                idx=idx+1;
            end % tt (conditions)
            SPM = spm_contrasts(SPM,1:length(SPM.xCon));
            save('SPM.mat', 'SPM','-v7.3');
            SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
            save(fullfile(glmDir, subj_name{s}, 'SPM_light.mat'), 'SPM');

            % rename contrast images and spmT images
            conName = {'con','spmT'};
            for i = 1:length(SPM.xCon)
                for n = 1:numel(conName)
                    oldName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%2.4d.nii',conName{n},i));
                    newName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                    movefile(oldName{i}, newName{i});
                end % conditions (n, conName: con and spmT)
            end % i (contrasts)
        end % sn 
    case 'GLM:mdtb:contrast_task'
        %%% Calculating contrast images for tasks.
        % 'SPM_light' is created in this step (xVi is removed as it slows
        % down code for FAST GLM). The contrast for the task is created as
        % the average of the beta values for the conditions of a task
        % Example: sc1_sc2_mdtb_task_switch('GLM:mdtb:contrast_task', 'sn', [3])
        
        sn             = returnSubjs;  %% list of subjects
        glm            = 7;            %% The glm number :)
        experiment_num = 1;
        con_vs         = 'average_1';  %% set it to 'rest' or 'average_1' or 'average_2' (depending on the contrast you want)
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num', 'con_vs'})
        
        % gt the task info
        C   = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc  = getrow(C,C.StudyNum == experiment_num);
        
        conNames = unique(Cc.condNames);
        
        % in 'sc1_sc2_taskConds.txt' file, instruct is not coded as a
        % task/condition name. So I will have to add that to the list of
        % names
        conNames = ['Instruct'; conNames];
        
        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'
        
        %%% setting directory paths I need
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        
        for s = sn
            fprintf('******************** calculating contrasts for %s ********************\n', subj_name{s});
            load(fullfile(glmDir, subj_name{s}, 'SPM.mat'))
            
            SPM  = rmfield(SPM,'xCon');
            T    = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            
            % t contrast for tasks
            utask = unique(T.task);
            
            idx = 1;
            for tt = 1:length(utask) % 0 is "instruct" regressor
                switch con_vs
                    case 'rest' % contrast against rest
                        con                  = zeros(1,size(SPM.xX.X,2));
                        con(:,logical(T.task == utask(tt))) = 1;
                        con                  = con/abs(sum(con));
                    case 'average' % contrast against the average of all the tasks or all the other tasks???
                        con                                = zeros(1,size(SPM.xX.X,2));
                        con(:,T_rd.task == task(tt))       = 1;
                        con                                = con/abs(sum(con));
                        con                                = bsxfun(@minus, con, 1/length(T.cond));
                    case 'average_1' % contrast vs the average of all the tasks
                        con        = zeros(1,size(SPM.xX.X, 2));
                        con(1,logical(T.task == utask(tt))) = 1./sum(logical(T.task == utask(tt)));
                        con(1,logical(T.task ~= utask(tt))) = -1./sum(logical(T.task ~= utask(tt)));
                    case 'average_2' % contrast vs average of all the tasks except for instructions
                        con        = zeros(1,size(SPM.xX.X, 2));
                        % TO TRY below - no instructions as contrasts
                        con(1,logical(T.task == utask(tt))) = 1./sum(logical(T.task == utask(tt)));
                        con(1,logical(T.task ~= utask(tt) .* (T.inst == 0))) = -1./sum(logical((T.task ~= utask(tt)) .* (T.inst == 0)));
                end
                % fix the name!!!!!!!!!!!!!!!!
                name = sprintf('%s-%s_taskCon',char(unique(conNames(Cc.taskNum == utask(tt)))), con_vs);
                
                SPM.xCon(idx) = spm_FcUtil('Set',name, 'T', 'c',con',SPM.xX.xKXs);
                idx=idx+1;
            end % tt (tasks)
            SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
            save('SPM.mat', 'SPM','-v7.3');
            SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
            save(fullfile(glmDir, subj_name{s}, 'SPM_light.mat'), 'SPM');

            % rename contrast images and spmT images
            conName = {'con','spmT'};
            for i = 1:length(SPM.xCon)
                for n = 1:numel(conName)
                    oldName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%2.4d.nii',conName{n},i));
                    newName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                    movefile(oldName{i}, newName{i});
                end % conditions (n, conName: con and spmT)
            end % i (contrasts)
        end % sn        
    case 'GLM:mdtb:contrast_transitions_id'
        %%% Calculating contrast images for all transitions. There will be
        %%% 256 (16*16) transitions for sc1
        % 'SPM_light' is created in this step (xVi is removed as it slows
        % down code for FAST GLM)
        % Example: sc1_sc2_mdtb_task_switch('GLM:mdtb:contrast_transitions_id', 'sn', [3])
        
        sn             = returnSubjs;   %% list of subjects
        glm            = 8;             %% The glm number :)
        experiment_num = 1;             %% set to 1 for sc1 and 2 for sc2
        con_vs         = 'average_1';   %% contrasts were calculated vs 'rest' or 'average'
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num'})
        
        experiment = sprintf('sc%d', experiment_num); %% experiment number is converted to 'sc1' or 'sc2'
        
        %%% setting directory paths I need
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        
%         load(fullfile(glmDir, 'all_trans_task_names.mat'));
        for s = sn
            fprintf('******************** calculating contrasts for %s ********************\n', subj_name{s});
            load(fullfile(glmDir, subj_name{s}, 'SPM.mat'))
            
            SPM  = rmfield(SPM,'xCon');
            T    = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            
            if isfield(T, 'deriv') % derivatives included
                T_rd = getrow(T, T.deriv == 0);
            elseif ~isfield(T, 'deriv') % derivatives not included
                T_rd = T;
            end % are derivatives included?
            
            % load in trans_info
            load(fullfile(glmDir, 'trans_info.mat'));
            
            idx = 1;
            for tt = 1:max(t.run)
                for it = (unique(t.instOrder_run(t.run==tt)))'
                    
                    con = zeros(1,size(SPM.xX.X,2));
                    con(:,T_rd.run == tt & T_rd.instOrder == it & T_rd.inst == 1) = 1;
                     if con==0
                         keyboard;
                     end
                    con  = con/abs(sum(con));
                    name = sprintf('transition_%d-%s', idx, con_vs);
                    
                    SPM.xCon(idx) = spm_FcUtil('Set',name, 'T', 'c',con',SPM.xX.xKXs);
                    idx=idx+1;
                end
            end % tt run

            SPM = spm_contrasts(SPM,[1:length(SPM.xCon)]);
            save('SPM.mat', 'SPM','-v7.3');
            SPM = rmfield(SPM,'xVi'); % 'xVi' take up a lot of space and slows down code!
            save(fullfile(glmDir, subj_name{s}, 'SPM_light.mat'), 'SPM');

            % rename contrast images and spmT images
            conName = {'con','spmT'};
            for i = 1:length(t.instOrder_all)
                for n = 1:numel(conName)
                    oldName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%2.4d.nii',conName{n},i));
                    newName{i} = fullfile(glmDir, subj_name{s}, sprintf('%s_%s.nii',conName{n},SPM.xCon(i).name));
                    movefile(oldName{i}, newName{i});
                end % conditions (n, conName: con and spmT)
            end % i (contrasts)
        end % sn
    case 'GLM:mdtb:transition_info'
        % gets the task transition info. Uses one subject to get the
        % transition info for all the subjects, assuming that task
        % transitions are the same across all the subjects.
        % Example: sc1_sc2_mdtb_task_switch('GLM:mdtb:transition_info')
        
        sn             = 3;  %% the subject you want to use to get the transition info
        experiment_num = 1;
        glm            = 8;
        
        vararginoptions(varargin, {'experiment_num', 'glm'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directory paths
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        
        T = load(fullfile(glmDir, subj_name{sn}, 'SPM_info'));
        
        if isfield(T, 'deriv') % derivatives were included in the glm
            t = getrow(T, T.deriv == 0 & T.inst == 1); % getting all the non-derivative regressors and regressors corresponding to instructions
            t = rmfield(t,{'SN','TN', 'deriv', 'task'});         % removing unnecessary fields
        elseif ~isfield(T, 'deriv') % derivatives were not included in the glm
            t = getrow(T, T.inst == 1);
            t = rmfield(t,{'SN','TN','task'});         % removing unnecessary fields
        end % checking if derivatives were included in the glm

        t.instOrder_all = (1:length(t.inst))';
        t.instOrder_run = t.instOrder;
        t               = rmfield(t, {'instOrder'});
        
        save(fullfile(glmDir, 'trans_info'), 't', '-v7.3');
        
    case 'SURF:mdtb:map_con'
        % projects individual contrast map volume files for the conditions
        % to the workbench surface.
        % Example: sc1_sc2_mdtb('SURF:mdtb:map_con', 'sn', 2, 'glm', 8, 'experiment_num', 1)
    
        sn             = returnSubjs; %% list of subjects
        atlas_res      = 32;          %% set it to 32 or 164
        experiment_num = 1;           %% enter 1 for sc1 and 2 for sc2
        glm            = 7;           %% glm number
        con_vs         = 'average_1'; %% set it to 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'atlas_res', 'experiment_num', 'glm', 'con_vs'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting glm and surfaceWB directory
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        wbDir  = 'surfaceWB';
        glmSurfDir = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm));
        dircheck(glmSurfDir);
        
        for s = sn
            fprintf('******************** start mapping contrasts to surface for %s ********************\n', subj_name{s});
            subjSurfDir = fullfile(baseDir, 'sc1', wbDir, 'data', subj_name{s});
            dircheck(fullfile(glmSurfDir, subj_name{s})); %% directory to save the contrast maps
            
            T  = load(fullfile(glmDir, subj_name{s}, 'SPM_info.mat'));
            conNames = unique(T.TN);
            for h = 1:2 % two hemispheres
                white   = fullfile(subjSurfDir,sprintf('%s.%s.white.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                pial    = fullfile(subjSurfDir,sprintf('%s.%s.pial.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                C1      = gifti(white);
                C2      = gifti(pial);
                for ic = 1:length(conNames)
                    
                    conMapName      = sprintf('con_%s-%s', conNames{ic}, con_vs);
                    images{1}       = fullfile(glmDir, subj_name{s}, sprintf('%s.nii', conMapName));
                    column_name{1}  = sprintf('con_%s-%s.nii', conNames{ic}, con_vs);
                    outfile         = fullfile(glmSurfDir, subj_name{s}, sprintf('%s.%s.con_%s-%s.func.gii', ...
                        subj_name{s}, hemI{h}, conNames{ic}, con_vs));
                    G               = surf_vol2surf(C1.vertices,C2.vertices,images,'column_names', column_name, ...
                        'anatomicalStruct',hemName{h});
                    save(G, outfile);
                    fprintf('******************** mapping to surface for %s hemi %s contrast %s done********************\n', subj_name{s}, ...
                        hemI{h}, conMapName);
                end % ic (condition/contrast)   
            end % hemi
        end % sn  
    case 'SURF:mdtb:map_con_task'
        % projects individual contrast map volume files for tasks to
        % WorkBench surface
        % Example: sc1_sc2_mdtb_task_switch('SURF:mdtb:map_con_task', 'sn', [3])
    
        sn         = returnSubjs; %% list of subjects
        atlas_res  = 32;          %% set it to 32 or 164
        experiment_num = 1;           %% enter 1 for sc1 and 2 for sc2
        glm        = 8;           %% glm number
        con_vs     = 'average_1'; %% set it to 'rest' or 'average_1' or 'average_2'
        
        vararginoptions(varargin,{'sn', 'atlas_res', 'experiment_num', 'glm', 'con_vs'});
        
        % gt the task info
        C   = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc  = getrow(C,C.StudyNum == experiment_num);
        
        conNames = unique(Cc.condNames);
        
        % in 'sc1_sc2_taskConds.txt' file, instruct is not coded as a
        % task/condition name. So I will have to add that to the list of
        % names
        conNames = ['Instruct'; conNames];

        experiment = sprintf('sc%d', experiment_num);
        
        % setting directory paths
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        wbDir  = 'surfaceWB';
        
        glmSurfDir = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm));
        dircheck(glmSurfDir);
        
        for s = sn
            fprintf('******************** start mapping contrasts to surface for %s ********************\n', subj_name{s});
            subjSurfDir = fullfile(baseDir, 'sc1', wbDir, 'data', subj_name{s});
            dircheck(fullfile(glmSurfDir, subj_name{s})); %% directory to save the contrast maps
            
            % get the task names
            conNames = unique(Cc.taskNames);
            
            for h = 1:2 % two hemispheres
                white   = fullfile(subjSurfDir,sprintf('%s.%s.white.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                pial    = fullfile(subjSurfDir,sprintf('%s.%s.pial.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                C1      = gifti(white);
                C2      = gifti(pial);
                for ic = 1:length(conNames)
%                     if ~ strcmp(conNames{ic}, 'rest')
                        conMapName      = sprintf('con_%s-%s_taskCon', conNames{ic}, con_vs);
                        images{1}       = fullfile(glmDir, subj_name{s}, sprintf('%s.nii', conMapName));
                        column_name{1}  = sprintf('con_%s-%s_taskCon.nii', conNames{ic}, con_vs);
                        
                        outfile    = fullfile(glmSurfDir, subj_name{s}, sprintf('%s.%s.con_%s-%s_taskCon.func.gii', subj_name{s}, hemI{h}, conNames{ic}, con_vs));
                        G          = surf_vol2surf(C1.vertices,C2.vertices,images,'column_names', column_name, ...
                            'anatomicalStruct',hemName{h});
                        save(G, outfile);
                        fprintf('******************** mapping to surface for %s hemi %s contrast %s done********************\n', subj_name{s}, ...
                            hemI{h}, conMapName);
%                     end
                    
                end % ic (condition/contrast)   
            end % hemi
        end % sn
    case 'SURF:mdtb:map_con_utransitions_names'
        % projects individual contrast map volume files for the unique
        % transitions with the task names to workbench surface
        % Example: sc1_sc2_mdtb_task_switch('SURF:mdtb:map_con_utransitions_names', 'sn', [3])
    
        sn         = returnSubjs;   %% list of subjects
        atlas_res  = 32;            %% set it to 32 or 164
        experiment_num = 1;             %% enter 1 for sc1 and 2 for sc2
        glm        = 8;             %% glm number
        con_vs     = 'rest';        %% set it to 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'atlas_res', 'experiment_num', 'glm', 'con_vs'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directory paths
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        wbDir  = 'surfaceWB';
        
        load(fullfile(glmDir, 'all_trans_task_names.mat'));
        glmSurfDir = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm));
        dircheck(glmSurfDir);
        
        for s = sn
            fprintf('******************** start mapping contrasts to surface for %s ********************\n', subj_name{s});
            subjSurfDir = fullfile(baseDir, experiment, wbDir, 'data', subj_name{s});
            dircheck(fullfile(glmSurfDir, subj_name{s})); %% directory to save the contrast maps
            
            for tt = 1:size(all_trans, 1)
                for h = 1:2 % two hemispheres
                    white   = fullfile(subjSurfDir,sprintf('%s.%s.white.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                    pial    = fullfile(subjSurfDir,sprintf('%s.%s.pial.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                    C1      = gifti(white);
                    C2      = gifti(pial);
                    conMapName      = sprintf('con_transition_%s_%s-%s', all_trans{tt, 1}, all_trans{tt, 2}, con_vs);
                    images{1}       = fullfile(glmDir, subj_name{s}, sprintf('%s.nii', conMapName));
                    column_name{1}  = sprintf('con_transition_%s_%s-%s.nii', all_trans{tt, 1}, all_trans{tt, 2}, con_vs);

                    outfile    = fullfile(glmSurfDir, subj_name{s}, sprintf('%s.%s.con_transition_%s_%s-%s.func.gii', subj_name{s}, hemI{h}, all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
                    G          = surf_vol2surf(C1.vertices,C2.vertices,images,'column_names', column_name, ...
                                    'anatomicalStruct',hemName{h});
                    save(G, outfile);
                    fprintf('******************** mapping to surface for %s hemi %s contrast %s done********************\n', subj_name{s}, ...
                        hemI{h}, conMapName);
                end % hemi
            end % tt (transitions)
        end % sn    
    case 'SURF:mdtb:map_con_transitions_id'
        % projects individual contrast map volume files for all the
        % transitions (256 for sc1) to workbench surface
        % Example: , 'sn', [3])
    
        sn             = returnSubjs;   %% list of subjects
        atlas_res      = 32;     %% set it to 32 or 164
        experiment_num = 1;      %% enter 1 for sc1 and 2 for sc2
        glm            = 8;      %% glm number
        con_vs         = 'average_2'; %% set it to 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'atlas_res', 'experiment_num', 'glm', 'con_vs'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directory paths
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        wbDir  = 'surfaceWB';
        
        load(fullfile(glmDir, 'trans_info.mat'));
        glmSurfDir = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm));
        dircheck(glmSurfDir);
        
        for s = sn
            fprintf('******************** start mapping contrasts to surface for %s ********************\n', subj_name{s});
            subjSurfDir = fullfile(baseDir, experiment, wbDir, 'data', subj_name{s});
            dircheck(fullfile(glmSurfDir, subj_name{s})); %% directory to save the contrast maps
            
            for tt = 1:length(t.instOrder_all)
                for h = 1:2 % two hemispheres
                    white   = fullfile(subjSurfDir,sprintf('%s.%s.white.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                    pial    = fullfile(subjSurfDir,sprintf('%s.%s.pial.%dk.surf.gii',subj_name{s},hemI{h}, atlas_res));
                    C1      = gifti(white);
                    C2      = gifti(pial);
                    conMapName      = sprintf('con_transition_%d-%s', t.instOrder_all(tt), con_vs);
                    images{1}       = fullfile(glmDir, subj_name{s}, sprintf('%s.nii', conMapName));
                    column_name{1}  = sprintf('con_transition_%d-%s.nii', t.instOrder_all(tt), con_vs);

                    outfile    = fullfile(glmSurfDir, subj_name{s}, sprintf('%s.%s.con_transition_%d-%s.func.gii', subj_name{s}, hemI{h}, t.instOrder_all(tt), con_vs));
                    G          = surf_vol2surf(C1.vertices,C2.vertices,images,'column_names', column_name, ...
                                    'anatomicalStruct',hemName{h});
                    save(G, outfile);
                    fprintf('******************** mapping to surface for %s hemi %s contrast %s done********************\n', subj_name{s}, ...
                        hemI{h}, conMapName);
                end % hemi
            end % tt (transitions)
        end % sn
    case 'SURF:mdtb:groupmap_con'
        % creates group average contrast maps for task contrasts
        % Example: sc1_sc2_mdtb_task_switch('SURF:mdtb:groupmap_con', 'sn', [3])
    
        sn         = returnSubjs;   %% list of subjects
        atlas_res  = 32;     %% set it to 32 or 164
        experiment_num = 1;      %% enter 1 for sc1 and 2 for sc2
        glm        = 8;      %% glm number
        replaceNaN = 1;      %% replacing NaNs
        con_vs     = 'average_1'; %% contrast was calculated against 'rest' or 'average'        
        smooth     = 1;      %% add smoothing
        kernel     = 1;      %% for smoothing
        which      = 'task'; %% 'task' for glm8 and 'cond' for glm7
        
        vararginoptions(varargin,{'sn', 'atlas_res', 'experiment_num', 'glm', 'replaceNaN', 'con_vs', 'smooth', 'kernel', 'which'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc       = getrow(C, C.StudyNum == experiment_num);
        switch which
            case 'task' % task for glm8
                conNames = unique(Cc.taskNames);
            case 'cond' % condition for glm7
                conNames = unique(Cc.condNames);
        end %% do you want the group maps for tasks or conditions
        
        % in 'sc1_sc2_taskConds.txt' file, instruct is not coded as a
        % task/condition name. So I will have to add that to the list of
        % names
        conNames = ['Instruct'; conNames];
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directory
        wbDir  = 'surfaceWB';        
        
        % go to the directory where fs_LR atlas is.
        groupSurfDir     = fullfile(baseDir, experiment, wbDir, 'data', sprintf('group%dk', atlas_res));
        groupSurfDir_glm = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), sprintf('group%dk', atlas_res));
        dircheck(groupSurfDir_glm);
        cd(groupSurfDir);
        
        for h = 1:2 % two hemispheres
            for cc = 1:length(conNames)
                for s = 1:length(sn)
                    %%% make the group metric file for each contrasts
%                     infilenames{s}   = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)},sprintf('%s.%s.con_%s-%s_taskCon.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs));
                    infilenames{s}   = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)},sprintf('%s.%s.con_%s-%s.func.gii', subj_name{sn(s)}, hemI{h}, conNames{cc}, con_vs));
                    
                    if smooth
                        surfFile    = fullfile(groupSurfDir,sprintf('fs_LR.32k.%s.inflated.surf.gii',hemI{h}));
                        surf_smooth(infilenames{s},'surf',surfFile,'kernel',kernel); % smooth outfilenames - it will prefix an 's'
                    end
%                     s_infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_%s-%s_taskCon.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs));
                    s_infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_%s-%s.func.gii', subj_name{sn(s)}, hemI{h}, conNames{cc}, con_vs));
                
                end % sn
%                 outfilenames    = fullfile(groupSurfDir_glm,sprintf('%s.con_%s-%s_taskCon.func.gii',hemI{h},taskNames{cc}, con_vs));
%                 summaryname     = fullfile(groupSurfDir_glm,sprintf('%s.group.con_%s-%s_taskCon.func.gii',hemI{h},taskNames{cc}, con_vs));
                
                outfilenames    = fullfile(groupSurfDir_glm,sprintf('%s.con_%s-%s.func.gii',hemI{h},conNames{cc}, con_vs));
                summaryname     = fullfile(groupSurfDir_glm,sprintf('%s.group.con_%s-%s.func.gii',hemI{h},conNames{cc}, con_vs));
                
                surf_groupGiftis(infilenames, 'outfilenames', {outfilenames}, 'groupsummary', summaryname, 'replaceNaNs', replaceNaN);
                if smooth % also save the smoothed versions
%                     s_outfilenames    = fullfile(groupSurfDir_glm,sprintf('s%s.con_%s-%s_taskCon.func.gii', hemI{h},taskNames{cc}, con_vs));
                    s_outfilenames    = fullfile(groupSurfDir_glm,sprintf('s%s.con_%s-%s.func.gii', hemI{h},conNames{cc}, con_vs));
%                     s_summaryname     = fullfile(groupSurfDir_glm,sprintf('s%s.group.con_%s-%s_taskCon.func.gii', hemI{h},taskNames{cc}, con_vs));
                    s_summaryname     = fullfile(groupSurfDir_glm,sprintf('s%s.group.con_%s-%s.func.gii', hemI{h},conNames{cc}, con_vs));
                    surf_groupGiftis(s_infilenames, 'outfilenames', {s_outfilenames}, 'groupsummary', s_summaryname, 'replaceNaNs', replaceNaN);
                end
                
                % delet the smoothed contrasts for each subject
                for s = 1:length(sn)
                    delete(fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_%s-%s_taskCon.func.gii', subj_name{sn(s)}, hemI{h}, conNames{cc}, con_vs)));
                end
                
                fprintf('******************** group average contrast for %s vs %s is created! ********************\n\n', conNames{cc}, con_vs);
            end % contrasts(ic)
        end % hemi(h)
    case 'SURF:mdtb:groupmap_con_task'
        % creates group average contrast maps for task contrasts
        % Example: sc1_sc2_mdtb_task_switch('SURF:mdtb:groupmap_con_task', 'sn', [3])
    
        sn         = returnSubjs;   %% list of subjects
        atlas_res  = 32;     %% set it to 32 or 164
        experiment_num = 1;      %% enter 1 for sc1 and 2 for sc2
        glm        = 8;      %% glm number
        replaceNaN = 1;      %% replacing NaNs
        con_vs     = 'rest'; %% contrast was calculated against 'rest' or 'average'        
        smooth     = 1;      %% add smoothing
        kernel     = 1;      %% for smoothing
        
        vararginoptions(varargin,{'sn', 'atlas_res', 'experiment_num', 'glm', 'replaceNaN', 'con_vs', 'smooth', 'kernel'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc       = getrow(C, C.StudyNum == experiment_num);
        taskNames = unique(Cc.taskNames);
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directory
        wbDir  = 'surfaceWB';
        
        % go to the directory where fs_LR atlas is.
        groupSurfDir     = fullfile(baseDir, experiment, wbDir, 'data', sprintf('group%dk', atlas_res));
        groupSurfDir_glm = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), sprintf('group%dk', atlas_res));
        dircheck(groupSurfDir_glm);
        cd(groupSurfDir);
        
        for h = 1:2 % two hemispheres
            for cc = 1:length(taskNames)
                for s = 1:length(sn)
                    %%% make the group metric file for each contrasts
%                     infilenames{s}   = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)},sprintf('%s.%s.con_%s-%s_taskCon.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs));
                    infilenames{s}   = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)},sprintf('%s.%s.con_%s-%s.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs));
                    
                    if smooth
                        surfFile    = fullfile(groupSurfDir,sprintf('fs_LR.32k.%s.inflated.surf.gii',hemI{h}));
                        surf_smooth(infilenames{s},'surf',surfFile,'kernel',kernel); % smooth outfilenames - it will prefix an 's'
                    end
%                     s_infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_%s-%s_taskCon.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs));
                    s_infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_%s-%s.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs));
                
                end % sn
%                 outfilenames    = fullfile(groupSurfDir_glm,sprintf('%s.con_%s-%s_taskCon.func.gii',hemI{h},taskNames{cc}, con_vs));
%                 summaryname     = fullfile(groupSurfDir_glm,sprintf('%s.group.con_%s-%s_taskCon.func.gii',hemI{h},taskNames{cc}, con_vs));
                
                outfilenames    = fullfile(groupSurfDir_glm,sprintf('%s.con_%s-%s.func.gii',hemI{h},taskNames{cc}, con_vs));
                summaryname     = fullfile(groupSurfDir_glm,sprintf('%s.group.con_%s-%s.func.gii',hemI{h},taskNames{cc}, con_vs));
                
                surf_groupGiftis(infilenames, 'outfilenames', {outfilenames}, 'groupsummary', summaryname, 'replaceNaNs', replaceNaN);
                if smooth % also save the smoothed versions
%                     s_outfilenames    = fullfile(groupSurfDir_glm,sprintf('s%s.con_%s-%s_taskCon.func.gii', hemI{h},taskNames{cc}, con_vs));
                    s_outfilenames    = fullfile(groupSurfDir_glm,sprintf('s%s.con_%s-%s.func.gii', hemI{h},taskNames{cc}, con_vs));
%                     s_summaryname     = fullfile(groupSurfDir_glm,sprintf('s%s.group.con_%s-%s_taskCon.func.gii', hemI{h},taskNames{cc}, con_vs));
                    s_summaryname     = fullfile(groupSurfDir_glm,sprintf('s%s.group.con_%s-%s.func.gii', hemI{h},taskNames{cc}, con_vs));
                    surf_groupGiftis(s_infilenames, 'outfilenames', {s_outfilenames}, 'groupsummary', s_summaryname, 'replaceNaNs', replaceNaN);
                end
                
                % delet the smoothed contrasts for each subject
                for s = 1:length(sn)
                    delete(fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_%s-%s_taskCon.func.gii', subj_name{sn(s)}, hemI{h}, taskNames{cc}, con_vs)));
                end
                
                fprintf('******************** group average contrast for %s vs %s is created! ********************\n\n', taskNames{cc}, con_vs);
            end % contrasts(ic)
        end % hemi(h)
    case 'SURF:mdtb:groupmap_con_utransitions_names'
        % creates group average contrast maps for unique transitions with
        % task names
        % Example: sc1_sc2_mdtb_task_switch('SURF:mdtb:groupmap_con_utransitions_names', 'sn', [3])
    
        sn         = returnSubjs;   %% list of subjects
        atlas_res  = 32;     %% set it to 32 or 164
        experiment_num = 1;      %% enter 1 for sc1 and 2 for sc2
        glm        = 7;      %% glm number
        replaceNaN = 1;      %% replacing NaNs
        smooth     = 1;      %% add smoothing
        kernel     = 1;      %% for smoothing
        con_vs     = 'rest'; %% contrast was calculated against 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'atlas_res', 'experiment_num', 'glm', 'replaceNaN', 'convs', 'smooth', 'kernel'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directories
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        wbDir  = 'surfaceWB';
        
        load(fullfile(glmDir, 'all_trans_task_names.mat'));

        % go to the directory where the group data is 
        groupSurfDir = fullfile(baseDir, experiment, wbDir, 'data', sprintf('group%dk', atlas_res));
        cd(groupSurfDir);
        
        for tt = 1:size(all_trans, 1)
            for h = 1:2 % two hemispheres
                for s = 1:length(sn)
                    %%% make the group metric file for each contrasts
                    infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)},sprintf('%s.%s.con_transition_%s_%s-%s.func.gii', subj_name{s}, hemI{h}, all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
                    if smooth
                        surfFile    = fullfile(groupSurfDir,sprintf('fs_LR.32k.%s.inflated.surf.gii',hemI{h}));
                        surf_smooth(infilenames{s},'surf',surfFile,'kernel',kernel); % smooth outfilenames - it will prefix an 's'
                    end
                    s_infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_transition_%s_%s-%s.func.gii', subj_name{s}, hemI{h}, all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
                end % sn
                outfilenames    = fullfile(groupSurfDir,sprintf('%s.con_transition_%s_%s-%s.func.gii', hemI{h}, all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
                summaryname     = fullfile(groupSurfDir,sprintf('%s.group.con_transition_%s_%s-%s.func.gii', hemI{h}, all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
                
                surf_groupGiftis(infilenames, 'outfilenames', {outfilenames}, 'groupsummary', summaryname, 'replaceNaNs', replaceNaN);
                if smooth % also save the smoothed versions
                    s_outfilenames    = fullfile(groupSurfDir,sprintf('s%s.con_transition_%s_%s-%s.func.gii', hemI{h}, all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
                    s_summaryname     = fullfile(groupSurfDir,sprintf('s%s.group.con_transition_%s_%s-%s.func.gii', hemI{h}, all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
                    surf_groupGiftis(s_infilenames, 'outfilenames', {s_outfilenames}, 'groupsummary', s_summaryname, 'replaceNaNs', replaceNaN);
                end
                
                % delet the smoothed contrasts for each subject
                for s = 1:length(sn)
                    delete(fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_transition_%s_%s-%s.func.gii', subj_name{s}, hemI{h}, all_trans{tt, 1}, all_trans{tt, 2}, con_vs)));
                end
                
                fprintf('******************** group average contrast transition from %s to %s vs %s is created! ********************\n\n', all_trans{tt, 1}, all_trans{tt, 2}, con_vs);
            end % hemi(h)
        end % tt (transitions)
    case 'SURF:mdtb:groupmap_con_transitions_id'
        % creates group average contrast maps for all the transitions
        % Example: sc1_sc2_mdtb_task_switch('SURF:mdtb:groupmap_con_transitions_id', 'sn', [3])
     
        sn             = returnSubjs;   %% list of subjects
        atlas_res      = 32;     %% set it to 32 or 164
        experiment_num = 1;      %% enter 1 for sc1 and 2 for sc2
        glm            = 8;      %% glm number
        replaceNaN     = 1;      %% replacing NaNs
        smooth         = 1;      %% add smoothing
        kernel         = 1;      %% for smoothing
        con_vs         = 'rest'; %% contrast was calculated against 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'atlas_res', 'experiment_num', 'glm', 'replaceNaN', 'con_vs', 'convs', 'smooth', 'kernel'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directories
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        wbDir  = 'surfaceWB';
        
        load(fullfile(glmDir, 'trans_info.mat'));

        % go to the directory where the group data is 
        groupSurfDir = fullfile(baseDir, experiment, wbDir, 'data', sprintf('group%dk', atlas_res));
        cd(groupSurfDir);
        
        for tt = 1:length(t.instOrder_all)
            for h = 1:2 % two hemispheres
                for s = 1:length(sn)
                    %%% make the group metric file for each contrasts
                    infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)},sprintf('%s.%s.con_transition_%d-%s.func.gii', subj_name{sn(s)}, hemI{h}, t.instOrder_all(tt), con_vs));
                    if smooth
                        surfFile    = fullfile(groupSurfDir,sprintf('fs_LR.32k.%s.inflated.surf.gii',hemI{h}));
                        surf_smooth(infilenames{s},'surf',surfFile,'kernel',kernel); % smooth outfilenames - it will prefix an 's'
                    end
                    s_infilenames{s} = fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_transition_%d-%s.func.gii', subj_name{sn(s)}, hemI{h}, t.instOrder_all(tt), con_vs));
                end % sn
                outfilenames    = fullfile(groupSurfDir,sprintf('%s.con_transition_%d-%s.func.gii', hemI{h}, t.instOrder_all(tt), con_vs));
                summaryname     = fullfile(groupSurfDir,sprintf('%s.group.con_transition_%d-%s.func.gii', hemI{h}, t.instOrder_all(tt), con_vs));
                
                surf_groupGiftis(infilenames, 'outfilenames', {outfilenames}, 'groupsummary', summaryname, 'replaceNaNs', replaceNaN);
                if smooth % also save the smoothed versions
                    s_outfilenames    = fullfile(groupSurfDir,sprintf('s%s.con_transition_%d-%s.func.gii', hemI{h}, t.instOrder_all(tt), con_vs));
                    s_summaryname     = fullfile(groupSurfDir,sprintf('s%s.group.con_transition_%d-%s.func.gii', hemI{h}, t.instOrder_all(tt), con_vs));
                    surf_groupGiftis(s_infilenames, 'outfilenames', {s_outfilenames}, 'groupsummary', s_summaryname, 'replaceNaNs', replaceNaN);
                end % if smooth
                
                % delet the smoothed contrasts for each subject
                for s = 1:length(sn)
                    delete(fullfile(baseDir, experiment, wbDir, sprintf('glm%d', glm), subj_name{sn(s)}, sprintf('s%s.%s.con_transition_%d-%s.func.gii', subj_name{sn(s)}, hemI{h}, t.instOrder_all(tt), con_vs)));
                end
                
                fprintf('******************** group average contrast transition for %d vs %s is created! ********************\n\n', t.instOrder_all(tt), con_vs);
            end % hemi(h)
        end % tt (transitions)    
    case 'SURF:mdtb:noiseCeiling_get_utransitions_names_data'
        % Gets the data for calculating the noise ceiling
        % Uses the case for noise ceiling calculations
        % Creates surface gifti files for noise ceilings and saves them in
        % the workbench group directory
        % Example: sc1_sc2_mdtb_task_switch('SURF:mdtb:noiseCeiling_get_utransitions_names_data')
        
        experiment_num = 1;      %% sc1 or sc2;
        glm        = 7;
        atlas_res  = 32;     %% set it to 32 or 64
        smooth     = 1;      %% use the data with smoothing or without smoothing
        kernel     = 1;      %% the kernel used to smooth the data
        con_vs     = 'rest'; %% contrasts were calculated vs 'rest' or 'average'
        
        vararginoptions(varargin, {'experiment_num', 'glm', 'smooth', 'con_vs', 'smooth', 'kernel'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directories
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        wbDir  = 'surfaceWB';
        
        % load in all_transMat (contains task_before and task_after names)
        load(fullfile(glmDir, 'all_trans_task_names.mat'));
        
        nTrans = size(all_trans, 1);
        
        % Load in data
        for h = 1:2 % two hemispheres
            % group gifti for each task transition
            for trans=1:nTrans
                switch smooth % use data with or without smoothing
                    case 1
                        T = gifti(fullfile(baseDir, experiment, wbDir, 'data',...
                            sprintf('group%dk', atlas_res), sprintf('s%s.con_transition_%s_%s-%s.func.gii', hemI{h}, all_trans{trans, 1}, all_trans{trans, 2}, con_vs)));
                    case 0
                        T = gifti(fullfile(baseDir, experiment, wbDir, 'data',...
                            sprintf('group%dk', atlas_res), sprintf('%s.con_transition_%s_%s-%s.func.gii', hemI{h}, all_trans{trans, 1}, all_trans{trans, 2}, con_vs)));

                end % switch smooth
                groupData(:,:,trans) = T.cdata;
            end % trans
            
            % calculate low and high noise ceilings for hemI{h}
            [noise_low, noise_high] = sc1_sc2_mdtb_task_switch('SURF:mdtb:noiseCeiling_calculate', groupData);
            
            % creates surface gifti files for the noise ceilings
            G_high = surf_makeFuncGifti(noise_high, 'anatomicalStruct', hemName{h}, 'columnNames', {'noiseCeiling_high'});
            G_low  = surf_makeFuncGifti(noise_low, 'anatomicalStruct', hemName{h}, 'columnNames', {'noiseCeiling_low'});
            
            % save the gifti files
            switch smooth
                case 1 % with smoothing 
                    save(G_high, fullfile(baseDir, experiment, wbDir, 'data',...
                        sprintf('group%dk', atlas_res), sprintf('s%s.transition_names-%s_noiseCeiling_high.func.gii', hemI{h}, con_vs)));
                    save(G_low, fullfile(baseDir, experiment, wbDir, 'data',...
                        sprintf('group%dk', atlas_res), sprintf('s%s.transition_names-%s_noiseCeiling_low.func.gii', hemI{h}, con_vs)));
                case 0 % wihtout smoothing
                    save(G_high, fullfile(baseDir, experiment, wbDir, 'data',...
                        sprintf('group%dk', atlas_res), sprintf('%s.transition_names-%s_noiseCeiling_high.func.gii', hemI{h}, con_vs)));
                    save(G_low, fullfile(baseDir, experiment, wbDir, 'data',...
                        sprintf('group%dk', atlas_res), sprintf('%s.transition_names-%s_noiseCeiling_low.func.gii', hemI{h}, con_vs)));
            end % switch smooth
            
            fprintf('******************** low and high noise ceilings created for %s ********************', hemI{h});
        end % h
    case 'SURF:mdtb:noiseCeiling_get_transitions_id_data'
        % Gets the data for calculating the noise ceiling
        % Uses the case for noise ceiling calculations
        % Creates surface gifti files for noise ceilings and saves them in
        % the workbench group directory
        % Example: sc1_sc2_mdtb_task_switch('SURF:mdtb:noiseCeiling_get_transitions_id_data')
        
        experiment_num = 1;      %% sc1 or sc2;
        glm        = 8;
        atlas_res  = 32;     %% set it to 32 or 64
        con_vs     = 'average_2'; %% contrast vs 'rest' or 'average'
        smooth     = 1;      %% use the smoothed data?
        kernel     = 1;      %% smoothing kernel used
        
        vararginoptions(varargin, {'experiment_num', 'glm', 'atlas_res', 'con_vs', 'smooth', 'kernel'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directories
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        wbDir  = 'surfaceWB';
        
        % load in all_transMat (contains task_before and task_after names)
        load(fullfile(glmDir, 'trans_info.mat'));
        
        nTrans = length(t.instOrder_all);
        
        % Load in data
        for h = 1:2 % two hemispheres
            % group gifti for each task transition
            for trans=1:nTrans
                switch smooth % use the data with or without smoothing
                    case 1 % with smoothing
                        T = gifti(fullfile(baseDir, experiment, wbDir, 'data', sprintf('group%dk', atlas_res),...
                            sprintf('s%s.con_transition_%d-%s.func.gii', hemI{h}, t.instOrder_all(trans), con_vs)));
                    case 0 % without smoothing
                        T = gifti(fullfile(baseDir, experiment, wbDir, 'data', sprintf('group%dk', atlas_res),...
                            sprintf('%s.con_transition_%d-%s.func.gii', hemI{h}, t.instOrder_all(trans), con_vs)));
                end % switch smooth
                groupData(:,:,trans) = T.cdata;
            end % trans
            
            % calculate low and high noise ceilings for hemI{h}
            [noise_low, noise_high] = sc1_sc2_mdtb_task_switch('SURF:mdtb:noiseCeiling_calculate', groupData);
            
            % creates surface gifti files for the noise ceilings
            G_high = surf_makeFuncGifti(noise_high, 'anatomicalStruct', hemName{h}, 'columnNames', {'noiseCeiling_high'});
            G_low  = surf_makeFuncGifti(noise_low, 'anatomicalStruct', hemName{h}, 'columnNames', {'noiseCeiling_low'});
            
            % save the gifti files
            switch smooth % save giftis for data with or without smoothing
                case 1 % with smoothing
                    save(G_high, fullfile(baseDir, experiment, wbDir, 'data',...
                        sprintf('group%dk', atlas_res), sprintf('s%s.transition_id-%s_noiseCeiling_high.func.gii', hemI{h}, con_vs)));
                    save(G_low, fullfile(baseDir, experiment, wbDir, 'data',...
                        sprintf('group%dk', atlas_res), sprintf('s%s.transition_id-%s_noiseCeiling_low.func.gii', hemI{h}, con_vs)));
                case 0 % without smoothing
                    save(G_high, fullfile(baseDir, experiment, wbDir, 'data',...
                        sprintf('group%dk', atlas_res), sprintf('%s.transition_id-%s_noiseCeiling_high.func.gii', hemI{h}, con_vs)));
                    save(G_low, fullfile(baseDir, experiment, wbDir, 'data',...
                        sprintf('group%dk', atlas_res), sprintf('%s.transition_id-%s_noiseCeiling_low.func.gii', hemI{h}, con_vs)));
            end % switch smooth
            fprintf('******************** low and high noise ceilings created for %s hemi ********************\n\n', hemI{h});
        end % h
    case 'SURF:mdtb:noiseCeiling_calculate'
        % Calculates and returns the noise ceilings
        % Example: sc1_sc2_mdtb_task_switch('SURF:mdtb:noise_ceiling', groupData)
        
        Data = varargin{1}; %% group data for calculation of noise ceilings
        
        nNode  = size(Data,1);
        nSubj  = size(Data, 2);
        
        % preallocating the corr mats
        r_low  = zeros(nNode,nSubj);
        r_high = zeros(nNode,nSubj);

        % perform noise ceilings
        for s = 1:nSubj
            data1      = squeeze(Data(:, s, :));
            data2_low  = squeeze(nanmean(Data(:, ~ismember(1:nSubj, s), :), 2));
            data2_high = squeeze(nanmean(Data(:, :, :), 2));
%             % remove means
%             d1 = bsxfun(@minus,data1,mean(data1,2));
%             d2_l = bsxfun(@minus,data2_low,mean(data2_low,2));
%             d2_h = bsxfun(@minus,data2_high,mean(data2_high,2));
%             Pearson correlation formula
%             var1        = d1.^2;
%             var2_low    = d2_l.^2;
%             var2_high   = d2_h.^2;
%             cv_l        = var1.*var2_low;
%             cv_h        = var1.*var2_high;
%             r_low(:,i)  = sum(cv_l,2)./sqrt(sum(var1,2).*sum(var2_low,2));
%             r_high(:,i) = sum(cv_h,2)./sqrt(sum(var1,2).*sum(var2_high,2)); 
%             
            % This is not optimum!    
            for n = 1:nNode
                r_low(n, s)  = corr(data1(n, :)', data2_low(n, :)');
                r_high(n, s) = corr(data1(n, :)', data2_high(n, :)');
            end % n
        end % s (subjects)

        % overall noise ceiling across subjects
        noise_low   = nanmean(r_low,2);
        noise_high  = nanmean(r_high,2);
        
        % plot
        figure
        subplot(121);
        histogram(noise_low);
        subplot(122);
        histogram(noise_high);
        
        % output noise_low, noise_high
        varargout{1} = noise_low;
        varargout{2} = noise_high;
    case 'SURF:contrasts_single_gifti'
        % takes all the gifti files for the contrasts for each subject and
        % merge them all together in a single file. It can do it both on
        % individual level and group level.
        % Example: sc1_sc2_mdtb_backup('SURF:contrasts_single_gifti', 'sn', 3)
        
        sn          = returnSubjs;    %%
        glm         = 8;              %%
        experiment_num  = 2;              %%
        which       = 'task';         %%
        igroup      = 0;              %% do it for the group average files or for the individual subjects?
        con_vs      = 'average_1'; %% choose 'average_1', 'average_2', or 'rest', or 'rest_taskCon' (if glm7 and task)
%         atlas_res   = 32;           %% atlas resolution can either be 32 or 64
        replaceNaNs = 1;              %% set to 1 or 0 
        
        vararginoptions(varargin, {'sn', 'glm', 'experiment_num', 'which', 'con_vs', 'atlas_res', 'replaceNaNs'})
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc       = getrow(C, C.StudyNum == experiment_num);
        switch which
            case 'task' % task for glm8
                conNames = unique(Cc.taskNames);
            case 'cond' % condition for glm7
                conNames = unique(Cc.condNames);
        end %% do you want the group maps for tasks or conditions
        
        experiment = sprintf('sc%d', experiment_num);
        
        % in 'sc1_sc2_taskConds.txt' file, instruct is not coded as a
        % task/condition name. So I will have to add that to the list of
        % names
        conNames = ['Instruct'; conNames];
        
        % setting directories
        wbDir = fullfile(baseDir, experiment, 'surfaceWB');
        
        dircheck(wbDir);

        switch igroup % do it for group data or individual data
            case 1 % do it for the group data which is already smoothed?
            case 0 % do it for each subject separately
                for s = sn
                    dircheck(fullfile(wbDir, sprintf('glm%d', glm), subj_name{s}))
%                     subWbDir = fullfile(wbDir, subj_name{s});
                    % group all the contrast surface maps into a single
                    % file
                    %%% creating a single file for each hemisphere
                    for h = 1:2
                        for cc = 1:length(conNames)
                            infilenames{cc}   = fullfile(wbDir, sprintf('glm%d', glm), subj_name{s},sprintf('%s.%s.con_%s-%s.func.gii',...
                                subj_name{s}, hemI{h}, conNames{cc}, con_vs));
                            columnName{cc} = sprintf('%s-%s', conNames{cc}, con_vs);
                        end % cc (condition)
                        cd(fullfile(wbDir, sprintf('glm%d', glm), subj_name{s}));
                        outfilename = sprintf('%s.%s.con_%s-%s.func.gii', subj_name{s}, hemI{h}, which, con_vs);
                        surf_groupGiftis(infilenames, 'outfilenames', {outfilename}, 'outcolnames', columnName, 'replaceNaNs', replaceNaNs);
                        fprintf('a single gifti file for contrasts for %s hemi successfully created for %s\n', hemI{h}, subj_name{s})
                    end % h (hemi)  
                end % s (sn)
        end % switch igroup
        
    case 'SUIT:mdtb:suit_parcel2native'
        % maps each atlas of the suit into individual space
        % Example: sc1_sc2_mdtb_task_switch('SUIT:mdtb:suit_parcel2native', 'sn', [3])
        
        sn         = returnSubjs;
        parcelType = 'Buckner_7';                         %% set it to Buckner_7 or Buckner_17
        parcelDir  = fullfile(suitToolDir, 'atlasesSUIT'); %% directory where the nifti image is stored
        experiment_num = 1;
        
        vararginoptions(varargin, {'sn', 'parcelType', 'parcelDir', 'experiment_num'});
        
        experiment = sprintf('sc%d', experiment_num);
                
        spm('defaults','fmri');
        spm_jobman('initcfg');

        for s = sn
            J.Affine       = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'Affine_c_anatomical_seg1.mat'); %% the affine transformation already calculated using suit_normalize_dartel
            J.flowfield    = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'u_a_c_anatomical_seg1.nii');    %% the flow field already calculated
            J.resample{1}  = fullfile(parcelDir, sprintf('%sNetworks.nii', parcelType));                                          %% Buckner parcellation is suit space
            J.ref          = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'maskbrainSUITGrey.nii');        %% the reference image 
            
            suit_reslice_dartel_inv(J); %% creates sprintf('wi%sNetworks.nii', parcelType)
            
            fprintf('******************** %s transformation to native space for %s done! ********************\n\n', parcelType, subj_name{s}); 
        end % sn
    case 'SUIT:mdtb:reslice'
        % this case is used to reslice volumes into suit space
        % before you run all the cases for the group map, you have to run
        % this case to map all the contrast maps to suit space.
        % Example: sc1_sc2_mdtb_task_switch('SUIT:mdtb:reslice', 'sn', [3])
        
        sn         = returnSubjs;                   %% list of subjects
        experiment_num = 1;                      %% enter 1 for sc1 and 2 for sc2
        type       = 'con';                  %% enter the image you want to reslice to suit space
        glm        = 7;                      %% glm number
        mask       = 'cereb_prob_corr_grey'; %% the cerebellar mask to be used:'cereb_prob_corr_grey' or 'cereb_prob_corr' or 'dentate_mask'
        
        vararginoptions(varargin,{'sn', 'experiment_num', 'glm', 'type', 'mask'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directories
        glmDir        = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        imageDir      = fullfile(baseDir, experiment, 'imaging_data');
        glmSuitDir    = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm));
        dircheck(glmSuitDir);
        
        for s = sn
            dircheck(fullfile(glmSuitDir, subj_name{s}));
            glmSubjDir = fullfile(glmDir, subj_name{s});
            outDir     = fullfile(glmSuitDir, subj_name{s});
            
            job.subj.affineTr  = {fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, 'Affine_c_anatomical_seg1.mat')};
            job.subj.flowfield = {fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s},'u_a_c_anatomical_seg1.nii')};
            job.subj.mask      = {fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}, sprintf('%s.nii', mask))};
            job.vox            = [2 2 2];
            
            fprintf('******************** using dartel to reslice %s images to suit space for %s ********************\n', type, subj_name{s});
            
            switch type
                case 'beta'
                    images = 'beta_';
                    source = dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                    
                    job.subj.resample = {source.name};
                    suit_reslice_dartel(job);
                case 'con'  
                    % reslice all the contrasts (vs rest and vs average)
                    images = 'con';
                    source = dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                    
                    job.subj.resample = {source.name};
                    suit_reslice_dartel(job);
                case 'spmT'
                    images = 'spmT_';
                    source = dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                    
                    job.subj.resample = {source.name};
                    suit_reslice_dartel(job);
                case 'ResMS'
                    images = 'ResMS';
                    source = dir(fullfile(glmSubjDir,sprintf('*%s*',images))); % images to be resliced
                    cd(glmSubjDir);
                    
                    job.subj.resample = {source.name};
                    suit_reslice_dartel(job);
                case 'cerebellarGrey'
                    source = dir(fullfile(baseDir, experiment, suitDir,'anatomicals',subj_name{s},'c1anatomical.nii')); % image to be resliced
                    cd(fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{s}));
                    
                    job.subj.resample = {source.name};
                    suit_reslice_dartel(job);
                case 'time_series'
                    dircheck(fullfile(imageDir, subj_name{s}));
                    imageSubjDir = fullfile(imageDir, subj_name{s});
                    outDir       = fullfile(baseDir, experiment, suitDir, 'imaging_data_resliced', subj_name{s});
                    dircheck(outDir);
                    
                    images    = 'rrun';
                    
                    % there are 598 time points (images) that need to be
                    % resliced to the suit space.
                    for rr = 1:length(runLst)
                        for itp = 1:numTRs - numDummys
                            filenames{itp} = fullfile(imageSubjDir, sprintf('%s_%0.2d.nii,%d', images, rr, itp));
                        end % itp (time points)
                        job.subj.resample = filenames';
                        suit_reslice_dartel(job);
                    end % rr (runs)
            end % switch type
            
            if ~strcmp(type,'cerebellarGrey')
                if strcmp(type, 'time_series')
                    source = fullfile(imageSubjDir, '*wd*');
                    
                else
                    source = fullfile(glmSubjDir,'*wd*');
                end
                dircheck(fullfile(outDir));
                
                destination = fullfile(baseDir, experiment, suitDir, sprintf('glm%d',glm), subj_name{s});
                movefile(source, destination);
            end
            fprintf('******************** %s for glm %d have been resliced into suit space for %s ********************\n\n', type, glm, subj_name{s});
        end % sn
%         % plotting the flatmap
%         V = spm_vol(fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm), sprintf('wd%s_Instruct-rest.nii', type)));
%         D = suit_map2surf(V,'stats','nanmean');
%         suit_plotflatmap(D, 'cmap', colormap(jet(256)), 'cscale', [-3, 3]);
%         caxis([-3, 3]);
%         colorbar;
    case 'SUIT:mdtb:groupmap_con'
        % creates group average for the condition contrast maps.
        % you need to reslice all the images to suit space before running
        % this case
        % Example: sc1_sc2_mdtb_task_switch('SUIT:mdtb:groupmap_con', 'sn', [2, 3, 4, 6, 8, 9, 10, 12, 14, 15]);
        
        sn             = returnSubjs;        %% list of subjects
        experiment_num = 1;                  %% enter 1 for sc1 and 2 for sc2
        type           = 'con';              %% enter the image you want to reslice to suit space
        glm            = 8;                  %% glm number
        con_vs         = 'average_1';        %% is the contrast calculated vs 'rest' or 'average'
        which          = 'task';             %% you may choose 'cond' or 'task'
        
        vararginoptions(varargin,{'sn', 'experiment_num', 'glm', 'type', 'which'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc       = getrow(C, C.StudyNum == experiment_num);
        switch which
            case 'task' % task for glm8
                conNames = unique(Cc.taskNames);
            case 'cond' % condition for glm7
                conNames = unique(Cc.condNames);
        end %% do you want the group maps for tasks or conditions
        
        % in 'sc1_sc2_taskConds.txt' file, instruct is not coded as a
        % task/condition name. So I will have to add that to the list of
        % names
        conNames = ['Instruct'; conNames];
        
        experiment = sprintf('sc%d', experiment_num);
        
        % Setting directories
        glmSuitDir      = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm));
        glmSuitGroupDir = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm), 'group');
        dircheck(glmSuitDir);
        dircheck(glmSuitGroupDir);
        
        % preallocating!
        maps = []; %% the structure that will have alll?! the info I will need (hopefully :))
        for cc = 1:length(conNames)
            for s = 1:length(sn)                
                infilename{s} = fullfile(glmSuitDir, subj_name{sn(s)}, sprintf('wd%s_%s-%s.nii', type, conNames{cc}, con_vs));
                
                % load in individual cerebellar maps
                V_tmp = spm_vol(infilename{s});
                X     = spm_read_vols(V_tmp);
                X_vec = X(:);
                
                maps_tmp.conNames           = cellstr(conNames{cc});
                maps_tmp.conBaseName        = cellstr(con_vs);
                maps_tmp.subjectMaps{1}     = X;
                maps_tmp.subjectMapsVec     = X_vec';
                maps_tmp.subjectMapNames{1} = infilename{s};
                maps_tmp.glm                = 7;
                maps_tmp.sn                 = sn(s);
                maps_tmp.vol                = V_tmp;
                
                maps  = addstruct(maps, maps_tmp);
            end % sn
            outfilename = fullfile(glmSuitGroupDir, sprintf('group_%s_%s-%s.nii', type, conNames{cc}, con_vs));
            opt.dmtx    = 1;
            % calculate the average across subjects
            spm_imcalc(infilename', outfilename, 'nanmean(X)', opt);
            
            % saving the gifti file for group map
            V = spm_vol(fullfile(glmSuitGroupDir, sprintf('group_%s_%s-%s.nii', type, conNames{cc}, con_vs)));
            D = suit_map2surf(V,'stats','nanmean');
            
            G = surf_makeFuncGifti(D, 'anatomicalStruct', 'Cerebellum', 'columnNames', {sprintf('group_%s_%s-%s', type, conNames{cc}, con_vs)});
            
            save(G, fullfile(glmSuitGroupDir, sprintf('Cereb.group.con_%s-%s.func.gii', conNames{cc}, con_vs)));
            save(fullfile(glmSuitGroupDir, sprintf('indMaps_%s_%s-vs-%s.mat', type, conNames{cc}, con_vs)), 'maps', '-v7.3');
            fprintf('******************** %s group average for %s vs %s is created! ********************\n\n', type, conNames{cc}, con_vs);
            
        end % contrasts (cc)
        save(fullfile(glmSuitGroupDir, sprintf('indMaps_%s-vs-%s.mat', type, con_vs)), 'maps', '-v7.3');
    case 'SUIT:mdtb:groupmap_con_task'
        % creates group map for task contrasts
        % Example: sc1_sc2_mdtb_task_switch('SUIT:mdtb:groupmap_con_task', 'sn', [3]);
        
        sn         = returnSubjs;                   %% list of subjects
        experiment_num = 1;                      %% enter 1 for sc1 and 2 for sc2
        type       = 'con';                  %% enter the image you want to reslice to suit space
        glm        = 7;                      %% glm number
        con_vs     = 'rest';                 %% is the contrast calculated vs 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'experiment_num', 'glm', 'type', 'mask'});
        
        % load in task information
        C        = dload(fullfile(baseDir,'sc1_sc2_taskConds.txt'));
        Cc       = getrow(C, C.StudyNum == experiment_num);
        taskNames = unique(Cc.taskNames);
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directories
        glmSuitDir      = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm));
        glmSuitGroupDir = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm), 'group');
        
        dircheck(glmSuitDir);
        dircheck(glmSuitGroupDir);
        
        % preallocating!
        maps = []; %% the structure that will have alll?! the info I will need (hopefully :))
        for cc = 1:length(taskNames)
            for s = 1:length(sn)                
                infilename{s} = fullfile(glmSuitDir, subj_name{sn(s)}, sprintf('wd%s_%s-%s_taskCon.nii', type, taskNames{cc}, con_vs));
                
                % load in individual cerebellar maps
                V_tmp = spm_vol(infilename{s});
                X     = spm_read_vols(V_tmp);
                X_vec = X(:);
                
                maps_tmp.conNames           = cellstr(taskNames{cc});
                maps_tmp.conBaseName        = cellstr(con_vs);
                maps_tmp.subjectMaps{1}     = X;
                maps_tmp.subjectMapsVec     = X_vec';
                maps_tmp.subjectMapNames{1} = infilename{s};
                maps_tmp.glm                = 7;
                maps_tmp.sn                 = sn(s);
                maps_tmp.vol                = V_tmp;
                
                maps  = addstruct(maps, maps_tmp);
            end % sn
            outfilename = fullfile(glmSuitGroupDir, sprintf('group_%s_%s-%s_taskCon.nii', type, taskNames{cc}, con_vs));
            opt.dmtx    = 1;
            % calculate the average across subjects
            spm_imcalc(infilename', outfilename, 'nanmean(X)', opt);
            
            % saving the gifti file for group map
            V = spm_vol(fullfile(glmSuitGroupDir, sprintf('group_%s_%s-%s_taskCon.nii', type, taskNames{cc}, con_vs)));
            D = suit_map2surf(V,'stats','nanmean');
            
            G = surf_makeFuncGifti(D, 'anatomicalStruct', 'Cerebellum', 'columnNames', {sprintf('group_%s_%s-%s', type, taskNames{cc}, con_vs)});
            
            save(G, fullfile(glmSuitGroupDir, sprintf('Cereb.group.con_%s-%s_taskCon.func.gii', taskNames{cc}, con_vs)));
            save(fullfile(glmSuitGroupDir, sprintf('indMaps_%s_%s-vs-%s.mat', type, taskNames{cc}, con_vs)), 'maps', '-v7.3');
            fprintf('******************** %s group average for %s vs %s is created! ********************\n\n', type, taskNames{cc}, con_vs);
        end % contrasts (cc)
        save(fullfile(glmSuitGroupDir, sprintf('indMaps_%s-vs-%s.mat', type, con_vs)), 'maps', '-v7.3');
    case 'SUIT:mdtb:groupmap_con_utransitions_names'
        % creates group map for unique transitions with task names
        % Example: sc1_sc2_mdtb_task_switch('SUIT:mdtb:groupmap_con_utransitions_names', 'sn', [3])
        
        sn         = returnSubjs;   %% list of subjects
        experiment_num = 1;      %% enter 1 for sc1 and 2 for sc2
        glm        = 7;      %% glm number
        con_vs     = 'rest'; %% contrast was calculated against 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'atlas_res', 'experiment_num', 'glm', 'replaceNaN', 'convs', 'smooth', 'kernel'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directories
        glmDir          = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm)); % where the all_trans matrix is stored
        glmSuitDir      = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm));
        glmSuitGroupDir = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm), 'group');
        
        dircheck(glmSuitDir);
        dircheck(glmSuitGroupDir);
        
        % load in all_trans.mat for the task names before and after
        load(fullfile(glmDir, 'all_trans_task_names.mat'));
        
        % preallocating!
        maps = []; %% the structure that will have alll?! the info I will need (hopefully :))
        for tt = 1:size(all_trans, 1)
            for s = 1:length(sn)
                infilename{s} = fullfile(glmSuitDir, subj_name{sn(s)}, sprintf('wdcon_transition_%s_%s-%s', all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
                
                % load in individual cerebellar maps
                V_tmp = spm_vol(infilename{s});
                X     = spm_read_vols(V_tmp);
                X_vec = X(:);
                
                maps_tmp.conNames           = cellstr(sprintf('%s_%s', all_trans{tt, 1}, all_trans{tt, 2}));
                maps_tmp.conBaseName        = cellstr(con_vs);
                maps_tmp.subjectMaps{1}     = X;
                maps_tmp.subjectMapsVec     = X_vec';
                maps_tmp.subjectMapNames{1} = infilename{s};
                maps_tmp.glm                = 7;
                maps_tmp.sn                 = sn(s);
                maps_tmp.vol                = V_tmp;
                
                maps  = addstruct(maps, maps_tmp);
            end % sn
            outfilename = fullfile(glmSuitGroupDir, sprintf('group_con_transition_%s_%s-%s.nii', all_trans{tt, 1}, all_trans{tt, 2}, con_vs));
            opt.dmtx    = 1;
            % calculate the average across subjects
            spm_imcalc(infilename', outfilename, 'nanmean(X)', opt);
            
            % saving the gifti file for group map
            V = spm_vol(fullfile(glmSuitGroupDir, sprintf('group_con_transition_%s_%s-%s.nii', all_trans{tt, 1}, all_trans{tt, 2}, con_vs)));
            D = suit_map2surf(V,'stats','nanmean');
            
            G = surf_makeFuncGifti(D, 'anatomicalStruct', 'Cerebellum', 'columnNames', {sprintf('group_con_transition_%s_%s-%s', all_trans{tt, 1}, all_trans{tt, 2}, con_vs)});
            
            save(G, fullfile(glmSuitGroupDir, sprintf('Cereb.group.con_transition_%s_%s-%s.func.gii', all_trans{tt, 1}, all_trans{tt, 2}, con_vs)));
            save(fullfile(glmSuitGroupDir, sprintf('indMaps_con_transitions_%s_%s-vs-%s.mat', all_trans{tt, 1}, all_trans{tt, 2}, con_vs)), 'maps', '-v7.3');
            fprintf('******************** con group average for transitions from %s to %s vs %s is created! ********************\n\n', all_trans{tt, 1}, all_trans{tt, 2}, con_vs);
        end % tt 
    case 'SUIT:mdtb:groupmap_con_transitions_id'
        % creates group map for all the transitions (256 for sc1)
        % creates group map for unique transitions with task names
        % Example: sc1_sc2_mdtb_task_switch('SUIT:mdtb:groupmap_con_transitions_id', 'sn', [3])
        
        sn         = returnSubjs;   %% list of subjects
        experiment_num = 1;             %% enter 1 for sc1 and 2 for sc2
        glm        = 7;             %% glm number
        con_vs     = 'rest';        %% contrast was calculated against 'rest' or 'average'
        
        vararginoptions(varargin,{'sn', 'atlas_res', 'experiment_num', 'glm', 'replaceNaN', 'convs', 'smooth', 'kernel'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directories
        glmDir          = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm)); % where the all_trans matrix is stored
        glmSuitDir      = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm));
        glmSuitGroupDir = fullfile(baseDir, experiment, suitDir, sprintf('glm%d', glm), 'group');
        
        dircheck(glmSuitDir);
        dircheck(glmSuitGroupDir);
        
        % load in all_trans.mat for the task names before and after
        load(fullfile(glmDir, 'trans_info.mat'));
        
        % preallocating!
        for tt = 1:length(t.instOrder_all)           
            maps = []; %% the structure that will have alll?! the info I will need (hopefully :))
            for s = 1:length(sn)
                infilename{s} = fullfile(glmSuitDir, subj_name{sn(s)}, sprintf('wdcon_transition_%d-%s.nii', t.instOrder_all(tt), con_vs));
                
                % load in individual cerebellar maps
                V_tmp = spm_vol(infilename{s});
                X     = spm_read_vols(V_tmp);
                
                X_vec = X(:);
                
                maps_tmp.conNames           = cellstr(sprintf('transition_%d', t.instOrder_all(tt)));
                maps_tmp.conBaseName        = cellstr(con_vs);
                maps_tmp.subjectMaps{1}     = X;
                maps_tmp.subjectMapsVec     = X_vec';
                maps_tmp.subjectMapNames{1} = infilename{s};
                maps_tmp.glm                = 7;
                maps_tmp.sn                 = sn(s);
                maps_tmp.vol                = V_tmp;
                
                maps  = addstruct(maps, maps_tmp);
            end % sn            
            outfilename = fullfile(glmSuitGroupDir, sprintf('group_con_transition_%d-%s.nii', t.instOrder_all(tt), con_vs));
            opt.dmtx    = 1;
            % calculate the average across subjects
            spm_imcalc(infilename', outfilename, 'nanmean(X)', opt);
            
            % saving the gifti file for group map
            % saving the gifti file for group map
            V = spm_vol(fullfile(glmSuitGroupDir, sprintf('group_con_transition_%d-%s.nii', t.instOrder_all(tt), con_vs)));
            D = suit_map2surf(V,'stats','nanmean');
            
            G = surf_makeFuncGifti(D, 'anatomicalStruct', 'Cerebellum', 'columnNames', {sprintf('group_con_transition_%d-%s', t.instOrder_all(tt), con_vs)});
            
            save(G, fullfile(glmSuitGroupDir, sprintf('Cereb.group.con_transition_%d-%s.func.gii', t.instOrder_all(tt), con_vs)));
            
            save(fullfile(glmSuitGroupDir, sprintf('indMaps_con_transitions_id_%d-vs-%s.mat', t.instOrder_all(tt), con_vs)), 'maps', '-v7.3');
            
            fprintf('******************** con group average for transition %d vs %s is created! ********************\n\n', t.instOrder_all(tt), con_vs);
        end % tt    
    case 'SUIT:mdtb:noiseCeiling_get_transitions_id_data'
        % calculating noise ceilings for the cerebellum
        % Example: sc1_sc2_mdtb_task_switch('SUIT:mdtb:noiseCeiling_get_transitions_id_data')
        
        experiment_num = 1;      %% sc1 or sc2;
        glm        = 8;
        con_vs     = 'rest'; %% contrast vs 'rest' or 'average'
%         smooth     = 1;      %% use the smoothed data?
%         kernel     = 1;      %% smoothing kernel used
        
        vararginoptions(varargin, {'experiment_num', 'glm', 'atlas_res', 'con_vs', 'smooth', 'kernel'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % setting directories
        glmDir = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        suitGroupDir  = fullfile(baseDir, experiment, 'suit', sprintf('glm%d', glm), 'group');
        
        % load in all_transMat (contains task_before and task_after names)
        load(fullfile(glmDir, 'trans_info.mat'));
        
        nTrans = length(t.instOrder_all);
        
        % Load in data
        % suitmaps for each task transition
        for trans=1:nTrans
            fprintf('load in group data for transition %0.3d\n', trans)
            load(fullfile(suitGroupDir, sprintf('indMaps_con_transitions_id_%d-vs-%s', trans, con_vs)));
            
%             tmp = maps.subjectMaps; % this is a temp variable that will have the data for all the subjects.
%             for s = 1:length(maps.sn)
%                 groupData(:, :, :, s, trans) = tmp{s};
%             end
            groupData2(:, :, trans) = maps.subjectMapsVec';
        end % trans
        
        % calculate low and high noise ceilings for hemI{h}
%         [noise_low, noise_high] = sc1_sc2_mdtb_task_switch('SUIT:mdtb:noiseCeiling_calculate', groupData);
        [noise_low, noise_high] = sc1_sc2_mdtb_task_switch('SURF:mdtb:noiseCeiling_calculate', groupData2);
        
        myTmp     = maps.vol;
        V_low     = myTmp(1);
        V_low     = rmfield(V_low, 'fname');
        V_low.dat = noise_low;
        
        D_low = suit_map2surf(V_low,'stats','nanmean');
        
        V_high     = myTmp(1);
        V_high     = rmfield(V_high, 'fname');
        V_high.dat = noise_high;
        
        D_high = suit_map2surf(V_high,'stats','nanmean');
        
        
        % creates surface gifti files for the noise ceilings
        G_high = surf_makeFuncGifti(D_high, 'anatomicalStruct', 'Cerebellum', 'columnNames', {'noiseCeiling_high'});
        G_low  = surf_makeFuncGifti(D_low, 'anatomicalStruct', 'Cerebellum', 'columnNames', {'noiseCeiling_low'});
        
        % save the gifti files
        switch smooth % save giftis for data with or without smoothing
            case 0 % without smoothing
                save(G_high, fullfile(suitGroupDir, sprintf('Cereb.transition_id-%s_noiseCeiling_high.func.gii', con_vs)));
                save(G_low, fullfile(suitGroupDir, sprintf('Cereb.transition_id-%s_noiseCeiling_low.func.gii', con_vs)));
        end % switch smooth
        fprintf('******************** low and high noise ceilings created for %s hemi ********************\n\n', hemI{h});
    case 'SUIT:mdtb:noiseCeiling_calculate'
        % Calculates and returns the noise ceilings
        % Example: sc1_sc2_mdtb_task_switch('SURF:mdtb:noise_ceiling', groupData)
        
        Data = varargin{1}; %% group data for calculation of noise ceilings
        
%         [nx, ny, nz]  = size(Data{});
        nSubj  = size(Data, 4);
        nx = size(Data, 1);
        ny = size(Data, 2);
        nz = size(Data, 3);
        
        % preallocating the corr mats
        r_low  = zeros(nx, ny, nz, nSubj);
        r_high = zeros(nx, ny, nz, nSubj);

        % perform noise ceilings
        for s = 1:nSubj
            data1      = squeeze(Data(:, :, :, s, :));
            data2_low  = squeeze(nanmean(Data(:, :, :, ~ismember(1:nSubj, s), :), 4));
            data2_high = squeeze(nanmean(Data(:, :, :, :, :), 4));
%             % remove means
%             d1 = bsxfun(@minus,data1,mean(data1,2));
%             d2_l = bsxfun(@minus,data2_low,mean(data2_low,2));
%             d2_h = bsxfun(@minus,data2_high,mean(data2_high,2));
%             Pearson correlation formula
%             var1        = d1.^2;
%             var2_low    = d2_l.^2;
%             var2_high   = d2_h.^2;
%             cv_l        = var1.*var2_low;
%             cv_h        = var1.*var2_high;
%             r_low(:,i)  = sum(cv_l,2)./sqrt(sum(var1,2).*sum(var2_low,2));
%             r_high(:,i) = sum(cv_h,2)./sqrt(sum(var1,2).*sum(var2_high,2)); 
%             
            % This is not optimum!    
            for ix = 1:nx
                for iy = 1:ny
                    for iz = 1:nz
                        r_low(ix, iy, iz, s)  = corr(squeeze(data1(ix, iy, iz, :)), squeeze(data2_low(ix, iy, iz, :)));
                        r_high(ix, iy, iz, s) = corr(squeeze(data1(ix, iy, iz, :)), squeeze(data2_high(ix, iy, iz, :)));
                    end % iz
                end % iy
            end % ix
        end % s (subjects)

        % overall noise ceiling across subjects
        noise_low   = nanmean(r_low,4);
        noise_high  = nanmean(r_high,4);
        
        % plot
        figure
        subplot(121);
        histogram(noise_low(:));
        subplot(122);
        histogram(noise_high(:));
        
        % output noise_low, noise_high
        varargout{1} = noise_low;
        varargout{2} = noise_high;
      
    case 'Visualize:mdtb:suit'
        % takes in a vector (or map) for a subject and transfer it to suit space and plot it
        % on the flatmap
        % Example: sc1_sc2_mdtb_task_switch('Visualize:mdtb:suit', corrmap);

        subj       = 2;           % the suit anatomical data for the subject is used in the mapping 
        experiment_num = 1;
        data       = varargin{1}; % the vector map you want to transfer to suit flatmap
        
        vararginoptions(varargin, {'sn', 'experiment_num'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        % Determine the voxels we want to resample in SUIT space
        V = spm_vol(fullfile(baseDir, experiment,suitDir,'anatomicals','cerebellarGreySUIT.nii'));
        X = spm_read_vols(V);
        
        grey_threshold = 0.1; % gray matter threshold
        
        linIn1     = find(X > grey_threshold);
        [i1,j1,k1] = ind2sub(V.dim,linIn1');
        [x1,y1,z1] = spmj_affine_transform(i1, j1, k1, V.mat);
        
        %%% Determine voxel locations from the original ROI
        load(fullfile(baseDir, experiment, regDir, 'data', subj_name{subj}, 'regions_cerebellum_grey.mat')); % 'regions' are defined in 'ROI_define'
        
        Vmask      = spm_vol(fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{subj}, 'maskbrainSUITGrey.nii')); %% cerebellum grey matter mask
        [i3,j3,k3] = spmj_affine_transform(R{1}.data(:,1),R{1}.data(:,2),R{1}.data(:,3),inv(Vmask.mat));
        linIn3     = sub2ind(Vmask.dim, round(i3), round(j3), round(k3));
        
        %%% transform SUIT coords into anatomical space of the individual
        flowfield    = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{subj}, 'u_a_c_anatomical_seg1.nii');
        affine       = fullfile(baseDir, experiment, suitDir, 'anatomicals', subj_name{subj}, 'Affine_c_anatomical_seg1.mat');
        [Def, Aff]   = spmdefs_get_dartel(flowfield, affine);
        [x2, y2, z2] = spmdefs_transform(Def, Aff, x1, y1, z1);
        [i2, j2, k2] = spmj_affine_transform(x2, y2, z2, inv(Vmask.mat));
        
        %%% resample the weights into SUIT space for the lp vector
        
        Vout             = Vmask;
        Vout.dat         = zeros(Vout.dim);
        Vout.dat(linIn3) = data;
        Vout.dt          = [64 0];
        Vout.pinfo       = [1 0 0]';
        
        DataSUIT(ip,:) = spm_sample_vol(Vout,i2,j2,k2,1);
        
        V.dat   = zeros(V.dim);
        Vres(ip) = V;
        Vres(ip).dat(linIn1) = DataSUIT(ip,:);  % Offset by one to account for 1 being medial wall
        Vres(ip).fname       = sprintf('data_%2.2d.nii',ip);
        Vres(ip).pinfo       = [1 0 0]';
        
        myV = Vres(ip);
        
        %%% Now map the Lp vector to surface-based representation
        D = suit_map2surf(myV,'stats','nanmean');
        
        figure; suit_plotflatmap(D , 'cmap', colormap(jet(256)), 'cscale', [min(D(:)), max(D(:))]);
        caxis([min(D(:)), max(D(:))]);
        colorbar;
        title(sprintf('cerebellar regions with high corr with cortex parcel %0.2d for %s', ip, subj_name{subj}));
        
        keyboard;
        
    case 'Houskeeping:renameSPM'     % rename SPM directories
        % Example: sc1_sc2_mdtb_task_switch('HOUSEKEEPING_renameSPM', 'experiment_num', 2, 'glm', 8)
        
        sn             = returnSubjs;
        experiment_num = 1;
        glm            = 8;
        
        vararginoptions(varargin, {'sn', 'experiment_num', 'glm'});
        
        experiment = sprintf('sc%d', experiment_num);
        
        glmDir     = fullfile(baseDir, experiment, sprintf('GLM_firstlevel_%d', glm));
        imagingDir = fullfile(baseDir, 'sc1/imaging_data'); %% all the imaging files are in sc1
        
        for s = sn
            fprintf('******************** changing directories for %s ********************\n', subj_name{s});
            newGLMDir   = fullfile(glmDir,subj_name{s});
            newRawDir   = fullfile(imagingDir,subj_name{s});
            
            % load SPM file
            load(fullfile(newGLMDir, 'SPM.mat'));
            
            SPM         = spmj_move_rawdata(SPM,newRawDir);
            SPM.swd     = fullfile(newGLMDir);
            save(fullfile(newGLMDir,'SPM.mat'),'SPM','-v7.3');
            varargout{1} = SPM;
        end % s (sn)
end
end

% Local functions
function dircheck(dir)
if ~exist(dir,'dir')
    warning('%s doesn''t exist. Creating one now. You''re welcome! \n',dir);
    mkdir(dir);

end
end