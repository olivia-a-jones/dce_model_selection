%% Requires fm_toolbox on path
%% Requires SPM on path

%% Preprocess some DCE-MRI data.
data_dir = 'E:/DCE_Model_Selection/';
ids_file = fullfile(data_dir,'Patient_IDs.txt');
hct_file = fullfile(data_dir,'Hematocrit.txt');

id_list = strsplit(ids_file);
hct_list = strsplit(hct_file);

for person=1:length(id_list)
    id = char(id_list(person))
    participant_dir = fullfile(data_dir,id);
    vsize = [1.5 1.5 2]

    % Make sure VFA files are organised:
    % directory: data_dir/participant_id/vfa/
    % filenames: [flip_angle]_deg.nii, 

    %% VFA
    vfadir = fullfile(participant_dir, 'vfa')
    for fa = [2,6,10,15]
        deg = double(load_nii(fullfile(vfadir, [num2str(fa) 'deg.nii'])).img);
        deg_avg = mean(deg(:,:,:,2:6),4);
        deg_avg_nii = make_nii(deg_avg, vsize)
        save_nii(deg_avg_nii, fullfile(vfadir, [num2str(fa) 'deg.nii']))
        cd(rawdir)
    end

    % Make sure dce injection file is organised:
    % directory: data_dir/participant_id/dynamic_series/
    % filename: dce_inj.nii

    dyndir = fullfile(participant_dir, 'dynamic_series')
    deg12 = double(load_nii(fullfile(dyndir, 'dce_inj.nii')).img);
    deg12_avg = mean(deg12(:,:,:,2:6),4);
    deg12_avg_nii = make_nii(deg12_avg, vsize)
    save_nii(deg12_avg_nii, fullfile(vfadir, ['12deg.nii']).img)

   %% B1 map
    
    % Make sure multi TR files are organised:
    % directory: data_dir/participant_id/vfa/
    % filenames: TR=25ms, imgsTR1.nii ; TR=125ms, imgsTR2.nii

    imgsTR1 = double(load_nii(fullfile(vfadir, 'imgsTR1.nii')).img);
    imgsTR2 = double(load_nii(fullfile(vfadir, 'imgsTR2.nii')).img);
    
    TR1 = 25.0; % ms
    TR2 = 125.0; % ms
       
    sImgsTR1 = smooth3(imgsTR1,'box',[5 5 3]);
    sImgsTR2 = smooth3(imgsTR2,'box',[5 5 3]);
    
    r = imgsTR1./imgsTR2;
    sr = sImgsTR1./sImgsTR2;
    n = TR2/TR1;
    
    FA_map = acosd((r.*n - 1)./(n - r))/60; % B1 map
    sFA_map = acosd((sr.*n - 1)./(n - sr))/60; % smoothed B1 map
    
    sImgsTR2_nii = make_nii(sImgsTR2, vsizeB1);
    save_nii(sImgsTR2_nii,fullfile(vfadir,'sImgsTR2.nii'));
    centre_header_file((vfadir,'sImgsTR2.nii'));
    
    sImgsTR1_nii = make_nii(sImgsTR1, vsizeB1);
    save_nii(sImgsTR1_nii,(vfadir,'sImgsTR1.nii'));
    centre_header_file('sImgsTR1.nii');
    
    sB1_nii = make_nii(sFA_map*100, vsizeB1);
    save_nii(sB1_nii,fullfile(vfadir,'sB1.nii'));
    centre_header_file(fullfile(vfadir,'sB1.nii'));

    %% Make Madym XTR files  
    for i = 1:length(VFAs)
    filename = [num2str(VFAs(i)) 'deg.nii'];
    filepaths(i) = {fullfile(participant_dir, 'vfa', filename)};
    end
    
    run_madym_MakeXtr(...
    'T1_vols', [filepaths], ...
    'working_directory', participant_dir, ...
    'dynamic_basename', 'dynamic_series/rdyn_', ...
    'sequence_format', '%01u', ...Format for converting dynamic series index to string, eg %01u
    'sequence_start', 1, ...Start index for dynamic series file names
    'sequence_step', 1, ...Step between indexes of filenames in dynamic series
    'n_dyns', 160, ...
    'make_t1', 1, ...
    'make_dyn', 1, ...  
    'temp_res', 7.64, ... %either set this or specify dyn_times
    'TR', 2.4, ... 
    'FA', 10, ... 
    'VFAs', VFAs, ... 
    'dummy_run', 0);

    %% Make T1 map
    run_madym_T1(...
        'cmd_exe', [local_madym_root 'madym_T1'],...
        'T1_vols', [filepaths],... Cell array of variable flip angle file paths
  	    'TR', '2.4',... TR in msecs, required if directly fitting (otherwise will be taken from FA map headers);
        'method', 'VFA_B1',...T1 method to use to fit, see notes for options
        'B1_name', B1_input,...Path to B1 correction map
        'B1_correction', true, ... Apply B1 correction
        'B1_scaling', 100, ... Scaling factor to use with B1 map
        'output_dir', [fullfile(participant_dir,'t1_map')], ...Output path, will use temp dir if empty;
        'overwrite', true,...Set overwrite existing analysis in output dir ON
        'dummy_run', false ...Don't run any thing, just print the cmd we'll run to inspect
        );
    %% Unpack and realign dynamic series
    % Make sure dce injection file is organised:
    % directory: data_dir/participant_id/dynamic_series/
    % filename: dce_inj.nii

    % split dynamic series for alignment
    dce = double(load_nii(fullfile(dyndir,'dce_inj.nii')).img);
    dims = size(dce);
    for d = 1:dims(4)
        dce_dyn = dce(:,:,:,d);
        dyn_nii = make_nii(dce_dyn, vsize);
        dyn_nii_name = ['dyn_' num2str(d) '.nii'];
        save_nii(dyn_nii, fullfile(dyndir,dyn_nii_name));
        centre_header_file(fullfile(dyndir,dyn_nii_name));
    end
     % Align
    SPM_align(id);
    clear jobs;
end