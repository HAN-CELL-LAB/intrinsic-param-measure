# Extract intrinsic parameters from whole-cell recordings 

## Installation and dependencies

- Currently was run on `MATLAB R2021a`, previous versions may also work 
- Required toolboxes from `MATLAB` (just go to `Apps -> Get More Apps` if not already installed)
    - `Signal Processing Toolbox`
    - `Curve Fitting Toolbox`
- Then run `matlab_install.m` (check the above dependencies)
    - Will prompt whether to save paths, generally say `No` so as to not have conflict from your other project 
- Usage: before using the files, either
    - Run `matlab_startup.m` in `MATLAB`
    - Or just add necessary paths, with `addpath(genpath('matlab-src'))` at the beginning of your scripts 

## Demo data 

- The `data` folder contains a `demo-data.tar.gz`, for demonstration purposes
- Just extract the compressed file inside, e.g. in `bash/shell`
    ``` bash
    cd data/ 
    tar xvf demo-data.tar.gz
    cd .. # go up to fiddle the repo
    ```

## Conversion to necessary `mat` file

- In `MATLAB`, run `itx2mat` on command window (see also `A_convert_itx2mat.m`) to convert `itx` files to `mat` files. 
    - The structure of the `itx` file should be from the `Patchmaster -> Export -> Igor ITX`. 
    - And this script could only handle one channel for now in the `itx` file (though can be adapted to more). 
- Each `mat` file would need to include a struct array variable, 
    which I call `recordings` in the output of `itx2mat`:

    ``` text 
    recordings: struct array containing each sweep, of only 1 channel 
                each struct contains these fields: 
        + id:   unique segment id
        + dt:   time sampling interval (second) 
        + Fs:   sampling rate (Hz)
        + data: vector containing recorded data 
    ```
- The purpose of this is to detect membrane potential parameters like spike numbers, thresholds and AHP, 
    hence the data here are assumed to be membrane potential

## Usage

- For detailed configurations and usage, see `B_demo.m`
- Of course, for different cells and different experimental settings, one would have to adjust configurations
- The configurations in that script are defined in the 4 structs: 
    - `stim_info`
    - `spike_params` 
    - `AP_params` 
    - `mAHP_params` 
- The functions needed are:

    ``` MATLAB
    % t: time in ms
    % Vm: membrane potential vector in mV
    % dt: sampling time in ms

    % only spiking rate and AP locations
    spike_properties = calculate_spike_properties(t, Vm, spike_params, stim_info); 

    % AP properties including threshold, AP widths, estimated rise/fall times and fAHP
    % do this for each spike if there are spikes from `spike_properties` output 
    AP_properies = arrayfun(@(AP_ind) calculate_AP_properties(Vm, AP_ind, dt, AP_params), ...
        spike_properties.spike_inds, 'uni', 1);
    AP_properies = structarray_to_struct(AP_properies);

    % mAHP estimations  
    ind_start_mAHP = ceil((stim_info.pulse_off_time)/dt); % need to know when stim is off 
    mAHP_properties = calculate_mAHP(Vm, ind_start_mAHP, dt, mAHP_params);
    mAHP_methods = fieldnames(mAHP_properties);

    % Additional passive properties
    Vrest_properties = calculate_Vrest(t, Vm, stim_info, Vrest_calc_time_window);
    Rin_properties = calculate_Rin(t, Vm, stim_info, Rin_calc_time_window);
    ```
