function [allData] = load2ColorData(data_files, saveToDir)
    scan_lengths = zeros(length(data_files),1);
    num_fwd_gfp = scan_lengths;
    num_fwd_tdt = scan_lengths;
    num_rev_gfp = scan_lengths;
    num_rev_tdt = scan_lengths;
    
    num_OneFP = scan_lengths;
    num_TwoFP = scan_lengths;
    fwd_pk_ampl_gfp = [];
    fwd_pk_ampl_tdt = [];
    rev_pk_ampl_gfp = [];
    rev_pk_ampl_tdt = [];
    OneFP_ampl = [];
    OneFP_color = [];
    TwoFP_ampl = [];
    TwoFP_color = [];
    
    for i = 1:length(data_files)
        load([saveToDir data_files{i}], 'out_dat');
        
        scan_lengths(i) = out_dat.scan_length_minutes;
        
        num_fwd_gfp(i) = out_dat.fwd_peaks_color1(1).count;
        num_fwd_tdt(i) = out_dat.fwd_peaks_color2(1).count;
        num_rev_gfp(i) = out_dat.rev_peaks_color1(1).count;
        num_rev_tdt(i) = out_dat.rev_peaks_color2(1).count;
        
        num_OneFP(i) = out_dat.probe1_1FP(1).count + out_dat.probe1_1FP(1).count;
        num_TwoFP(i) = out_dat.probe1_2FP(1).count + out_dat.probe2_2FP(1).count;
        
        fwd_pk_ampl_gfp = [fwd_pk_ampl_gfp; [out_dat.fwd_peaks_color1(1).pks out_dat.fwd_peaks_color1(2).pks]];
        fwd_pk_ampl_tdt = [fwd_pk_ampl_tdt; [out_dat.fwd_peaks_color2(1).pks out_dat.fwd_peaks_color2(2).pks]];
        rev_pk_ampl_gfp = [rev_pk_ampl_gfp; [out_dat.rev_peaks_color1(1).pks out_dat.rev_peaks_color1(2).pks]];
        rev_pk_ampl_tdt = [rev_pk_ampl_tdt; [out_dat.rev_peaks_color2(1).pks out_dat.rev_peaks_color2(2).pks]];
        
        OneFP_ampl = [OneFP_ampl; out_dat.OneFP_ampl];
        OneFP_color = [OneFP_color; out_dat.OneFP_color];
        TwoFP_ampl = [TwoFP_ampl; out_dat.TwoFP_ampl];
        TwoFP_color = [TwoFP_color; out_dat.TwoFP_color];
    end
    
    allData.scan_lengths = scan_lengths;
    allData.num_fwd_gfp = num_fwd_gfp;
    allData.num_fwd_tdt = num_fwd_tdt;
    allData.num_rev_gfp = num_rev_gfp;
    allData.num_rev_tdt = num_rev_tdt;
    allData.num_OneFP = num_OneFP;
    allData.num_TwoFP = num_TwoFP;
    allData.fwd_pk_ampl_gfp = fwd_pk_ampl_gfp;
    allData.fwd_pk_ampl_tdt = fwd_pk_ampl_tdt;
    allData.rev_pk_ampl_gfp = rev_pk_ampl_gfp;
    allData.rev_pk_ampl_tdt = rev_pk_ampl_tdt;
    allData.OneFP_ampl = OneFP_ampl;
    allData.OneFP_color = OneFP_color;
    allData.TwoFP_ampl = TwoFP_ampl;
    allData.TwoFP_color = TwoFP_color;
    
end

