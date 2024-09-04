%% File Directories
rawDataDir = 'Two-ColorDiFC_RawData\';
saveToDir = 'Two-colorDiFC_Data\';

%% Control Mouse Data
ctrlDataFiles = {'NSG_Control_1' 'NSG_Control_2' 'NSG_Control_3' 'NSG_Control_4'};
processMouse(rawDataDir, ctrlDataFiles, saveToDir)

%% MM Phantom Data (Single Cells)
MMPhantom_GFP = 'Phantom_MM_GFP';
MMPhantom_tdTomato = 'Phantom_MM_TDT';
processMouse(rawDataDir, {MMPhantom_GFP MMPhantom_tdTomato}, saveToDir);

%% MM Mouse Data
MM_files = {'MM_Mouse1_1' 'MM_Mouse1_2' 'MM_Mouse1_3' 'MM_Mouse1_4' 'MM_Mouse1_5' 'MM_Mouse1_6' 'MM_Mouse1_7' 'MM_Mouse1_8'...
    'MM_Mouse2_1' 'MM_Mouse2_2' 'MM_Mouse2_3' 'MM_Mouse2_4' 'MM_Mouse2_5'...
    'MM_Mouse3_1' 'MM_Mouse3_2' 'MM_Mouse3_3' 'MM_Mouse3_4'...
    'MM_Mouse4_1' 'MM_Mouse4_2' 'MM_Mouse4_3' 'MM_Mouse4_4'};
processMouse(rawDataDir, MM_files, saveToDir);

%% 4T1 Phantom Data (Clusters)
BC4T1_Phantom_2cCluster = {'Phantom_4T1_2Color_Clusters'};
BC4T1_Phantom_both_cluster = {'Phantom_4T1_Both_GFP_TDT_1' 'Phantom_4T1_Both_GFP_TDT_2' 'Phantom_4T1_Both_GFP_TDT_3' 'Phantom_4T1_Both_GFP_TDT_4' };
BC4T1_Phantom_GFP = {'Phantom_4T1_GFP_1' 'Phantom_4T1_GFP_2' 'Phantom_4T1_GFP_3'};
BC4T1_Phantom_TDT = {'Phantom_4T1_TDT_1' 'Phantom_4T1_TDT_2' 'Phantom_4T1_TDT_3' 'Phantom_4T1_TDT_4'};

processMouse(rawDataDir, BC4T1_Phantom_2cCluster, saveToDir);
processMouse(rawDataDir, BC4T1_Phantom_both_cluster, saveToDir);
processMouse(rawDataDir, BC4T1_Phantom_GFP, saveToDir);
processMouse(rawDataDir, BC4T1_Phantom_TDT, saveToDir);

%%
function [] = processMouse(data_path, data_files, saveToDir)
    for i = 1:length(data_files)
        saveDir = [data_path data_files{i}];
        Two_Color_DiFC_Processing(saveDir, 5, 5, 0, 1, saveToDir);
    end
end