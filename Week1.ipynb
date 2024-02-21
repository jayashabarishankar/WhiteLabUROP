import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import mofax as mofa

# Load the data for all phosphoproteomics and global proteomics
phosphoproteomics_data = pd.read_csv("data/Data/globalphospho(pST)datasets")
global_proteomics_data = pd.read_csv("data/Data/globalproteomicsdatasets")

# Perform analysis on files for individual runs
# Assuming "individual runs" data is available in separate files, adjust the file paths accordingly

# Load the bridge data for normalization?

# Separate individual runs from bridge data
individual_runs = [phosphoproteomics_data, global_proteomics_data] 
                   
for run_num, run_data in individual_runs:
                   
    run_data = run_data.set_index("Sample_ID") 
    bridge_data = bridge_data.set_index("Sample_ID")
    common_samples = run_data.index.intersection(bridge_data.index)
    run_data = run_data.loc[common_samples]
    bridge_data = bridge_data.loc[common_samples]
    
    combined_data = pd.concat([run_data, bridge_data], axis=1)
    
    scaler = StandardScaler()
    standardized_data = scaler.fit_transform(combined_data)
    

    mofa_result = mofa.fit_model(standardized_data)

    mofa.plot_factors(mofa_result)
    plt.title(f"MOFA Analysis - Individual Run {run_num}")
    plt.show()
