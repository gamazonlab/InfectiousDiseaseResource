#Using FinnGen and UKBB Predixcan Results as Imput, Processes Results Based on Trait and 
#Single Tissue of Interest to Obtain Observed Pvalues for QQ Plot R Script

import pandas as pd
import glob

#First read in all disease traits from BioVU study
files = glob.glob(r'#Absolute path to BioVU results file directory/*.txt')
print(files)
dfs=[]
for file in files:
    print(file)
    df = pd.read_csv(file, delimiter='\t', index_col=False)
    df['Trait']=file.split('\\')[7].rstrip('_SingleTissuePrediXcan.txt')
    dfs.append(df)
print(dfs)

#Subset each trait of interest into its own dataframe
influenza = dfs[6]
viral_pneumonia = dfs[5]
bacterial_pneumonia = dfs[4]
encephalitis = dfs[3]
meningitis = dfs[2]

#First trait influenza
#Import first FinnGen Trait of Interest Influenza
files2 = glob.glob(r'#Absolute path to all influenza Predixcan Results'+'/*.csv')
dfs2 = []
for file2 in files2:
    df2 = pd.read_csv(file2, index_col=False, header=0)
    df2['Tissue']=file2.split('\\')[7].rstrip('.csv')
    dfs2.append(df2)
print(dfs2)

#Subset to lung tissue
influenza_finn = pd.concat(dfs2)

influenza_finn = influenza_finn.loc[influenza_finn['Tissue'] == 'Lung']
influenza_lung = influenza.loc[influenza['tissue'] == 'Lung']

#Subset to Top BioVU Associations
genes = influenza_lung['gene']
influenza_finn2 = influenza_finn.loc[influenza_finn['gene_name'].isin(genes)]

#Write to csv files
influenza_finn2.to_csv(r'#Absolute path to output file location\Influenza_Lung_Finngen.csv', index=False, encoding='utf-8')

#Second trait viral pneumonia
#Import first FinnGen Trait of Interest Influenza
files2 = glob.glob(r'#Absolute path to all viral pneumonia Predixcan Results'+'/*.csv')
dfs2 = []
for file2 in files2:
    df2 = pd.read_csv(file2, index_col=False, header=0)
    df2['Tissue']=file2.split('\\')[7].rstrip('.csv')
    dfs2.append(df2)
print(dfs2)

#Subset to lung tissue
viralpneumo_finn = pd.concat(dfs2)
viralpneumo_finn = viralpneumo_finn.loc[viralpneumo_finn['Tissue'] == 'Lung']
viralpneumo_lung = viral_pneumonia.loc[viral_pneumonia['tissue'] == 'Lung']

#Subset to Top BioVU Associations
genes = viralpneumo_lung['gene']
viralpneumo_finn2 = viralpneumo_finn.loc[viralpneumo_finn['gene_name'].isin(genes)]

#Write to csv files
viralpneumo_finn2.to_csv(r'#Absolute path to output file location\Viral_Pneumonia_Lung_Finngen.csv', index=False, encoding='utf-8')

#Third trait bacterial pneumonia
#Import first FinnGen Trait of Interest Influenza
files2 = glob.glob(r'#Absolute path to all Bacterial Pneumonia Predixcan Results'+'/*.csv')
dfs2 = []
for file2 in files2:
    df2 = pd.read_csv(file2, index_col=False, header=0)
    df2['Tissue']=file2.split('\\')[7].rstrip('.csv')
    dfs2.append(df2)
print(dfs2)

#Subset to lung tissue
bactpneumo_finn = pd.concat(dfs2)
bactpneumo_finn = bactpneumo_finn.loc[bactpneumo_finn['Tissue'] == 'Lung']
bactpneumo_lung = bacterial_pneumonia.loc[bacterial_pneumonia['tissue'] == 'Lung']

#Subset to Top BioVU Associations
genes = bactpneumo_lung['gene']
bactpneumo_finn2 = bactpneumo_finn.loc[bactpneumo_finn['gene_name'].isin(genes)]

#Write to csv files
bactpneumo_finn2.to_csv(r'#Absolute path to output file location\Bacterial_Pneumonia_Lung_Finngen.csv', index=False, encoding='utf-8')

#Fourth trait encephalitis
#Import first FinnGen Trait of Interest Influenza
files2 = glob.glob(r'#Absolute path to all Encephalitis Predixcan Results'+'/*.csv')
dfs2 = []
for file2 in files2:
    df2 = pd.read_csv(file2, index_col=False, header=0)
    df2['Tissue']=file2.split('\\')[7].rstrip('.csv')
    dfs2.append(df2)
print(dfs2)

#Subset to Brain tissue
encephalitis_finn = pd.concat(dfs2)
brain = ['Brain_Anterior_cingulate_cortex_BA24','Brain_Caudate_basal_ganglia','Brain_Cerebellar_Hemisphere','Brain_Cerebellum','Brain_Cortex','Brain_Frontal_Cortex_BA9','Brain_Hippocampu','Brain_Hypothalamu','Brain_Nucleus_accumbens_basal_ganglia','Brain_Putamen_basal_ganglia']
brain2 = ['Brain_Anterior_cingulate_cortex_BA24','Brain_Caudate_basal_ganglia','Brain_Cerebellar_Hemisphere','Brain_Cerebellum','Brain_Cortex','Brain_Frontal_Cortex_BA9','Brain_Hippocampus','Brain_Hypothalamus','Brain_Nucleus_accumbens_basal_ganglia','Brain_Putamen_basal_ganglia']
encephalitis_finn3 = encephalitis_finn.loc[encephalitis_finn['Tissue'].isin(brain)]
encephalitis_brain = encephalitis.loc[encephalitis['tissue'].isin(brain)]

#Subset to Top BioVU Associations
genes = encephalitis_brain['gene']
encephalitis_finn2 = encephalitis_finn.loc[encephalitis_finn['gene_name'].isin(genes)]

#Write to csv files
encephalitis_finn2.to_csv(r'#Absolute path to output file location\Encephalitis_Brain_Finngen.csv', index=False, encoding='utf-8')


#Subsets FinnGen to Brain Region
encephalitis_regions_finn = []
for i in range(len(brain)):
    encephalitis_finn2 = encephalitis_finn.loc[encephalitis_finn['Tissue'] == brain[i]]
    encephalitis_regions_finn.append(encephalitis_finn2)

#Subsets BioVU to Brain Region
encephalitis_regions_biovu = []
for i in range(len(brain)):
    encephalitis_brain2 = encephalitis.loc[encephalitis['tissue'] == brain2[i]]
    encephalitis_regions_biovu.append(encephalitis_brain2)

#Loop
final_list =[]
for i in range(len(brain)):
    p = encephalitis_regions_finn[i]
    print(p)
    gene_list = []
    gns = p["gene_name"]
    filterer = encephalitis_regions_biovu[i]
    gnsFilter = filterer["gene"]
    out = p[p["gene_name"].isin(gnsFilter)]
    final_list.append(out)

#Region Specific to csv
final_list[0].to_csv(r'#Absolute path to output file location\Encephalitis_Anterior_Finngen.csv', index=False, encoding='utf-8')
final_list[1].to_csv(r'#Absolute path to output file location\Encephalitis_Caudate_Finngen.csv', index=False, encoding='utf-8')
final_list[2].to_csv(r'#Absolute path to output file location\Encephalitis_CerebellarHemi_Finngen.csv', index=False, encoding='utf-8')
final_list[3].to_csv(r'#Absolute path to output file location\Encephalitis_Cerebellum_Finngen.csv', index=False, encoding='utf-8')
final_list[4].to_csv(r'#Absolute path to output file location\Encephalitis_Cortex_Finngen.csv', index=False, encoding='utf-8')
final_list[5].to_csv(r'#Absolute path to output file location\Encephalitis_FrontalCortex_Finngen.csv', index=False, encoding='utf-8')
final_list[6].to_csv(r'#Absolute path to output file location\Encephalitis_Hipppocampus_Finngen.csv', index=False, encoding='utf-8')
final_list[7].to_csv(r'#Absolute path to output file location\Encephalitis_Hypothalamus_Finngen.csv', index=False, encoding='utf-8')
final_list[8].to_csv(r'#Absolute path to output file location\Encephalitis_NucleusAccumbens_Finngen.csv', index=False, encoding='utf-8')
final_list[9].to_csv(r'#Absolute path to output file location\Encephalitis_Putamen_Finngen.csv', index=False, encoding='utf-8')


#genes = encephalitis_brain['gene']
#encephalitis_finn2 = encephalitis_finn.loc[encephalitis_finn['gene_name'].isin(genes)]


#Fifth trait meningitis
#Import first FinnGen Trait of Interest Influenza
files2 = glob.glob(r'#Absolute path to all Meningitis Predixcan Results'+'/*.csv')
dfs2 = []
for file2 in files2:
    df2 = pd.read_csv(file2, index_col=False, header=0)
    df2['Tissue']=file2.split('\\')[7].rstrip('.csv')
    dfs2.append(df2)
print(dfs2)

#Subset to Brain tissue
meningitis_finn = pd.concat(dfs2)
brain = ['Brain_Anterior_cingulate_cortex_BA24','Brain_Caudate_basal_ganglia','Brain_Cerebellar_Hemisphere','Brain_Cerebellum','Brain_Cortex','Brain_Frontal_Cortex_BA9','Brain_Hippocampu','Brain_Hypothalamu','Brain_Nucleus_accumbens_basal_ganglia','Brain_Putamen_basal_ganglia']
brain2 = ['Brain_Anterior_cingulate_cortex_BA24','Brain_Caudate_basal_ganglia','Brain_Cerebellar_Hemisphere','Brain_Cerebellum','Brain_Cortex','Brain_Frontal_Cortex_BA9','Brain_Hippocampus','Brain_Hypothalamus','Brain_Nucleus_accumbens_basal_ganglia','Brain_Putamen_basal_ganglia']
meningitis_finn = meningitis_finn.loc[meningitis_finn['Tissue'].isin(brain)]
meningitis_brain = meningitis.loc[meningitis['tissue'].isin(brain2)]

#Subset to Top BioVU Associations
genes = meningitis_brain['gene']
meningitis_finn2 = meningitis_finn.loc[meningitis_finn['gene_name'].isin(genes)]

#Write to csv files
meningitis_finn2.to_csv(r'#Absolute path to output file location\Meningitis_Brain_Finngen.csv', index=False, encoding='utf-8')

#Subsets FinnGen to Brain Region
meningitis_regions_finn = []
for i in range(len(brain)):
    meningitis_finn2 = meningitis_finn.loc[meningitis_finn['Tissue'] == brain[i]]
    meningitis_regions_finn.append(meningitis_finn2)

#Subsets BioVU to Brain Region
meningitis_regions_biovu = []
for i in range(len(brain)):
    meningitis_brain2 = meningitis.loc[meningitis['tissue'] == brain2[i]]
    meningitis_regions_biovu.append(meningitis_brain2)

#Loop
final_list =[]
for i in range(len(brain)):
    p = meningitis_regions_finn[i]
    print(p)
    gene_list = []
    gns = p["gene_name"]
    filterer = meningitis_regions_biovu[i]
    gnsFilter = filterer["gene"]
    out = p[p["gene_name"].isin(gnsFilter)]
    final_list.append(out)

#Region Specific to csv
final_list[0].to_csv(r'#Absolute path to output file location\Meningitis_Anterior_Finngen.csv', index=False, encoding='utf-8')
final_list[1].to_csv(r'#Absolute path to output file location\Meningitis_Caudate_Finngen.csv', index=False, encoding='utf-8')
final_list[2].to_csv(r'#Absolute path to output file location\Meningitis_CerebellarHemi_Finngen.csv', index=False, encoding='utf-8')
final_list[3].to_csv(r'#Absolute path to output file location\Meningitis_Cerebellum_Finngen.csv', index=False, encoding='utf-8')
final_list[4].to_csv(r'#Absolute path to output file location\Meningitis_Cortex_Finngen.csv', index=False, encoding='utf-8')
final_list[5].to_csv(r'#Absolute path to output file location\Meningitis_FrontalCortex_Finngen.csv', index=False, encoding='utf-8')
final_list[6].to_csv(r'#Absolute path to output file location\Meningitis_Hipppocampus_Finngen.csv', index=False, encoding='utf-8')
final_list[7].to_csv(r'#Absolute path to output file location\Meningitis_Hypothalamus_Finngen.csv', index=False, encoding='utf-8')
final_list[8].to_csv(r'#Absolute path to output file location\Meningitis_NucleusAccumbens_Finngen.csv', index=False, encoding='utf-8')
final_list[9].to_csv(r'#Absolute path to output file location\Meningitis_Putamen_Finngen.csv', index=False, encoding='utf-8')


