#Remove any rows from the FinnGen trait gwas that does not have an rsid

import pandas as pd

filename = '#Absolute path to FinnGen gwas trait file'

df = pd.read_csv(filename, compression='gzip', sep='\t')
rsids = df['rsids']


rs_names = rsids[rsids.str.contains('rs', na=False)]
df3 = df[df['rsids'].isin(rs_names)]


df3.to_csv('#Output file name .txt', sep='\t', index=False)