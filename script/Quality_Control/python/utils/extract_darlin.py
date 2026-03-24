import pandas as pd
import numpy as np
import os

def extract(adata, output, whitelist_path, darlin):
    whitelist = [line.strip() for line in open(whitelist_path).readlines()]
    output_path = output + '/' + darlin
    os.makedirs(output_path, exist_ok=True)
    adata_sub = adata[:, f'DARLIN-{darlin}']
    umi = adata_sub.X.toarray().ravel()
    df = pd.DataFrame({
    'barcode': adata_sub.obs_names,
    'umi_count': umi
    })
    df = df[df['umi_count'] != 0]
    df['xbc'] = df['barcode'].str[8:16]
    df['ybc'] = df['barcode'].str[:8]
    df['x'] = df['xbc'].apply(lambda bc: whitelist.index(bc) if bc in whitelist else -1)
    df['y'] = df['ybc'].apply(lambda bc: whitelist.index(bc) if bc in whitelist else -1)
    df.to_csv(output_path + '/' + darlin + '.csv', index=False)

    return output_path + '/' + darlin + '.csv'