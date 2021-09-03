from vienna_interface import parse_structured_score, parse_score
from rRNA_23S import seq_23S, struct_23S
from rRNA_16S import seq_16S, struct_16S
import numpy as np
import json
import os
import seaborn as sns
from tqdm import tqdm
import pandas as pd

# ensure uppercase
seq_16S = seq_16S.upper()
seq_23S = seq_23S.upper()

# Look, this is obviously a terrible, bizarre way to coerce a json-serialized
# dictionary into a dataframe. I cannot honestly imagine why I did this a
# couple years ago. It does work, though.
def load_cache(jsonfile):
    with open(jsonfile) as f:
        x,y = [], []
        cache = json.loads(f.read())
        for k,v in cache.items():
            y.append(v)
            x.append(k)
        return x, y

x,y = load_cache('cache_bp_23S.json')
df_23S_mc = pd.DataFrame.from_dict([{'seq': k, 'delta': v} for k,v in zip(x,y)])
x,y = load_cache('cache_bp_16S.json')
df_16S_mc = pd.DataFrame.from_dict([{'seq': k, 'delta': v} for k,v in zip(x,y)])


def count_muts(seq, wt):
    return len([0 for c, c_ in zip(seq, wt) if c != c_])

df_16S_mc['muts'] = df_16S_mc.apply(lambda x: count_muts(x.seq.upper(), seq_16S.upper()), axis=1)
df_23S_mc['muts'] = df_23S_mc.apply(lambda x: count_muts(x.seq.upper(), seq_23S.upper()), axis=1)

df = pd.read_csv('pilot_sequences.csv')
df['muts'] = df.apply(lambda x: count_muts(x.Sequence.upper(), seq_23S.upper()) if x.Subunit == '23S' else count_muts(x.Sequence, seq_16S.upper()), axis=1)

# produce sequences with these mutation numbers
print('delta;seq;muts')
for mut_num in df[(df.Round == 'Pilot') & (df.Subunit == '23S')].muts.unique():
    num = len(df[(df.Round == 'Pilot') & (df.Subunit == '23S') & (df.muts == mut_num)])
    # emit one sequence per sequence from the pilot round with that number mutations
    for ii in range(num):
        try: 
            print('{};{};{}'.format(
                df_23S_mc[df_23S_mc.muts == mut_num].nsmallest(num, 'delta').delta.values[ii],
                df_23S_mc[df_23S_mc.muts == mut_num].nsmallest(num, 'delta').seq.values[ii],
                df_23S_mc[df_23S_mc.muts == mut_num].nsmallest(num, 'delta').muts.values[ii]
            ))
        except:
            print("stumped by mut_num", mut_num)

for mut_num in df[(df.Round == 'Pilot') & (df.Subunit == '16S')].muts.unique():
    num = len(df[(df.Round == 'Pilot') & (df.Subunit == '16S') & (df.muts == mut_num)])
    for ii in range(num):
        print('{};{};{}'.format(
            df_16S_mc[df_16S_mc.muts == mut_num].nsmallest(num, 'delta').delta.values[ii],
            df_16S_mc[df_16S_mc.muts == mut_num].nsmallest(num, 'delta').seq.values[ii],
            df_16S_mc[df_16S_mc.muts == mut_num].nsmallest(num, 'delta').muts.values[ii]
        ))
