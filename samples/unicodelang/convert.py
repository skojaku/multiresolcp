import pandas as pd
import numpy as np

country_code = pd.read_csv('ent.unicodelang.coutry.code', sep='\t', header=None)
lang_code = pd.read_csv('ent.unicodelang.language.code', sep='\t', header=None)


df = pd.read_csv('out.unicodelang', sep='\t', header=None)

c = country_code.values[df[0].values-1]
l = lang_code.values[df[1].values-1]
df = pd.DataFrame({'source':pd.Series(c.flatten()), 'target':pd.Series(l.flatten())})

df.to_csv('sample.dat', '\t', index=False)

print(df)
