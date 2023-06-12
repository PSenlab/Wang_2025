
#Installing mamba
cd /data/$USER
export TMPDIR=/lscratch/$SLURM_JOB_ID
echo "create fresh mamba env at /data directory"
wget https://gcc02.safelinks.protection.outlook.com/?url=https%3A%2F%2Fgithub.com%2Fconda-forge%2Fminiforge%2Freleases%2Flatest%2Fdownload%2FMambaforge-Linux-x86_64.sh&amp;data=05%7C01%7Clin.wang2%40nih.gov%7C3ede61309d4742fc693508dab60025aa%7C14b77578977342d58507251ca2dc2b06%7C0%7C0%7C638022408988084650%7CUnknown%7CTWFpbGZsb3d8eyJWIjoiMC4wLjAwMDAiLCJQIjoiV2luMzIiLCJBTiI6Ik1haWwiLCJXVCI6Mn0%3D%7C3000%7C%7C%7C&amp;sdata=rNXPJ53yZLBufSkeony%2FjRQ7PV9eJpDVQdMPoiJJGTc%3D&amp;reserved=0
bash Mambaforge-Linux-x86_64.sh -p /data/$USER/conda -b rm Mambaforge-Linux-x86_64.sh source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh
export MAMBA_NO_BANNER=1
mamba activate base
mamba update --all
mamba clean --all --yes

#Running Spatial DE
(base)mamba create -n spatialde
(base)mamba activate spatialde
(spatialde)pip install spatialde
(spatialde)pip install naivede
(spatialde)pip install matplotlib
(spatialde)python
   import matplotlib.pyplot as plt
   import pandas as pd
   import NaiveDE
   import SpatialDE
>>>counts = pd.read_csv('/path/to/sample_count.csv', index_col=0)
>>>counts = counts.T[counts.sum(0) >= 3].T
>>>sample_info = pd.read_csv('/path/to/sample_info.csv', index_col=0)
>>>counts = counts.loc[sample_info.index]
>>>sample_info.head(5)
>>>norm_expr = NaiveDE.stabilize(counts.T).T
>>>resid_expr = NaiveDE.regress_out(sample_info, norm_expr.T, 'np.log(total_counts)').T
>>>sample_resid_expr = resid_expr.sample(n=10109, axis=1, random_state=1,replace=True)
>>>X = sample_info[['X', 'Y']]
>>>results = SpatialDE.run(X, sample_resid_expr)
>>>results.head().T
>>>results.to_csv('/path/to/sample.csv', index=False)
>>>sign_results = results.query('qval < 0.05')
>>>sign_results['l'].value_counts()
>>>histology_results, patterns = SpatialDE.aeh.spatial_patterns(X, resid_expr, sign_results, C=4, l=520, verbosity=1)
>>>histology_results.head()
>>>for i in range(4):
    plt.subplot(1, 4, i + 1)
    plt.scatter(sample_info['X'], sample_info['Y'], c=patterns[i],s=3);
    plt.axis('equal')
    plt.title('Pattern {} - {} genes'.format(i, histology_results.query('pattern == @i').shape[0] ))
    plt.colorbar(ticks=[]);


>>>import os
>>>if not os.path.exists('/path/to/results'):
       os.makedirs('/results')
>>>with open('results/result.txt', 'w') as file:
       file.write(result)



#step 1#
$ module load python
$ pip install spatialde
$ python
pip install naivede
pip install matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import NaiveDE
import SpatialDE

#Spatially visualize select DE genes
figsize(5, 5)
for i, g in enumerate(["Cabp7","Opalin","Gpr88","Lamp5","Hcrt","Ccn3","Prkcd"]):
    plt.subplot(2, 4, i + 1)
    plt.scatter(sample_info['Y'], sample_info['X'], c=norm_expr[g],s=5);
    plt.title(g)
    plt.axis('equal')
    
plt.colorbar(ticks=[1, 4, 7])
plt.show()

fig, axs = plt.subplots(1, 1, figsize=(6, 6))
for i, g in enumerate(["Cabp7"]):
    plt.subplot(1, 1, i + 1)
    plt.scatter(sample_info['X'], sample_info['Y'], c=norm_expr[g], cmap=cmap,s=20);
    plt.title(g)
    plt.axis([0,2000,0,2000])
    plt.gca().set_aspect('equal', 'box')

plt.show()
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

fig, axs = plt.subplots(1, 1, figsize=(6, 6))
colors = ['silver','firebrick','red']
n_bins = [100]
cmap_name = 'my_list'
for n_bin, ax in zip(n_bins, axs.flat):
    cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)

    
plt.colorbar(ticks=[1, 2.5, 4])
plt.show()

#Spatial visualization analysis of select genes
for i, g in enumerate(["Ccn3"]):
    plt.subplot(1, 1, i + 1)
    plt.scatter(sample_info['X'], sample_info['Y'], c=norm_expr[g],s=6);
    plt.title(g)
    plt.axis('equal')

for i, g in enumerate(["Prkcd"]):
    plt.subplot(1, 1, i + 1)
    plt.scatter(sample_info['X'], sample_info['Y'], c=norm_expr[g],s=8);
    plt.title(g)
    plt.axis('equal')


for i, g in enumerate(["Opalin"]):
    plt.subplot(1, 1, i + 1)
    plt.scatter(sample_info['X'], sample_info['Y'], c=norm_expr[g],s=8);
    plt.title(g)
    plt.axis('equal')
    
plt.colorbar(ticks=[1, 1.75, 2.5])
plt.show()

for i, g in enumerate(["Gpr88"]):
    plt.subplot(1, 1, i + 1)
    plt.scatter(sample_info['X'], sample_info['Y'], c=norm_expr[g],s=8);
    plt.title(g)
    plt.axis('equal')
    
plt.colorbar(ticks=[1, 2, 3])
plt.show()

for i, g in enumerate(["Lamp5"]):
    plt.subplot(1, 1, i + 1)
    plt.scatter(sample_info['X'], sample_info['Y'], c=norm_expr[g],s=4);
    plt.title(g)
    plt.axis('equal')

plt.colorbar(ticks=[1, 2, 3])
plt.show()

for i, g in enumerate(["Hcrt"]):
    plt.subplot(1, 1, i + 1)
    plt.scatter(sample_info['X'], sample_info['Y'], c=norm_expr[g],s=8);
    plt.title(g)
    plt.axis('equal')
