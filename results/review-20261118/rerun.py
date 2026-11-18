import hashlib
from pathlib import Path
import pandas as pd
from mitoage_pipelines import paths
from mitoage_pipelines.longevity import rates, structure_pgls
from mitoage_pipelines.adaptations import longevity
from mitoage_pipelines.mitoport import cybrid, cybrid_energy
r=paths.ROOT
print('DATA ROOT',r,flush=True)
print('SPECIES SHA256',hashlib.sha256((paths.BUILD/'species.parquet').read_bytes()).hexdigest(),flush=True)
print('H7 rates',flush=True)
rates.main(paths.BUILD,paths.LONGEVITY_OUT/'composition.tsv',paths.LONGEVITY_OUT/'fourfold.tsv',paths.DATED_TREE,paths.LONGEVITY_OUT)
print('H12 structural regressions',flush=True)
structure_pgls.main(paths.BUILD,paths.LONGEVITY_OUT/'composition.tsv',paths.LONGEVITY_OUT/'fourfold.tsv',paths.DATED_TREE,paths.MITOPORT_IFACE,paths.LONGEVITY_OUT)
print('H26 traits and pairs',flush=True)
longevity.main_traits(paths.BUILD,paths.DATED_TREE,paths.ADAPT_LONGEVITY_OUT)
longevity.run_pairs(paths.BUILD,paths.DATED_TREE,paths.ADAPT_LONGEVITY_OUT)
print('H16 and H17 readouts',flush=True)
cal,summary=cybrid.summarise(pd.read_csv(paths.CYBRID_OUT/'scores.tsv',sep='\t'),cybrid.read_measurements())
cal.to_csv(paths.CYBRID_OUT/'calibration.tsv',sep='\t',index=False)
summary.to_csv(paths.CYBRID_OUT/'summary.tsv',sep='\t',index=False)
cybrid.readout(paths.CYBRID_OUT)
print('H18 partial collection: unavailable tasks must remain unassessed',flush=True)
cybrid_energy.collect(paths.CYBRID_OUT)
cybrid.readout_energy(paths.CYBRID_OUT)
print('ALL AVAILABLE RERUNS COMPLETED',flush=True)
