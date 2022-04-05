# About

This repository analyzes viral genomes using [Nextstrain](https://nextstrain.org) to understand how SARS-CoV-2, the virus that is responsible for the COVID-19 pandemic, evolves and spreads. 
This workflow was copied from [github.com/nextstrain/ncov](https://github.com/nextstrain/ncov) on April 4, 2022.

This is an updated version of the workflow used in the manuscript ["Rapid and parallel adaptive mutations in spike S1 drive clade success in SARS-CoV-2"](https://www.cell.com/cell-host-microbe/pdf/S1931-3128(22)00148-2.pdf).
The updates improve quality control of the sequences included in the tree. 
The version of the workflow used in the manuscript is tagged here https://github.com/blab/sarscov2-adaptive-evolution/releases/tag/chm-publication.

# Running

After installing Nextstrain dependencies this workflow can be run as:
```
snakemake -j 1 -p --profile adaptive_evolution_profiles/adaptive-evolution
```

However, running via cluster or AWS will be faster:

Run a normal build (~10000 sequences) with:
```
nextstrain build --aws-batch --cpus 16 --memory 64GiB --detach . --set-threads tree=16 --profile adaptive_evolution_profiles --configfile adaptive_evolution_profiles/builds.yaml
```

Or a high density build (~20000 sequences) with:
```
nextstrain build --aws-batch --cpus 16 --memory 64GiB --detach . --set-threads tree=16 --profile adaptive_evolution_profiles/ --configfile adaptive_evolution_profiles/high_density_builds.yaml
```

Or a region-specific build (with ~10000 sequences)
```
nextstrain build --aws-batch --cpus 16 --memory 64GiB --detach . --set-threads tree=16 --profile adaptive_evolution_profiles --configfile adaptive_evolution_profiles/regional_builds.yaml
```