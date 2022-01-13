# About

This repository analyzes viral genomes using [Nextstrain](https://nextstrain.org) to understand how SARS-CoV-2, the virus that is responsible for the COVID-19 pandemic, evolves and spreads. This workflow is copied from [github.com/nextstrain/ncov](https://github.com/nextstrain/ncov).

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