# legend-simflow

![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/legend-exp/legend-simflow?logo=git)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![GitHub issues](https://img.shields.io/github/issues/legend-exp/legend-simflow?logo=github)
![GitHub pull requests](https://img.shields.io/github/issues-pr/legend-exp/legend-simflow?logo=github)
![License](https://img.shields.io/github/license/legend-exp/legend-simflow)

## Rationale

The workflow metadata (e.g. rules for generating simulation macros,
post-processing metadata) is stored in
[legend-simflow-config](https://github.com/legend-exp/legend-simflow-config).

The simulated data is organized in tiers:

- `ver`: stands for "vertices"
- `raw`: contains the output of the Geant4 simulations

*To be documented...*

## Setup

[legend-prodenv](https://github.com/legend-exp/legend-prodenv) should be used
to collect software dependencies, instantiate and manage multiple production
environments. Snakemake can be installed by following instructions in
[`legend-prodenv/README.md`](https://github.com/legend-exp/legend-prodenv).

To instantiate a new production cycle:

```
git clone git@github.com:legend-exp/legend-prodenv
cd legend-prodenv && source setup.sh
simprod-init-cycle <path-to-cycle-directory>
```

Before proceeding with the production, `<path-to-cycle-directory>/config.json`
should be customized.

### Running jobs in the LEGEND container with [container-env](https://github.com/oschulz/container-env)

In `config.json`:
```
"execenv": [
    "MESHFILESPATH=$_/inputs/simprod/MaGe/data/legendgeometry/stl_files",
    "MAGERESULTS=$_/inputs/simprod/MaGe/data/legendgeometry",
    "cenv", "legendsw"
]
```
where the `legendsw` environment has been previously created with
```
cenv --create legendsw <path-to-container>
```

## Production

Run a production by using one of the provided profiles (recommended):
```
snakemake --profile workflow/profiles/<profile-name>
```
Find some useful Snakemake command-line options at the bottom of this page.

> **Warning**
>
> Geant4 macro files are marked as
> [`ancient`](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#ignoring-timestamps)
> as a workaround to the fact that they might have been re-generated (i.e. they have a more recent
> creation time) but with the same content as before (i.e. there is no need
> to re-run the simulations). Since Snakemake removes output files before re-generating,
> there is no simple way to avoid overwriting files with unchanged content.
> The workaround lets user increase simulation statistics by running new jobs (identical to
> those already on disk). However, **If the macro content is updated in the metadata, users will need
> to manually remove the output simulation files or force execution.**

### Benchmarking runs

This workflow implements the possibility to run special "benchmarking" runs in
order to evaluate the speed of simulations, for tuning the number of events to
simulate for each simulation run.

*To be documented...*

## NERSC-specific instructions

### Setup

As an alternative to installing Snakemake through legend-prodenv's tools,
[NERSC's Mamba can be
used](https://docs.nersc.gov/jobs/workflow/snakemake/#building-an-environment-containing-snakemake).

To make proper use of LEGEND's shifter containers, special permissions must be
set on the production cycle directory (see
[docs](https://docs.nersc.gov/development/shifter/faq-troubleshooting/#invalid-volume-map)):
```
setfacl -R -m u:nobody:X <path-to-cycle-directory>
```

### Production

Start the production on the interactive node:
```
snakemake --profile workflow/profiles/nersc-interactive
```

Start the production on the batch nodes:
```
snakemake --profile workflow/profiles/nersc-batch
```

## Useful Snakemake CLI options

```
usage: snakemake [OPTIONS] -- [TARGET ...]

  --dry-run, --dryrun, -n
                        Do not execute anything, and display what would be done. If you have a very large workflow,
                        use --dry-run --quiet to just print a summary of the DAG of jobs.
  --profile PROFILE     Name of profile to use for configuring Snakemake.
  --cores [N], -c [N]   Use at most N CPU cores/jobs in parallel. If N is omitted or 'all', the limit is set to the
                        number of available CPU cores.
  --config [KEY=VALUE ...], -C [KEY=VALUE ...]
                        Set or overwrite values in the workflow config object.
  --configfile FILE [FILE ...], --configfiles FILE [FILE ...]
                        Specify or overwrite the config file of the workflow.
  --touch, -t           Touch output files (mark them up to date without really changing them) instead of running
                        their commands.
  --keep-going, -k      Go on with independent jobs if a job fails.
  --forceall, -F        Force the execution of the selected (or the first) rule and all rules it is dependent on
                        regardless of already created output.
  --shadow-prefix DIR   Specify a directory in which the 'shadow' directory is created. If not supplied, the value
                        is set to the '.snakemake' directory relative to the working directory.
  --list, -l            Show available rules in given Snakefile. (default: False)
  --list-target-rules, --lt
                        Show available target rules in given Snakefile.
  --summary, -S         Print a summary of all files created by the workflow.
  --printshellcmds, -p  Print out the shell commands that will be executed. (default: False)
  --quiet [{progress,rules,all} ...], -q [{progress,rules,all} ...]
                        Do not output certain information. If used without arguments, do not output any progress or
                        rule information. Defining 'all' results in no information being printed at all.
```
