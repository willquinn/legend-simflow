# legend-simflow

![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/legend-exp/legend-simflow?logo=git)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![GitHub issues](https://img.shields.io/github/issues/legend-exp/legend-simflow?logo=github)
![GitHub pull requests](https://img.shields.io/github/issues-pr/legend-exp/legend-simflow?logo=github)
![License](https://img.shields.io/github/license/legend-exp/legend-simflow)

End-to-end Snakemake workflow to run Monte Carlo simulations of signal and background
signatures in the LEGEND experiment and produce probability-density functions.
Configuration metadata (e.g. rules for generating simulation macros or
post-processing settings) is stored at
[legend-simflow-config](https://github.com/legend-exp/legend-simflow-config).

## Workflow steps

1. Macro files are generated and writted to disk according to rules
   defined in the metadata. These macros are in one-to-one correspondence
   with simulation jobs in the next tiers
1. Tier `ver` building: run simulations that generate Monte Carlo
   event vertices needed to (some) simulations in the next tier
1. Tier `raw` building: run full event simulations.
1. *To be documented...*

## Setup

1. Set up [legend-prodenv](https://github.com/legend-exp/legend-prodenv)
   to collect software dependencies, instantiate and manage multiple production
   environments. Snakemake can be installed by following instructions in
   [`legend-prodenv/README.md`](https://github.com/legend-exp/legend-prodenv/blob/main/README.md).

1. Instantiate a new production cycle:
   ```console
   > git clone git@github.com:legend-exp/legend-prodenv
   > cd legend-prodenv && source setup.sh
   > simprod-init-cycle <path-to-cycle-directory>
   ```

1. Customize `<path-to-cycle-directory>/config.json` (see next section).

## The configuration file

*To be documented...*

- A dedicated configuration should be defined in `setups` for each different experiment.
  The same name is used in the metadata. The experiment for which a production should be
  run is defined in the `Snakefile`.
- The `benchmark` section is used to configure a benchmarking run (see below)
- The `runcmd` section defines for each tier the command to be executed in order to
  produce output files. The `scripts/MaGe.sh` wrapper, which detects problems in the
  simulation output and exits with a proper code, should be used (and kept up-to-date)
  instead of the `MaGe` command. This allows Snakemake to better detect job failures.
- the `execenv` section defines the software environment (container) where all jobs
  should be executed (see below).

## Production

Run a production by using one of the provided profiles (recommended):

```console
> cd <path-to-cycle-directory>
> snakemake --profile workflow/profiles/<profile-name>
```

If no system-specific profiles are provided, the `default` profile should be used.

On a system providing [Scratch space](https://en.wikipedia.org/wiki/Scratch_space),
like NERSC, the `--shadow-prefix` option should be set to point to it:

```console
> snakemake --profile workflow/profiles/<profile-name> --shadow-prefix <path-to-scratch-area>
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

> **Note**
> 
> To make proper use of LEGEND's shifter containers, special permissions must be
> set on the production cycle directory (see
> [docs](https://docs.nersc.gov/development/shifter/faq-troubleshooting/#invalid-volume-map)):
> ```
> setfacl -R -m u:nobody:X <path-to-cycle-directory>
> ```

### Production

Start the production on the interactive node:
```
snakemake --profile workflow/profiles/nersc-interactive --shadow-prefix "$PSCRATCH"
```

Start the production on the batch nodes (via SLURM):
```
snakemake --profile workflow/profiles/nersc-batch --shadow-prefix "$PSCRATCH"
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

## Running jobs in the LEGEND software container

### with [container-env](https://github.com/oschulz/container-env)

In `config.json`:
```json
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
and the environment variables are needed for MaGe to work.

### with NERSC Shifter

In `config.json`:
```json
"execenv": [
    "shifter",
    "--image",
    "legendexp/legend-software:latest",
    "--volume $_/inputs/simprod/MaGe:/private",
    "--env MESHFILESPATH=/private/data/legendgeometry/stl_files",
    "--env MAGERESULTS=/private/data/legendgeometry"
]
```
