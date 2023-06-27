# legend-simflow

![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/legend-exp/legend-simflow?logo=git)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![GitHub issues](https://img.shields.io/github/issues/legend-exp/legend-simflow?logo=github)
![GitHub pull requests](https://img.shields.io/github/issues-pr/legend-exp/legend-simflow?logo=github)
![License](https://img.shields.io/github/license/legend-exp/legend-simflow)

### Usage

```
snakemake -j --configfile config.json
snakemake -j --configfile config.json -- print_stats
```

### Useful Snakemake CLI options

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
