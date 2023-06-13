from pathlib import Path

from scripts.utils import utils, simjobs, patterns, aggregate

# NOTE: must set config file via --configfile

if not config:
    raise RuntimeError("you must set a config file with --configfile")

utils.subst_vars_in_snakemake_config(workflow, config)

experiment = "l200a"
setup = config["setups"][experiment]
setup["experiment"] = experiment
swenv = utils.runcmd(setup["execenv"])


rule all:
    """Run the entire pdf production."""
    input:
        # list here all output pdf files by parsing filelist
        aggregate.gen_list_of_all_simid_outputs(setup, tier="ver"),
        aggregate.gen_list_of_all_simid_outputs(setup, tier="raw"),


wildcard_constraints:
    tier="\w+",
    simid="[-\w]+",
    jobid="\w+",


# since the number of generated macros for the 'output' field
# must be deduced at runtime from the JSON configuration, we need here to
# generate a separate rule for each 'simid'
simconfigs = aggregate.collect_simconfigs(setup, ["ver", "raw"])

for tier, simid, n_macros in simconfigs:

    rule:
        f"""Generates all needed simulation macros ({n_macros})
        for {simid} in tier '{tier}'.
        """
        localrule: True
        input:
            **patterns.macro_gen_inputs(setup, tier=tier, simid=simid),
        output:
            patterns.input_simjob_filenames(setup, n_macros, tier=tier, simid=simid),
        params:
            tier=tier,
            simid=simid,
            setup=setup,
        threads: 1
        script:
            "scripts/generate_macros.py"


rule:
    """Run a single simulation job for the 'ver' tier.
    Uses wildcards `simid` and `jobid`.
    """
    input:
        patterns.input_simjob_filename(setup, tier="ver"),
    output:
        patterns.output_simjob_filename(setup, tier="ver"),
    log:
        patterns.log_file_path(setup, "ver"),
    threads: 1
    shell:
        patterns.run_command(setup, "ver")

rule:
    """Run a single simulation job for the 'raw' tier.
    Uses wildcards `simid` and `jobid`.
    """
    input:
        macro=patterns.input_simjob_filename(setup, tier="raw"),
        verfile=lambda wildcards: patterns.smk_ver_filename_for_raw(setup, wildcards),
    output:
        patterns.output_simjob_filename(setup, tier="raw"),
    log:
        patterns.log_file_path(setup, "raw"),
    threads: 1
    shell:
        patterns.run_command(setup, "raw")
