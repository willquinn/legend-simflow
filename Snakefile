# TODO: simid lists
# TODO: snakemake wrapper scripts

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
setup.setdefault("benchmark", {"enabled": False, "n_primaries": 5000})


wildcard_constraints:
    tier="\w+",
    simid="[-\w]+",
    jobid="\w+",


rule gen_all_macros:
    """Aggregate and produce all the macro files."""
    input:
        aggregate.gen_list_of_all_macros(setup, tier="ver"),
        aggregate.gen_list_of_all_macros(setup, tier="raw"),


rule gen_all_tier_raw:
    """Aggregate and produce all the 'raw' tier files."""
    input:
        aggregate.gen_list_of_all_simid_outputs(setup, tier="ver"),
        aggregate.gen_list_of_all_simid_outputs(setup, tier="raw"),
    default_target: True


# since the number of generated macros for the 'output' field
# must be deduced at runtime from the JSON configuration, we need here to
# generate a separate rule for each 'simid'
simconfigs = aggregate.collect_simconfigs(setup, ["ver", "raw"])

for tier, simid, n_macros in simconfigs:

    rule:
        f"""Generates all needed simulation macros ({n_macros})
        for {simid} in tier '{tier}'. No wildcards are used.
        """
        localrule: True
        input:
            **patterns.macro_gen_inputs(setup, tier, simid),
        output:
            patterns.input_simjob_filenames(setup, n_macros, tier=tier, simid=simid),
        params:
            tier=tier,
            simid=simid,
            setup=setup,
        threads: 1
        message:
            f"Generating macros for '{tier}.{simid}'"
        script:
            "scripts/generate_macros.py"


rule tier_ver:
    """Run a single simulation job for the 'ver' tier.
    Uses wildcards `simid` and `jobid`.

    Warning
    -------
    The macro file is marked as "ancient" as a workaround to the fact that
    it might have been re-generated (i.e. it effectively has a more recent
    creation time) but with the same content as before (i.e. there is no need
    to re-run the simulation). If the macro content is updated, users will need
    to manually remove the output simulation files or force execution.
    """
    message:
        "Producing output file for job 'ver.{simid}.{jobid}'"
    input:
        ancient(patterns.input_simjob_filename(setup, tier="ver")),
    output:
        protected(patterns.output_simjob_filename(setup, tier="ver")),
    log:
        patterns.log_file_path(setup, "ver"),
    benchmark:
        patterns.benchmark_file_path(setup, "ver")
    shadow:
        "minimal"
    threads: 1
    shell:
        patterns.run_command(setup, "ver")


rule tier_raw:
    """Run a single simulation job for the 'raw' tier.
    Uses wildcards `simid` and `jobid`.

    Warning
    -------
    The macro file is marked as "ancient" as a workaround to the fact that
    it might have been re-generated (i.e. it effectively has a more recent
    creation time) but with the same content as before (i.e. there is no need
    to re-run the simulation). If the macro content is updated, users will need
    to manually remove the output simulation files or force execution.
    """
    message:
        "Producing output file for job 'raw.{simid}.{jobid}'"
    input:
        macro=ancient(patterns.input_simjob_filename(setup, tier="raw")),
        verfile=lambda wildcards: patterns.smk_ver_filename_for_raw(setup, wildcards),
    output:
        protected(patterns.output_simjob_filename(setup, tier="raw")),
    log:
        patterns.log_file_path(setup, "raw"),
    benchmark:
        patterns.benchmark_file_path(setup, "raw")
    shadow:
        "minimal"
    threads: 1
    shell:
        patterns.run_command(setup, "raw")


rule print_stats:
    """Prints a table with summary runtime information for each `simid`.
    No wildcards are used.
    """
    input:
        aggregate.gen_list_of_all_simid_outputs(setup, tier="ver"),
        aggregate.gen_list_of_all_simid_outputs(setup, tier="raw"),
    params:
        setup=setup,
    script:
        "scripts/print_stats.py"
