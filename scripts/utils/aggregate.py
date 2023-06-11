import json

from . import simjobs


def collect_simconfigs(setup, tiers):
    cfgs = []
    for tier in tiers:
        with (simjobs.template_macro_dir(setup, tier) / "simconfig.json").open() as f:
            for sid, val in json.load(f).items():
                cfgs.append((tier, sid, simjobs.get_simid_n_macros(setup, tier, sid)))

    return cfgs


def gen_list_of_all_simids(setup, tier):
    with (simjobs.template_macro_dir(setup, tier) / "simconfig.json").open() as f:
        return json.load(f).keys()


def gen_list_of_all_macros(setup, tier):
    mlist = []
    for sid in gen_list_of_all_simids(setup, tier):
        mlist += simjobs.gen_list_of_simid_inputs(setup, tier, sid)

    return mlist


def gen_list_of_all_simid_outputs(setup, tier):
    mlist = []
    for sid in gen_list_of_all_simids(setup, tier):
        mlist += simjobs.gen_list_of_simid_outputs(setup, tier, sid)

    return mlist
