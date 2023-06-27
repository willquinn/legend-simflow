from __future__ import annotations

import copy
import os
import string
from pathlib import Path


def subst_vars_impl(x, var_values, ignore_missing=False):
    if isinstance(x, str):
        if "$" in x:
            if ignore_missing:
                return string.Template(x).safe_substitute(var_values)
            else:
                return string.Template(x).substitute(var_values)
        else:
            return x
    if isinstance(x, dict):
        for key in x:
            value = x[key]
            new_value = subst_vars_impl(value, var_values, ignore_missing)
            if new_value is not value:
                x[key] = new_value
        return x
    if isinstance(x, list):
        for i in range(len(x)):
            value = x[i]
            new_value = subst_vars_impl(value, var_values, ignore_missing)
            if new_value is not value:
                x[i] = new_value
        return x
    else:
        return x


def subst_vars(props, var_values={}, use_env=False, ignore_missing=False):
    if use_env:
        combined_var_values = dict(iter(os.environ.items()))
        combined_var_values.update(copy.copy(var_values))
    else:
        combined_var_values = var_values

    return subst_vars_impl(props, combined_var_values, ignore_missing)


def subst_vars_in_snakemake_config(workflow, config):
    # TODO: better way of handling this?
    config_filename = workflow.overwrite_configfiles[0]
    return subst_vars(
        config,
        var_values={"_": str(Path(config_filename).parent)},
        use_env=True,
        ignore_missing=False,
    )
