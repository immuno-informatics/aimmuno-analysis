"""."""

import random
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

# Root TransPHLA-AOMP directory
transphla_dir = Path("/path/to/TransPHLA-AOMP")

datasets_dir = transphla_dir.joinpath("Dataset")

output_dir = datasets_dir

ext_file = datasets_dir.joinpath("external_set.csv")
indep_file = datasets_dir.joinpath("independent_set.csv")


def train_file(fold):
    """File name."""
    return datasets_dir.joinpath(f"train_data_fold{fold}.csv")


def val_file(fold):
    """File name."""
    return datasets_dir.joinpath(f"val_data_fold{fold}.csv")


def ext_rand_file(suffix=""):
    """File name."""
    return output_dir.joinpath(f"external_set_rand{suffix}.csv")


def indep_rand_file(suffix=""):
    """File name."""
    return output_dir.joinpath(f"independent_set_rand{suffix}.csv")


def train_rand_file(fold, suffix=""):
    """File name."""
    return output_dir.joinpath(f"train_data_fold{fold}_rand{suffix}.csv")


def val_rand_file(fold, suffix=""):
    """File name."""
    return output_dir.joinpath(f"val_data_fold{fold}_rand{suffix}.csv")


data_file_col_pep = "peptide"  # column name
data_file_col_hla_seq = "HLA_sequence"  # column name
data_file_col_pep_len = "length"  # column name

# Setting random seed for reproducibility
seed = 19961231
random.seed(seed)
np.random.seed(seed)


def load_src_data(fold):
    """."""
    ext = pd.read_csv(ext_file, index_col=0)
    indep = pd.read_csv(indep_file, index_col=0)
    train = pd.read_csv(train_file(fold), index_col=0)
    val = pd.read_csv(val_file(fold), index_col=0)
    return train, val, indep, ext


def randomize_hlas(data, all_data_list):
    """."""
    all_data = (
        pd.concat(all_data_list)
        .loc[:, [data_file_col_pep, data_file_col_hla_seq]]
        .drop_duplicates()
    )

    u_hlas = all_data[data_file_col_hla_seq].unique()

    all_data = all_data.groupby(data_file_col_pep)[data_file_col_hla_seq].apply(list)

    new_hlas = []
    for pep in tqdm(data[data_file_col_pep]):
        poss_hlas = all_data[pep]
        hla_choice = u_hlas[~np.isin(u_hlas, poss_hlas)]
        new_hlas.append(np.random.choice(hla_choice))

    data[data_file_col_hla_seq] = new_hlas

    return data


def randomize_peptides(data, all_data_list):
    """."""
    all_data = (
        pd.concat(all_data_list)
        .loc[:, [data_file_col_pep, data_file_col_hla_seq]]
        .drop_duplicates()
    )

    u_peps = set(all_data[data_file_col_pep])

    all_data = all_data.groupby(data_file_col_hla_seq)[data_file_col_pep].apply(set)

    data_hlas = data.groupby(data_file_col_hla_seq)[data_file_col_pep].apply(list)
    data_hlas_n = data_hlas.str.len()

    for hla, p_n in zip(tqdm(data_hlas.index), data_hlas_n):
        poss_peps = all_data[hla]
        pep_choice = sorted(u_peps - poss_peps)
        sel_peps = random.sample(pep_choice, k=p_n)
        data.loc[data[data_file_col_hla_seq] == hla, data_file_col_pep] = sel_peps

    data[data_file_col_pep_len] = data[data_file_col_pep].str.len().to_list()

    return data


# Process only the first fold (as originally in the notebook)
fold = 0

# HLA sequences randomization

# Randomize every data file

# the output files have a `_rand` suffix already attached to the end
suffix = "_hla_all"

print(f"Processing: {suffix}")

train, val, indep, ext = load_src_data(fold)

all_data = [train, val, indep, ext]
train = randomize_hlas(train, all_data)
train.to_csv(train_rand_file(fold, suffix=suffix))

all_data = [train, val, indep, ext]
val = randomize_hlas(val, all_data)
val.to_csv(val_rand_file(fold, suffix=suffix))

all_data = [train, val, indep, ext]
indep = randomize_hlas(indep, all_data)
indep.to_csv(indep_rand_file(suffix=suffix))

all_data = [train, val, indep, ext]
ext = randomize_hlas(ext, all_data)
ext.to_csv(ext_rand_file(suffix=suffix))

# Randomize only the independent and external data

# the output files have a `_rand` suffix already attached to the end
suffix = "_hla_test_only"

print(f"Processing: {suffix}")

train, val, indep, ext = load_src_data(fold)

train.to_csv(train_rand_file(fold, suffix=suffix))

val.to_csv(val_rand_file(fold, suffix=suffix))

all_data = [train, val, indep, ext]
indep = randomize_hlas(indep, all_data)
indep.to_csv(indep_rand_file(suffix=suffix))

all_data = [train, val, indep, ext]
ext = randomize_hlas(ext, all_data)
ext.to_csv(ext_rand_file(suffix=suffix))

# Randomize only the training and validation data

# the output files have a `_rand` suffix already attached to the end
suffix = "_hla_train_only"

print(f"Processing: {suffix}")

train, val, indep, ext = load_src_data(fold)

all_data = [train, val, indep, ext]
train = randomize_hlas(train, all_data)
train.to_csv(train_rand_file(fold, suffix=suffix))

all_data = [train, val, indep, ext]
val = randomize_hlas(val, all_data)
val.to_csv(val_rand_file(fold, suffix=suffix))

indep.to_csv(indep_rand_file(suffix=suffix))

ext.to_csv(ext_rand_file(suffix=suffix))

# Peptide sequences randomization

# Randomize every data file

# the output files have a `_rand` suffix already attached to the end
suffix = "_pep_all"

print(f"Processing: {suffix}")

train, val, indep, ext = load_src_data(fold)

all_data = [train, val, indep, ext]
train = randomize_peptides(train, all_data)
train.to_csv(train_rand_file(fold, suffix=suffix))

all_data = [train, val, indep, ext]
val = randomize_peptides(val, all_data)
val.to_csv(val_rand_file(fold, suffix=suffix))

all_data = [train, val, indep, ext]
indep = randomize_peptides(indep, all_data)
indep.to_csv(indep_rand_file(suffix=suffix))

all_data = [train, val, indep, ext]
ext = randomize_peptides(ext, all_data)
ext.to_csv(ext_rand_file(suffix=suffix))

# Randomize only the independent and external data

# the output files have a `_rand` suffix already attached to the end
suffix = "_pep_test_only"

print(f"Processing: {suffix}")

train, val, indep, ext = load_src_data(fold)

train.to_csv(train_rand_file(fold, suffix=suffix))

val.to_csv(val_rand_file(fold, suffix=suffix))

all_data = [train, val, indep, ext]
indep = randomize_peptides(indep, all_data)
indep.to_csv(indep_rand_file(suffix=suffix))

all_data = [train, val, indep, ext]
ext = randomize_peptides(ext, all_data)
ext.to_csv(ext_rand_file(suffix=suffix))

# Randomize only the training and validation data

# the output files have a `_rand` suffix already attached to the end
suffix = "_pep_train_only"

print(f"Processing: {suffix}")

train, val, indep, ext = load_src_data(fold)

all_data = [train, val, indep, ext]
train = randomize_peptides(train, all_data)
train.to_csv(train_rand_file(fold, suffix=suffix))

all_data = [train, val, indep, ext]
val = randomize_peptides(val, all_data)
val.to_csv(val_rand_file(fold, suffix=suffix))

indep.to_csv(indep_rand_file(suffix=suffix))

ext.to_csv(ext_rand_file(suffix=suffix))
