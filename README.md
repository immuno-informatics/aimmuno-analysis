# AImmuno: A pipeline for testing various MHC class I antigen binding prediction models

[![AImmuno analysis code](https://zenodo.org/badge/DOI/10.5281/zenodo.14536592.svg)](https://doi.org/10.5281/zenodo.14536592) [![AImmuno analysis data](https://zenodo.org/badge/DOI/10.5281/zenodo.14526105.svg)](https://doi.org/10.5281/zenodo.14526105)

> Analysis and data generation pipeline for the "[Artificial Intelligence Perceives Marginal Gains from MHC Haplotype Data in Antigen Presentation Predictions](https://www.google.com)" paper (in peer review).

Supporting data and results can be found in the [paper data repository](https://doi.org/10.5281/zenodo.14526105).

## Setup

Some of the pipeline files and file processing descriptions assume work in a Linux environment.

### Download The Code

Please follow these steps to obtain the code to run the pipeline:

- Clone this repository and `cd` into it.

  ```bash
  git clone https://github.com/immuno-informatics/aimmuno-analysis.git
  cd aimmuno-analysis
  ```

  or

- [Download](https://github.com/immuno-informatics/aimmuno-analysis/archive/refs/heads/main.zip) contents of this repository, unzip it, and `cd` into it.

  ```bash
  wget https://github.com/immuno-informatics/aimmuno-analysis/archive/refs/heads/main.zip
  unzip main.zip
  rm main.zip
  cd aimmuno-analysis-main
  ```

### Python Environment

**The pipeline requires a working [Python](https://www.python.org) installation**. We suggest to install [Conda](https://docs.anaconda.com/miniconda) to handle Python operations.

A Conda environment file (`environment.yml`) was prepared to make it more convenient to set-up the analysis. Use it to create a new Python environment before running any of the scripts/notebooks described below (you can change the name `aimmuno` to anything else):

```bash
conda env create --name aimmuno --file environment.yml --yes
conda activate aimmuno
```

### Required Databases

The analysis pipeline requires certain external data files to produce all results. Please follow the subsections below for details.

#### CARMEN Database Files

The analysis pipeline runs with certain core CARMEN database files available at [the main repository](https://doi.org/10.5281/zenodo.13928441). Please download the following file:

- `carmen-hla-sequences.parquet`&mdash;complete amino acid and pseudo sequences of HLA molecules referenced in [the CARMEN database](https://doi.org/10.5281/zenodo.13928441).

And subsequently move it to the `data` directory.

#### Initial Data and Results Related to the Paper

To recreate all the results and models from scratch, we need to download a set of files that were used as input for specific stages of the analysis. Please download the contents of the [paper data repository](https://doi.org/10.5281/zenodo.14526105) and unpack them to the `data` directory.

### External Code and Data

Part of the analysis utilizes the [TransPHLA-AOMP](https://github.com/a96123155/TransPHLA-AOMP) model and its data created by [Chu et al.](https://www.nature.com/articles/s42256-022-00459-7) Please follow these steps to set up the base code:

1. Clone the repository:

    ```bash
    git clone https://github.com/a96123155/TransPHLA-AOMP.git
    ```

    or

    [Download](https://github.com/a96123155/TransPHLA-AOMP/archive/refs/heads/master.zip) contents of the repository and unzip it:

    ```bash
    wget https://github.com/a96123155/TransPHLA-AOMP/archive/refs/heads/master.zip
    unzip master.zip
    rm master.zip
    ```

2. Unzip `TransPHLA-AOMP/Dataset/external_set.zip` and `TransPHLA-AOMP/Dataset/independent_set.zip`:

    ```bash
    unzip TransPHLA-AOMP/Dataset/external_set.zip -d TransPHLA-AOMP/Dataset
    unzip TransPHLA-AOMP/Dataset/independent_set.zip -d TransPHLA-AOMP/Dataset
    ```

## Running the Analysis

To be able to reproduce the entire analysis, we need to set-up and run specific scripts in a particular order. Please follow the instructions below.

**Always remember to review all paths in the scripts and make changes appropriate to your own system's environment.**

**Be advised that running all steps will replace some files and may lead to a slight change of results due to stochasticity in certain parts of the analysis.**

### TAPE Model

#### CARMEN Data

##### 1. Prepare Training and Testing Datasets From CARMEN Data

Open the `scripts/1_carmen_data_preparation.ipynb` Jupyter notebook and follow its structure and instructions. By changing which set of lines is left uncommented in cell no. 3 we can select the type of experimental dataset to be produced.

**This step will produce completely randomized datasets and alter the results for the machine learning models. If you want to exactly reproduce the paper's results, skip this step and use the datasets provided in [the data repository](https://doi.org/10.5281/zenodo.14526105) to train and test the models.**

##### 2. Create and Test Models Without Using HLA Sequences

Open the `scripts/2_training_and_testing_without_hla_carmen_data.ipynb` Jupyter notebook and follow its structure and instructions. By changing which line is left uncommented in cell no. 3 we can select the type of model to be produced.

##### 3. Create and Test Models Using HLA Sequences

Open the `scripts/3_training_and_testing_with_hla_carmen_data.ipynb` Jupyter notebook and follow its structure and instructions. By changing which line is left uncommented in cell no. 3 we can select the type of model to be produced.

#### Synthetic Data for TAPE

##### 4. Create Synthetic Datasets

Here, we create two datasets with randomly-generated sequences and modified pairs of amino acids: the "easy" dataset with matching pairs, and the "hard" dataset with mismatches. Use the `scripts/4_synthetic_data_preparation.py` script with the following parameters to produce them:

```bash
cd scripts
python 4_synthetic_data_preparation.py --sample_num 100_000 --decoys 1 --peptide_min 9 --peptide_max 12 --mhc_min 20 --mhc_max 40 --out_dir "../data/datasets/synthetic_easy"
python 4_synthetic_data_preparation.py --sample_num 100_000 --decoys 1 --peptide_min 9 --peptide_max 12 --mhc_min 20 --mhc_max 40 --out_dir "../data/datasets/synthetic_hard" --hard_case
cd ..
```

**This step will produce completely randomized datasets and alter the results for the machine learning models. If you want to exactly reproduce the paper's results, skip this step and use the datasets provided in [the data repository](https://doi.org/10.5281/zenodo.14526105) to train and test the models.**

##### 5. Create and Test Models Without Using HLA Sequences

Open the `scripts/5_training_and_testing_without_hla_synthetic_data.ipynb` Jupyter notebook and follow its structure and instructions. By changing which line is left uncommented in cell no. 3 we can select the type of model to be produced.

##### 6. Create and Test Models Using HLA Sequences

Open the `scripts/6_training_and_testing_with_hla_synthetic_data.ipynb` Jupyter notebook and follow its structure and instructions. By changing which line is left uncommented in cell no. 3 we can select the type of model to be produced.

##### 7. Test the Models and See the Influence of Using HLA Sequences

This is an optional convenience script to easily test and visually compare TAPE-based models that use HLA sequences with those without such input. Because of stochasticity, this script may produce slightly different results that those saved during training scripts.

To run this step, open the `scripts/7_testing_models.ipynb` Jupyter notebook and follow its structure and instructions. By changing which line is left uncommented in cell no. 3 we can select the type of model to be tested.

### TransPHLA Model

Following subsequent steps require making adjustments to a source [TransPHLA-AOMP](https://github.com/a96123155/TransPHLA-AOMP) script to avoid publishing edited versions of the original code. Please, make two copies of the `/path/to/TransPHLA-AOMP/Procedure Code/pHLAIformer.ipynb` Jupyter notebook and name them `pHLAIformer_src.ipynb` and `pHLAIformer_synthetic.ipynb` for clarity.

#### Source Data

##### 8. Prepare Additional Training and Testing Datasets From Source Data

Create additional modified datasets, with their peptide and/or HLA sequences randomly sampled, to be tested with the [TransPHLA-AOMP](https://github.com/a96123155/TransPHLA-AOMP) model. Run the `scripts/8_transphla_randomize_input.py` script:

```bash
python scripts/8_transphla_randomize_input.py
```

##### 9. Create and Test Models

Open the newly-created `pHLAIformer_src.ipynb` Jupyter notebook and make the following changes to the code (cell and line numbering accounts for subsequently added items):

1. Add a new cell after cell no. 2 and paste the following code:

    ```python
    # Root TransPHLA-AOMP directory
    transphla_dir = "/path/to/TransPHLA-AOMP"
    # Root AImmuno analysis directory
    aimmuno_dir = "/path/to/aimmuno"
    output_data_dir = aimmuno_dir + "/data"
    models_dir = output_data_dir + "/models"
    results_root_dir = output_data_dir + "/results"
    model_type = "transphla"
    ```

2. Add a new cell after cell no. 3 and paste the following code:

    ```python
    suffix = "_src"
    # suffix = "_rand_hla_train_only"
    # suffix = "_rand_hla_test_only"
    # suffix = "_rand_hla_all"
    # suffix = "_rand_pep_all"
    # suffix = "_rand_pep_test_only"
    # suffix = "_rand_pep_train_only"
    ```

3. Change the contents of cell no. 5 to the following:

    ```python
    hla_sequence = pd.read_csv(transphla_dir + "/Dataset/common_hla_sequence.csv")
    ```

4. Change line no. 6 of cell no. 10 to the following:

    ```python
    vocab = np.load(transphla_dir + "/TransPHLA-AOMP/vocab_dict.npy", allow_pickle=True).item()
    ```

5. Change line no. 1 of cell no. 11 to the following:

    ```python
    def data_with_loader(type_="train", fold=None, batch_size=1024, f_suffix="_src"):
    ```

6. Add a new line after line no. 1 of cell no. 11 and paste the following:

    ```python
    if f_suffix == "_src":
        f_suffix = ""
    ```

7. Change line no. 5 of cell no. 11 to the following:

    ```python
    data = pd.read_csv(transphla_dir + f"/Dataset/{type_}_set{f_suffix}.csv", index_col=0)
    ```

8. Change line no. 7 of cell no. 11 to the following:

    ```python
    data = pd.read_csv(transphla_dir + f"/Dataset/train_data_fold{fold}{f_suffix}.csv", index_col=0)
    ```

9. Change line no. 9 of cell no. 11 to the following:

    ```python
    data = pd.read_csv(transphla_dir + f"/Dataset/val_data_fold{fold}{f_suffix}.csv", index_col=0)
    ```

10. Add a new cell after cell no. 11 and paste the following code:

    ```python
    def save_results():
        fold_n = 0
        results_dir = results_root_dir + f"/{model_type}{suffix}"
        roc_col_x = "x"
        roc_col_y = "y"
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)
        # Save metrics (accuracy, F1 score, etc.)
        train_metrics = performances_to_pd(train_fold_metrics_list)
        train_metrics.to_csv(results_dir + f"/train_metrics_{model_type}{suffix}_layer{n_layers}_multihead{n_heads}.csv")
        val_metrics = performances_to_pd(val_fold_metrics_list)
        val_metrics.to_csv(results_dir + f"/val_metrics_{model_type}{suffix}_layer{n_layers}_multihead{n_heads}.csv")
        indep_metrics = performances_to_pd(independent_fold_metrics_list)
        indep_metrics.to_csv(results_dir + f"/indep_metrics_{model_type}{suffix}_layer{n_layers}_multihead{n_heads}.csv")
        ext_metrics = performances_to_pd(external_fold_metrics_list)
        ext_metrics.to_csv(results_dir + f"/ext_metrics_{model_type}{suffix}_layer{n_layers}_multihead{n_heads}.csv")
        # Save ROC curve coordinates and AUC scores
        result_ys = {"train": ys_train_fold_dict[fold_n], "val": ys_val_fold_dict[fold_n], "indep": ys_independent_fold_dict[fold_n], "ext": ys_external_fold_dict[fold_n]}
        for r_type, r in result_ys.items():
            r_fpr, r_tpr, _ = metrics.roc_curve(r[0], r[2])
            r_auc = metrics.roc_auc_score(r[0], r[2])
            roc = pd.DataFrame({roc_col_x: r_fpr, roc_col_y: r_tpr})
            roc.to_csv(results_dir + f"/{r_type}_roc_{model_type}{suffix}_layer{n_layers}_multihead{n_heads}.csv", index=False)
            with open(results_dir + f"/{r_type}_auc_{model_type}{suffix}_layer{n_layers}_multihead{n_heads}.txt", "w") as handle:
                handle.write(str(r_auc))
    ```

11. Change line no. 1 of cell no. 13 to the following:

    ```python
    independent_data, independent_pep_inputs, independent_hla_inputs, independent_labels, independent_loader = data_with_loader(type_="independent", fold=None, batch_size=batch_size, f_suffix=suffix)
    ```

12. Change line no. 2 of cell no. 13 to the following:

    ```python
    external_data, external_pep_inputs, external_hla_inputs, external_labels, external_loader = data_with_loader(type_="external", fold=None, batch_size=batch_size, f_suffix=suffix)
    ```

13. Change line no. 12 of cell no. 14 to the following:

    ```python
    train_data, train_pep_inputs, train_hla_inputs, train_labels, train_loader = data_with_loader(type_="train", fold=fold, batch_size=batch_size, f_suffix=suffix)
    ```

14. Change line no. 13 of cell no. 14 to the following:

    ```python
    val_data, val_pep_inputs, val_hla_inputs, val_labels, val_loader = data_with_loader(type_="val", fold=fold, batch_size=batch_size, f_suffix=suffix)
    ```

15. Change line no. 22 of cell no. 14 to the following:

    ```python
    dir_saver = models_dir + f"/model_{model_type}{suffix}
    ```

16. Change line no. 23 of cell no. 14 to the following:

    ```python
    path_saver = dir_saver + f"/model_{model_type}{suffix}_layer{n_layers}_multihead{n_heads}_fold{fold}.pkl"
    ```

17. Change line no. 34 of cell no. 14 to the following:

    ```python
    metrics_ep_avg = np.nansum(metrics_val[:4]) / 4
    ```

18. After line no. 67 of cell no. 14 add the following line (should be inside the first `for` loop):

    ```python
    save_results()
    ```

Now, by changing which line is left uncommented in cell no. 4 (`suffix` definitions) we can train and test different models.

#### Synthetic Data for TransPHLA-AOMP

Now, we use the same synthetic data generated earlier and create new [TransPHLA-AOMP](https://github.com/a96123155/TransPHLA-AOMP) models with it.

##### 10. Create and Test Models

Open the newly-created `pHLAIformer_synthetic.ipynb` Jupyter notebook and make the following changes to the code (cell and line numbering accounts for subsequently added items):

1. Add a new cell after cell no. 2 and paste the following code:

    ```python
    # Root TransPHLA-AOMP directory
    transphla_dir = "/path/to/TransPHLA-AOMP"
    # Root AImmuno analysis directory
    aimmuno_dir = "/path/to/aimmuno"
    output_data_dir = aimmuno_dir + "/data"
    models_dir = output_data_dir + "/models"
    results_root_dir = output_data_dir + "/results"
    input_data_root_dir = aimmuno_dir + "/data/datasets"
    dataset_easy_true_train_file = input_data_root_dir + "/synthetic_easy/true_dataset_train.csv"
    dataset_easy_true_test_file = input_data_root_dir + "/synthetic_easy/true_dataset_test.csv"
    dataset_easy_decoys_train_file = input_data_root_dir + "/synthetic_easy/decoys_dataset_train.csv"
    dataset_easy_decoys_test_file = input_data_root_dir + "/synthetic_easy/decoys_dataset_test.csv"
    dataset_hard_true_train_file = input_data_root_dir + "/synthetic_hard/true_dataset_train.csv"
    dataset_hard_true_test_file = input_data_root_dir + "/synthetic_hard/true_dataset_test.csv"
    dataset_hard_decoys_train_file = input_data_root_dir + "/synthetic_hard/decoys_dataset_train.csv"
    dataset_hard_decoys_test_file = input_data_root_dir + "/synthetic_hard/decoys_dataset_test.csv"
    dataset_src_col_pep = "peptide"
    dataset_src_col_len = "length"
    dataset_src_col_label = "label"
    dataset_src_col_hla_seq = "HLA_sequence"
    dataset_true_col_pep = "peptide"
    dataset_true_col_hla_seq = "mhc"
    dataset_true_col_label = "label"
    dataset_decoys_col_pep = "peptide_rand"
    dataset_decoys_col_hla_seq = "mhc_rand"
    dataset_decoys_col_label = "label"
    model_type = "transphla"
    ```

2. Add a new cell after cell no. 3 and paste the following code:

    ```python
    dataset_type = "_synthetic_easy"
    # dataset_type = "_synthetic_hard"
    ```

3. Change the contents of cell no. 5 to the following:

    ```python
    hla_sequence = pd.read_csv(transphla_dir + "/Dataset/common_hla_sequence.csv")
    ```

4. Change line no. 1 of cell no. 10 to the following:

    ```python
    pep_max_len = 13
    ```

    Note: this variable, as well as the `hla_max_len` variable, are dependent on the actual maximum lengths of peptide and HLA sequences, respectively, within the datasets. This configuration is tuned to the data published in the [data repository](https://doi.org/10.5281/zenodo.14526105) and should be changed if the data change. `pep_max_len` should be set to maximum peptide length + 1, while `hla_max_len` to exactly the maximum HLA sequence length.

5. Change line no. 6 of cell no. 10 to the following:

    ```python
    vocab = np.load(transphla_dir + "/TransPHLA-AOMP/vocab_dict.npy", allow_pickle=True).item()
    ```

6. Add a new cell after cell no. 10 and paste the following code:

    ```python
    if dataset_type == "_synthetic_easy":
        dataset_true_train = pd.read_csv(dataset_easy_true_train_file, index_col=0)
        dataset_true_test = pd.read_csv(dataset_easy_true_test_file, index_col=0)
        dataset_decoys_train = pd.read_csv(dataset_easy_decoys_train_file, index_col=0)
        dataset_decoys_test = pd.read_csv(dataset_easy_decoys_test_file, index_col=0)
    elif dataset_type == "_synthetic_hard":
        dataset_true_train = pd.read_csv(dataset_hard_true_train_file, index_col=0)
        dataset_true_test = pd.read_csv(dataset_hard_true_test_file, index_col=0)
        dataset_decoys_train = pd.read_csv(dataset_hard_decoys_train_file, index_col=0)
        dataset_decoys_test = pd.read_csv(dataset_hard_decoys_test_file, index_col=0)
    dataset_true_train = dataset_true_train.rename(columns={dataset_true_col_hla_seq: dataset_src_col_hla_seq})
    dataset_true_test = dataset_true_test.rename(columns={dataset_true_col_hla_seq: dataset_src_col_hla_seq})
    dataset_decoys_train = dataset_decoys_train[[dataset_decoys_col_pep, dataset_decoys_col_label, dataset_decoys_col_hla_seq]].rename(columns={dataset_decoys_col_pep: dataset_src_col_pep, dataset_decoys_col_hla_seq: dataset_src_col_hla_seq})
    dataset_decoys_test = dataset_decoys_test[[dataset_decoys_col_pep, dataset_decoys_col_label, dataset_decoys_col_hla_seq]].rename(columns={dataset_decoys_col_pep: dataset_src_col_pep, dataset_decoys_col_hla_seq: dataset_src_col_hla_seq})
    dataset_train_src = pd.concat([dataset_true_train, dataset_decoys_train], ignore_index=True, axis=0)
    dataset_test = pd.concat([dataset_true_test, dataset_decoys_test], ignore_index=True, axis=0)
    del dataset_true_train, dataset_true_test, dataset_decoys_train, dataset_decoys_test
    dataset_train_src[dataset_src_col_label] = dataset_train_src[dataset_src_col_label].replace({True: 1, False: 0})
    dataset_test[dataset_src_col_label] = dataset_test[dataset_src_col_label].replace({True: 1, False: 0})
    dataset_train_src[dataset_src_col_len] = dataset_train_src[dataset_src_col_pep].str.len()
    dataset_test[dataset_src_col_len] = dataset_test[dataset_src_col_pep].str.len()
    dataset_train, dataset_validation = train_test_split(dataset_train_src, random_state=42, test_size=0.2, stratify=dataset_train_src[dataset_src_col_label])
    del dataset_train_src
    dataset_train.reset_index(drop=True, inplace=True)
    dataset_validation.reset_index(drop=True, inplace=True)
    dataset_test.reset_index(drop=True, inplace=True)
    ```

7. Change line no. 3 of cell no. 12 to the following:

    ```python
    data = dataset_test
    ```

8. Change line no. 5 of cell no. 12 to the following:

    ```python
    data = dataset_train
    ```

9. Change line no. 7 of cell no. 12 to the following:

    ```python
    data = dataset_validation
    ```

10. Add a new cell after cell no. 12 and paste the following code:

    ```python
    def save_results():
        fold_n = 0
        results_dir = results_root_dir + f"/{model_type}{dataset_type}"
        roc_col_x = "x"
        roc_col_y = "y"
        if not os.path.exists(results_dir):
            os.makedirs(results_dir)
        # Save metrics (accuracy, F1 score, etc.)
        train_metrics = performances_to_pd(train_fold_metrics_list)
        train_metrics.to_csv(results_dir + f"/train_metrics_{model_type}{dataset_type}_layer{n_layers}_multihead{n_heads}.csv")
        val_metrics = performances_to_pd(val_fold_metrics_list)
        val_metrics.to_csv(results_dir + f"/val_metrics_{model_type}{dataset_type}_layer{n_layers}_multihead{n_heads}.csv")
        ext_metrics = performances_to_pd(external_fold_metrics_list)
        ext_metrics.to_csv(results_dir + f"/test_metrics_{model_type}{dataset_type}_layer{n_layers}_multihead{n_heads}.csv")
        # Save ROC curve coordinates and AUC scores
        result_ys = {"train": ys_train_fold_dict[fold_n], "val": ys_val_fold_dict[fold_n], "test": ys_external_fold_dict[fold_n]}
        for r_type, r in result_ys.items():
            r_fpr, r_tpr, _ = metrics.roc_curve(r[0], r[2])
            r_auc = metrics.roc_auc_score(r[0], r[2])
            roc = pd.DataFrame({roc_col_x: r_fpr, roc_col_y: r_tpr})
            roc.to_csv(results_dir + f"/{r_type}_roc_{model_type}{dataset_type}_layer{n_layers}_multihead{n_heads}.csv", index=False)
            with open(results_dir + f"/{r_type}_auc_{model_type}{dataset_type}_layer{n_layers}_multihead{n_heads}.txt", "w") as handle:
                handle.write(str(r_auc))
    ```

11. Comment out line no. 1 of cell no. 14.

12. Change line no. 22 of cell no. 15 to the following:

    ```python
    dir_saver = models_dir + f"/model_{model_type}{dataset_type}
    ```

13. Change line no. 23 of cell no. 15 to the following:

    ```python
    path_saver = dir_saver + f"/model_{model_type}{dataset_type}_layer{n_layers}_multihead{n_heads}_fold{fold}.pkl"
    ```

14. Change line no. 34 of cell no. 15 to the following:

    ```python
    metrics_ep_avg = np.nansum(metrics_val[:4]) / 4
    ```

15. Comment out line no. 54 of cell no. 15.

16. Comment out line no. 59 of cell no. 15.

17. Change line no. 62 of cell no. 15 to the following:

    ```python
    ys_train_fold_dict[fold], ys_val_fold_dict[fold], ys_external_fold_dict[fold] = ys_res_train, ys_res_val, ys_res_external
    ```

18. Change line no. 64 of cell no. 15 to the following:

    ```python
    loss_train_fold_dict[fold], loss_val_fold_dict[fold], loss_external_fold_dict[fold] = loss_res_train_list, loss_res_val_list, loss_res_external_list
    ```

19. After line no. 67 of cell no. 15 add the following line (should be inside the first `for` loop):

    ```python
    save_results()
    ```

20. Comment out lines no. 1 and 2 of cell no. 16.

Now, by changing which line is left uncommented in cell no. 4 (`dataset_type` definitions) we can train and test different models.

## Citing This Work

Please cite the described work and relevant paper as:

```bibtex
@article{aimmuno-paper,
  author  = {Surname, Name},
  title   = {{Artificial Intelligence Perceives Marginal Gains from MHC Haplotype Data in Antigen Presentation Predictions}},
  journal = {X},
  year    = {X},
  volume  = {X},
  number  = {X},
  pages   = {X--X},
  doi     = {X},
  url     = {}
}
```

Also, please add citations for other items associated with the described pipeline as:

```bibtex
@dataset{aimmuno-analysis-files,
  author    = {Mastromattei, Michele and Palkowski, Aleksander and Nourbakhsh, Aria and Alfaro, Javier Antonio and Zanzotto, Fabio Massimo},
  title     = {{Artificial intelligence perceives marginal gains from MHC haplotype data: analysis and datasets that challenge models}},
  month     = {dec},
  year      = {2024},
  publisher = {Zenodo},
  version   = {1.0.0},
  doi       = {10.5281/zenodo.14526105},
  url       = {https://doi.org/10.5281/zenodo.14526105}
}
@software{aimmuno-analysis-code,
  author    = {Mastromattei, Michele and Palkowski, Aleksander and Nourbakhsh, Aria and Alfaro, Javier Antonio and Zanzotto, Fabio Massimo},
  title     = {{AImmuno: A pipeline for testing various MHC class I antigen binding prediction models}},
  month     = {dec},
  year      = {2024},
  publisher = {Zenodo},
  version   = {v1.0.1},
  doi       = {10.5281/zenodo.14536592},
  url       = {https://doi.org/10.5281/zenodo.14536592}
}
```

## Acknowledgements

The described pipeline uses and/or references the following external libraries, packages, and other software:

- [Jupyter](https://jupyter.org)
- [Keras](https://keras.io)
- [Matplotlib](https://matplotlib.org)
- [NumPy](https://numpy.org)
- [pandas](https://pandas.pydata.org)
- [PyTorch](https://pytorch.org)
- [scikit-learn](https://scikit-learn.org)
- [SciPy](https://scipy.org)
- [TAPE](https://github.com/songlab-cal/tape)
- [tqdm](https://github.com/tqdm/tqdm)
- [TransPHLA-AOMP](https://github.com/a96123155/TransPHLA-AOMP)

We wish to thank all their contributors and maintainers!

This project has received funding from the European Union's Horizon 2020/H2020-SCI-FA-DTS-2020-1 research and innovation programme under the Knowledge At the Tip of Your fingers: Clinical Knowledge for Humanity project (no. 101017453).

## License and Disclaimer

Copyright (c) 2024 immuno-informatics

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

### Third-Party Software

The software, libraries, code, or data from third parties mentioned in the [Acknowledgements](#acknowledgements) section above may come with their own terms and conditions or licensing requirements. When using this third-party software, libraries, code, or data it's essential to adhere to these terms. Ensure you understand and can follow any relevant restrictions or terms and conditions prior to using them.
