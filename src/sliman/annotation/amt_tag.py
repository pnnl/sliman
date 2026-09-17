"""
sliman/annotation/amt_tag.py

Dylan Ross (dylan.ross@pnnl.gov)

    Utilities for constructing/interacting with AMT tag database
"""


import os
from typing import List, Dict, Tuple, Set
import csv
import sqlite3

import numpy as np

from sliman._util import tol_from_ppm


# AMT tag database schema
_AMT_TAG_DB_SCHEMA = """--sqlite3
CREATE TABLE Datasets (
    dataset_id INTEGER PRIMARY KEY,
    dataset_label TEXT NOT NULL
);--sqlite3

CREATE TABLE Peptides (
    peptide_id INTEGER PRIMARY KEY,
    sequence TEXT NOT NULL,
    z INT NOT NULL,
    mz REAL NOT NULL
);--sqlite3

CREATE TABLE Tags (
    dataset_id INT NOT NULL,
    peptide_id INT NOT NULL
);--sqlite3

CREATE TABLE Targets (
    target_id INTEGER PRIMARY KEY,
    target_labels TEXT NOT NULL,
    target_mz REAL NOT NULL,
    peptide_ids TEXT NOT NULL,
    dataset_ids TEXT NOT NULL
)
;"""


def init_amt_tag_db(amt_tag_db: str, overwrite: bool = False) -> None :
    """
    Initialize an AMT tag database
    
    Parameters
    ----------
    amt_tag_db : ``str``
        AMT tag database file
    overwrite : ``bool``, default=False
        if True, overwrite database file if it exists otherwise raise a RuntimeError
    """
    # check for already existing database file
    if os.path.isfile(amt_tag_db):
        if overwrite:
            os.remove(amt_tag_db)
        else:
            msg = f"AMT tag database file {amt_tag_db} exists but overwrite is False"
            raise RuntimeError(msg)
    # initialize the database
    con = sqlite3.connect(amt_tag_db)
    con.executescript(_AMT_TAG_DB_SCHEMA)
    con.commit()
    con.close()


def _fetch_existing_values(cur: sqlite3.Cursor
                           ) -> Tuple[Dict[str, int],
                                      Dict[Tuple[str, int], int],
                                      Set[Tuple[int, int]]] :
    """
    Fetch all of the existing Datasets and Peptides values from the AMT tag database
    """
    datasets = {
        dataset_label: dataset_id 
        for dataset_label, dataset_id
        in cur.execute("SELECT * FROM Datasets").fetchall()
    }
    peptides = {
        (sequence, z): peptide_id
        for peptide_id, sequence, z, _
        in cur.execute("SELECT * FROM Peptides").fetchall()
    }
    tags = set(cur.execute("SELECT * FROM Tags").fetchall())
    return datasets, peptides, tags


def _compute_mz(mz1: float, z: int) -> float :
    """ compute the m/z for the desired charge state (z) from the +1 m/z (mz1) """
    # only positive charge states are currently implemented
    assert z > 0, "only positive charge states currently implemented"
    # if +1 just give back the +1 m/z
    if z == 1:
        return mz1
    # compute other charge states from +1 m/z
    h = 1.00783 * (z - 1)
    return (mz1 + h) / z


def add_spec_lib_to_amt_tag_db(amt_tag_db: str, 
                               spec_lib_file: str, 
                               dataset_labels: List[str], 
                               charge_states: List[int]
                               ) -> None :
    """
    Add a spectral library (output from: ???) to the AMT tag database

    Takes all unique peptide sequences from the selected Dataset labels and recomputes
    m/zs for all of the specified charge states. Resulting targets are stored in the 
    AMT tag database file.

    Parameters
    ----------
    amt_tag_db : ``str``
        AMT tag database file
    spec_lib_file : ``str``
        Input spectral library file (output from: ???)
    dataset_labels : ``list(str)``
        list of dataset labels to include from the input spectral library
    charge_states : ``list(int)``
        list of charge states to include, these are all recomputed from the +1 m/z
        that is included in the spectral library
    """
    # TODO: convert these asserts to raised exceptions
    # only positive charge states are currently implemented
    for z in charge_states:
        assert z > 0, "only positive charge states currently implemented"
    # AMT tag database must exist
    assert os.path.isfile(amt_tag_db), f"AMT tag database file: {amt_tag_db} not found"
    # spec lib file must exist
    assert os.path.isfile(spec_lib_file), f"spec lib file: {spec_lib_file} not found"
    # connect to AMT tag DB
    con = sqlite3.connect(amt_tag_db)
    cur = con.cursor()
    # fetch existing values from the DB
    cached_datasets, cached_peptides, cached_tags = _fetch_existing_values(cur)
    # iterate through spectral library
    # grab the columns: "Dataset", "MH", "Peptide"
    # TODO: Should we store any additional metadata?
    with open(spec_lib_file, "r") as slf:
        _ = next(slf)
        rdr = csv.reader(slf, delimiter="\t")
        for dataset_label, _, _, _, _, _, _, _, mz1, sequence, *_ in rdr:
            if dataset_label in dataset_labels:
                # get a dataset identifier
                dataset_id = cached_datasets.get(dataset_label)
                if dataset_id is None:
                    cur.execute("INSERT INTO Datasets VALUES (?,?)", (None, dataset_label))
                    dataset_id = cur.lastrowid
                    # to satisfy type checking:
                    #assert dataset_id is not None
                    cached_datasets[dataset_label] = dataset_id
                # compute desired m/zs
                for z in charge_states:
                    mz = _compute_mz(float(mz1), z)
                    # get a peptide identifier
                    peptide_id = cached_peptides.get((sequence, z))
                    if peptide_id is None:
                        cur.execute("INSERT INTO Peptides VALUES (?,?,?,?)", (None, sequence, z, mz))
                        peptide_id = cur.lastrowid
                        # to satisfy type checking:
                        #assert peptide_id is not None
                        cached_peptides[(sequence, z)] = peptide_id
                    # map the dataset identifier to peptide identifier
                    if (dataset_id, peptide_id) not in cached_tags:
                        cur.execute("INSERT INTO Tags VALUES (?,?)", (dataset_id, peptide_id))
                        cached_tags.add((dataset_id, peptide_id))
    # clean up
    con.commit()
    con.close()


def _consolidate_group_ds_ids(group_ds_ids: List[str]) -> str :
    """ consolidate dataset IDs from the group to only include unique IDs """
    dataset_ids = []
    for group_ds_id in group_ds_ids:
        dataset_ids += group_ds_id.split(",")
    return ",".join(sorted(list(set(dataset_ids))))


def create_targets(amt_tag_db: str,
                   ppm: float,
                   sel_dsets: List[str],
                   clear_targets: bool = False,
                   ) -> None :
    """
    Create targets used for data analysis by combining peptides having m/z within 
    a specified PPM. Each target has a target_id, a set of (sequence, z) for 
    corresponding peptides as labels, mz, and all associated peptide IDs. Targets are 
    created from peptides from selected datasets

    Parameters
    ----------
    amt_tag_db : ``str``
        AMT tag database file
    ppm : ``float``
        ppm threshold for combining peptides by m/z
    sel_dsets : ``list(str)``
        list of selected dataset labels to create targets from
    clear_targets : ``bool``, default=False
        if True, then clear out the Targets table in the database before starting
    """
    # AMT tag database must exist
    assert os.path.isfile(amt_tag_db), f"AMT tag database file: {amt_tag_db} not found" 
    # connect to AMT tag DB
    con = sqlite3.connect(amt_tag_db)
    con.row_factory = sqlite3.Row  # enable accessing query result row elements by key
    cur = con.cursor()
    # optionally clear out the Targets table first
    if clear_targets:
        cur.execute("DELETE FROM Targets;")
    # 
    # create a map between dataset identifiers and dataset labels
    dslbl_to_dsid: Dict[str, str] = {
        row["dataset_label"]: str(row["dataset_id"])
        for row in cur.execute("SELECT dataset_id, dataset_label FROM Datasets").fetchall()
    }
    # get the peptides list
    qry_sel = """--sqlite3
    SELECT 
        peptide_id, 
        sequence, 
        z, 
        mz, 
        GROUP_CONCAT(dataset_id) AS dataset_ids 
    FROM 
        Peptides 
        JOIN Tags USING(peptide_id) 
        JOIN Datasets USING(dataset_id) 
    WHERE
        dataset_id IN ({})
    GROUP BY 
        peptide_id 
    ORDER BY 
        mz
    ;""".format(",".join([dslbl_to_dsid[lbl] for lbl in sel_dsets]))
    peps: List[sqlite3.Row] = cur.execute(qry_sel).fetchall()
    # peptide_id, sequence, z, mz
    # consolidate the peptides based on m/z
    i: int = 1
    group_mz: float = peps[0]["mz"] + tol_from_ppm(peps[0]["mz"], ppm)
    group_labels: List[str] = [f"{peps[0]["sequence"]}_{peps[0]["z"]}"] 
    group_mzs: List[float] = [peps[0]["mz"]]
    group_pep_ids: List[int] = [peps[0]["peptide_id"]]
    group_ds_ids: List[str] = [peps[0]["dataset_ids"]]
    while i < len(peps):
        if peps[i]["mz"] <= group_mz + tol_from_ppm(group_mz, ppm):
            group_labels.append(f"{peps[i]["sequence"]}_{peps[i]["z"]}")
            group_mzs.append(peps[i]["mz"])
            group_pep_ids.append(peps[i]["peptide_id"])
            group_ds_ids.append(peps[i]["dataset_ids"])
        else:
            # store the current group 
            cur.execute("INSERT INTO Targets VALUES (?,?,?,?,?)",
                        (None,                                       # target_id
                         "|".join(group_labels),                     # target_labels
                         np.mean(group_mzs),                         # target_mz
                         ",".join([str(_) for _ in group_pep_ids]),  # peptide_ids
                         _consolidate_group_ds_ids(group_ds_ids)))   # dataset_ids
            # reset the group_mz to the new peptide at i
            group_mz = peps[i]["mz"] + tol_from_ppm(peps[i]["mz"], ppm)
            group_labels = [f"{peps[i]["sequence"]}_{peps[i]["z"]}"]
            group_mzs = [peps[i]["mz"]]
            group_pep_ids = [peps[i]["peptide_id"]]
            group_ds_ids = [peps[i]["dataset_ids"]]
        # increment index
        i += 1
    # store the final group
    group_ds_ids = list(set(group_ds_ids))
    cur.execute("INSERT INTO Targets VALUES (?,?,?,?,?)",
                (None,                                       # target_id
                 "|".join(group_labels),                     # target_labels
                 np.mean(group_mzs),                         # target_mz
                 ",".join([str(_) for _ in group_pep_ids]),  # peptide_ids
                 _consolidate_group_ds_ids(group_ds_ids)))   # dataset_ids
    # clean up
    con.commit()
    con.close()
