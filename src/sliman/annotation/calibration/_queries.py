"""
sliman/annotation/calibration/_queries.py

Dylan Ross (dylan.ross@pnnl.gov)

    define a bunch of SQL queries for the calibration part of the
    analysis because it uses a database to coordinate intermediates
    and results from all of the steps
"""


PEPTIDE_DB_SCHEMA = """--sqlite3
-- database schema for peptide database
-- store the feature info 
CREATE TABLE Features (
    feat_id INTEGER PRIMARY KEY,
    smpl TEXT NOT NULL,
    mz REAL NOT NULL,
    dt REAL NOT NULL,
    dt_fwhm REAL NOT NULL,
    pk_ht REAL NOT NULL,
    pk_snr REAL NOT NULL
)
;--sqlite3
-- store calibrated CCS values 
-- separately (map to feat_id)
CREATE TABLE CalibratedCCS (
    feat_id INT NOT NULL,
    cal_ccs REAL NOT NULL
)
;--sqlite3
-- store individual peptide info
CREATE TABLE Peptides (
    pep_id INTEGER PRIMARY KEY,
    seq TEXT NOT NULL,
    z INT NOT NULL
)
;--sqlite3
-- store DL-predicted CCS values
-- separately (map to pep_id)
CREATE TABLE PredictedCCS (
    pep_id INT NOT NULL,
    pred_ccs REAL NOT NULL
)
;--sqlite3
-- allow for many-to-many mapping between 
-- peptides and features
CREATE TABLE FeatsPepsMapped (
    feat_id INT NOT NULL,
    pep_id INT NOT NULL
)
;--sqlite3
-- when calibrant selection is performed, store
-- the selected calibrants in this table
CREATE TABLE Calibrants (
    z INT NOT NULL,
    seq TEXT NOT NULL,
    mz REAL NOT NULL,
    dt REAL NOT NULL,
    ref_ccs REAL NOT NULL
)
;--sqlite3
-- store the pickle files for fitted calibrations 
-- directly within the database
CREATE TABLE Calibrations (
    z INT NOT NULL,
    pickle BLOB NOT NULL
)
;--sqlite3
-- for debugging purposes keep a list of what
-- steps have been performed
CREATE TABLE Log (
    entry TEXT NOT NULL
)
;"""

QRY_INS_LOG_ENTRY = """--sqlite3
INSERT INTO
    Log
VALUES
    (?)
;"""

QRY_SEL_PEPS_1 = """--sqlite3
SELECT
    pep_id,
    seq,
    z
FROM 
    Peptides
;"""

QRY_INSERT_FEAT = """--sqlite3
INSERT INTO
    Features
VALUES
    (?,?,?,?,?,?,?)
;"""

QRY_INSERT_PEP = """--sqlite3
INSERT INTO
    Peptides
VALUES
    (?,?,?)
;"""

QRY_INSERT_FEATPEPMAP = """--sqlite3
INSERT INTO
    FeatsPepsMapped
VALUES
    (?,?)
;"""

QRY_INS_PEP_DL_CCS = """--sqlite3
INSERT INTO 
    PredictedCCS
VALUES
    (?,?)
;"""

QRY_SEL_CANDIDATES = """--sqlite3
SELECT 
    seq, 
    mz,
    dt, 
    pred_ccs
FROM 
    Peptides
    LEFT JOIN 
        FeatsPepsMapped
        USING(pep_id)
    LEFT JOIN
        Features
        USING(feat_id)
    LEFT JOIN
        PredictedCCS
        USING(pep_id)
GROUP BY
    seq,
    z
HAVING
    pep_id IN
        (
            SELECT 
                pep_id
            FROM 
                FeatsPepsMapped
                LEFT JOIN 
                    Features
                    USING(feat_id)
            GROUP BY
                feat_id
            HAVING
                COUNT(*)=1
        )
    AND z=?
    AND pred_ccs IS NOT NULL
ORDER BY
    dt
;"""

QRY_SEL_ALL_FEATS = """--sqlite3
SELECT 
    seq,
    mz,
    dt, 
    pred_ccs
FROM 
    Peptides
    LEFT JOIN 
        FeatsPepsMapped
        USING(pep_id)
    LEFT JOIN
        Features
        USING(feat_id)
    LEFT JOIN
        PredictedCCS
        USING(pep_id)
GROUP BY
    seq,
    z
HAVING
    z=?
ORDER BY
    dt
;"""

QRY_SEL_FEATS = """--sqlite3
SELECT 
    feat_id,
    mz,
    dt,
    z
FROM 
    Features 
    LEFT JOIN 
        FeatsPepsMapped 
        USING(feat_id) 
    JOIN 
        Peptides 
        USING(pep_id)
GROUP BY
    feat_id
;"""

QRY_INS_CALIBRANT = """--sqlite3
INSERT INTO 
    Calibrants
VALUES 
    (?,?,?,?,?)
;"""

QRY_INS_CCS = """--sqlite3
INSERT INTO 
    CalibratedCCS 
VALUES 
    (?,?)
;"""

QRY_SEL_FEATS_2 = """--sqlite3
SELECT
    feat_id,
    smpl,
    mz,
    dt,
    pk_ht,
    cal_ccs
FROM
    Features
    JOIN
        CalibratedCCS
        USING(feat_id)
;"""

QRY_SEL_PEPS_2 = """--sqlite3
SELECT
    seq,
    z,
    pred_ccs
FROM 
    Peptides
    JOIN
        PredictedCCS
        USING(pep_id)
    JOIN
        FeatsPepsMapped
        USING(pep_id)
WHERE 
    feat_id=?
GROUP BY
    seq,
    z
;"""

QRY_SEL_CALIBRANTS = """--sqlite3
SELECT
    mz, 
    dt,
    ref_ccs
FROM
    Calibrants
WHERE
    z=?
;"""

QRY_INS_CALIBRATION = """--sqlite3
INSERT INTO
    Calibrations
VALUES
    (?,?)
;"""

QRY_SEL_CALIBRATIONS = """--sqlite3
SELECT
    z, 
    pickle
FROM 
    Calibrations
;"""

QRY_RESET_CALIBRANTS = """--sqlite3
DELETE FROM 
    Calibrants
;"""

QRY_RESET_CALIBRATIONS = """--sqlite3
DELETE FROM
    Calibrations
;"""

QRY_RESET_CALIBRATED_CCS = """--sqlite3
DELETE FROM 
    CalibratedCCS
;"""
