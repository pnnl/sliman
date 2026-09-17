"""
sliman/annotation/calibration/_plots.py

Dylan Ross (dylan.ross@pnnl.gov)

    Internal module with plotting functions related to reference-free CCS stuff
"""


import numpy as np
from matplotlib import pyplot as plt, rcParams


rcParams["font.size"] = 8
rcParams["font.family"] = "Roboto Condensed"


Z2C = {1: "#256EFF", 2: "#3DDC97", 3: "#46237A"}


def plot_candidate_mz_at_ccs(mzs, dts, ccss):
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    for z in [1, 2, 3]:
        ax.scatter(mzs[z], dts[z], marker='.', s=16, c=Z2C[z], linewidth=0, label=f"z={z}")
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.legend(frameon=False)
    ax.set_xlabel('m/z')
    ax.set_ylabel('arrival time (s)')
    plt.savefig("_figures/mz_at_calibration_candidates.png", dpi=400, bbox_inches="tight")
    plt.close()
    # plot the potential calibrant points (dt, ccs)
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    for z in [1, 2, 3]:
        ax.scatter(dts[z], ccss[z], marker='.', s=16, c=Z2C[z], linewidth=0, label=f"z={z}")
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.legend(frameon=False)
    ax.set_ylabel(r'DL-pred. CCS ($\AA^2$)')
    ax.set_xlabel('arrival time (s)')
    plt.savefig("_figures/at_ccs_calibration_candidates.png", dpi=400, bbox_inches="tight")
    plt.close()


def plot_candidate_at_hists(dts):
    fig, ax = plt.subplots(figsize=(1.5, 1.5))
    bins = np.arange(1., 2., 0.05)
    for z in [1, 2, 3]:
        ax.hist(dts[z], bins=bins, edgecolor=Z2C[z], label=f"z={z}", histtype="step", linewidth=1.0),
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.legend(frameon=False)
    ax.set_ylabel('#')
    ax.set_xlabel('arrival time (s)')
    plt.savefig("_figures/at_hist_calibration_candidates.png", dpi=400, bbox_inches="tight")
    plt.close()


def plot_calibrant_mz_at_ccs(mzs_all, dts_all, ccss_all, mzs, dts, ccss, sampled):
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    mzs_smp = {z: [a[_] for _ in sampled[z]] for z, a in mzs.items()}
    dts_smp = {z: [a[_] for _ in sampled[z]] for z, a in dts.items()}
    ccss_smp = {z: [a[_] for _ in sampled[z]] for z, a in ccss.items()}
    z2m = {1: "o", 2: "^", 3: "s"}
    z2sz = {1: 6, 2: 10, 3: 6}
    for z in [1, 2, 3]:
        ax.scatter(dts_all[z], mzs_all[z], marker='.', s=4, c=Z2C[z], linewidth=0, zorder=-1)
        #ax.scatter(dts_smp[z], mzs_smp[z], marker='x', s=8, c="k", linewidth=0.75, label=f"z={z}")
        ax.scatter(dts_smp[z], mzs_smp[z], marker=z2m[z], s=z2sz[z], edgecolor=Z2C[z], linewidth=0.75, 
                   color="w", label=f"z={z}")
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.legend(frameon=False)
    ax.set_ylabel('m/z')
    ax.set_xlabel('arrival time (s)')
    plt.savefig("_figures/mz_at_calibrants.png", dpi=400, bbox_inches="tight")
    plt.close()
    # plot the potential calibrant points (dt, ccs)
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    for z in [1, 2, 3]:
        ax.scatter(dts_all[z], ccss_all[z], marker='.', s=4, c=Z2C[z], linewidth=0, zorder=-1, alpha=0.3)
        #ax.scatter(dts_smp[z], ccss_smp[z], marker='x', s=8, c="k", linewidth=0.75, label=f"z={z}")
        ax.scatter(dts_smp[z], ccss_smp[z], marker=z2m[z], s=z2sz[z], edgecolor=Z2C[z], linewidth=0.75, 
                   color="w", label=f"z={z}")
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.legend(frameon=False)
    ax.set_ylabel(r'DL-pred. CCS ($\AA^2$)')
    ax.set_xlabel('arrival time (s)')
    plt.savefig("_figures/at_ccs_calibrants.png", dpi=400, bbox_inches="tight")
    plt.close()


def plot_calibrant_at_ccs(mzs, dts, ccss, cals):
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    ccss_cal = {z: cals[z].calibrated_ccs(mzs[z], dts[z]) for z in [1, 2, 3]}
    # plot the potential calibrant points (dt, ccs)
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    z2m = {1: "o", 2: "^", 3: "s"}
    z2sz = {1: 6, 2: 10, 3: 6}
    for z in [1, 2, 3]:
        dt = dts[z]
        idx = np.argsort(dt)
        ax.scatter(dt[idx], ccss[z][idx], marker=z2m[z], s=z2sz[z], edgecolor=Z2C[z], 
                   linewidth=0.75, color="#00000000")
        ax.plot(dt[idx], ccss_cal[z][idx], ls="--", lw=1, c=Z2C[z], label=f"z={z}")
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.legend(frameon=False)
    ax.set_ylabel(r'DL-pred. CCS ($\AA^2$)')
    ax.set_xlabel('arrival time (s)')
    plt.savefig("_figures/at_ccs_calibrated.png", dpi=400, bbox_inches="tight")
    plt.close()


def plot_filtered_at_ccs(dts_all, ccss_all, dts_filtered, ccss_filtered, figname):
    # plot the potential calibrant points (dt, ccs)
    fig, ax = plt.subplots(figsize=(3, 3))
    z2m = {1: "o", 2: "^", 3: "s"}
    z2sz = {1: 6, 2: 10, 3: 6}
    for z in [1, 2, 3]:
        ax.scatter(dts_all[z], ccss_all[z], marker='.', s=4, c="#969696", 
                   linewidth=0, zorder=-1, alpha=1)
        # ax.scatter(dts_smp[z], ccss_smp[z], marker='x', s=8, c="k", linewidth=0.75, label=f"z={z}")
        # ax.scatter(dts_filtered[z], ccss_filtered[z], marker=z2m[z], s=z2sz[z], edgecolor=Z2C[z], 
        #            linewidth=0.75, color="#00000000", label=f"z={z}")
        ax.scatter(dts_filtered[z], ccss_filtered[z], marker=".", s=4, c=Z2C[z], 
                   linewidth=0., label=f"z={z}")
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.legend(frameon=False)
    ax.set_ylabel(r'DL-pred. CCS ($\AA^2$)')
    ax.set_xlabel('arrival time (s)')
    plt.savefig(figname, dpi=400, bbox_inches="tight")
    plt.show()
    plt.close()


