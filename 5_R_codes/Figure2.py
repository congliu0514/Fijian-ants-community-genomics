import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

# Load the CSV file
file_path = "./Data/fsc_species_level_summary.csv"  
df = pd.read_csv(file_path)

# "T_m1_honeybee_t", "T_m2_heliconius_t", "T_m3_aphids_t", "T_m4_silkworm_t
selected_method = "T_m2_heliconius_t"
df_long = df[["Endemic_class", "NCUR_m1_honeybee", "NANC_m1_honeybee", selected_method]].copy()
df_long = df_long.rename(columns={selected_method: "T"})
df_long["T_method"] = selected_method  

df_long["Change"] = np.where(df_long["NCUR_m1_honeybee"] > df_long["NANC_m1_honeybee"], "Expansion", "Decline")
bins = np.arange(0, 3300, 150)
df_long["T_bin"] = pd.cut(df_long["T"], bins=bins, right=False, labels=bins[:-1])
grouped_summary = df_long.groupby(["Endemic_class", "T_bin", "Change"]).size().unstack(fill_value=0).reset_index()

species_groups = grouped_summary["Endemic_class"].unique()

fig, axes = plt.subplots(len(species_groups), 1, figsize=(10, len(species_groups) * 4), sharex=True)

if len(species_groups) == 1:
    axes = [axes]

for i, group in enumerate(species_groups):
    ax = axes[i]
    group_data = grouped_summary[grouped_summary["Endemic_class"] == group]

    expansion_counts = group_data["Expansion"] if "Expansion" in group_data else [0] * len(group_data)
    decline_counts = group_data["Decline"] if "Decline" in group_data else [0] * len(group_data)

    ax.bar(group_data["T_bin"], expansion_counts, width=150, label="Expansion", color="#2ca02c", edgecolor="black")
    if decline_counts.sum() > 0:  # Only plot the negative panel if there's a decline
        ax.bar(group_data["T_bin"], -decline_counts, width=150, label="Decline", color="#d62728", edgecolor="black")

    # Highlight the special time points (400, 1800, and 3000 years)
    ax.axvline(x=400, color="blue", linestyle="--", linewidth=1, label="400 years")
    ax.axvline(x=1800, color="orange", linestyle="--", linewidth=1, label="1800 years")
    ax.axvline(x=3000, color="purple", linestyle="--", linewidth=1, label="3000 years")

    ax.set_ylabel(f"{group} (Species Count)", fontsize=12)
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.grid(axis="y", linestyle="--", alpha=0.7)

    ax.set_ylim(-12, 15)
    
    ax.set_xlim(3300, -300)

fig.text(0.5, 0.04, "Time T (years)", ha="center", fontsize=14)
fig.text(0.04, 0.5, "Species Count", va="center", rotation="vertical", fontsize=14)

plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
output_pdf = "species_demographic_changes_T_m1_honeybee.pdf"
with PdfPages(output_pdf) as pdf:
    pdf.savefig(fig)

print(f"Plots saved to {output_pdf}")