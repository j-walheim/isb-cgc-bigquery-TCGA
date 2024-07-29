# %%
import os
import pandas as pd
from google.cloud import bigquery
from google.oauth2 import service_account
from dotenv import load_dotenv
import numpy as np
import matplotlib.pyplot as plt
from sksurv.nonparametric import kaplan_meier_estimator
from sksurv.linear_model import CoxPHSurvivalAnalysis
from sksurv.compare import compare_survival
from sksurv.util import Surv

# Load the environment variables from the .env file
load_dotenv()

gene = "PIK3CA"
disease = 'TCGA-LUAD'

def get_rnaseq_data(client, gene_name, case_barcode_list):
    query = f"""
    SELECT *
    FROM `isb-cgc-bq.TCGA.RNAseq_hg38_gdc_current`
    WHERE gene_name = '{gene_name}' AND case_barcode IN UNNEST({case_barcode_list})
    """
    return client.query(query).to_dataframe()


# Function to get clinical data (modified to include project_id)
def get_clinical_data(client, project_id=None):

    if project_id:
        query = f"""
        SELECT
            submitter_id,
            demo__vital_status as vital_status,
            demo__days_to_death as days_to_death,
            diag__days_to_last_follow_up as days_to_last_follow_up,
            proj__project_id as project_id
        FROM `isb-cgc-bq.TCGA_versioned.clinical_gdc_r37`
        WHERE proj__project_id = '{project_id}'
        """
    else:
        query = """
        SELECT
            submitter_id,
            demo__vital_status as vital_status,
            demo__days_to_death as days_to_death,
            diag__days_to_last_follow_up as days_to_last_follow_up,
            proj__project_id as project_id
        FROM `isb-cgc-bq.TCGA_versioned.clinical_gdc_r37`
        """
    return client.query(query).to_dataframe()


def get_somatic_mutation_data(client, gene_name, case_barcodes):
    query = f"""
    SELECT *
    FROM `isb-cgc-bq.TCGA.masked_somatic_mutation_hg38_gdc_current`
    WHERE Hugo_Symbol = '{gene_name}' AND case_barcode IN UNNEST({case_barcodes})
    """
    return client.query(query).to_dataframe()

credentials = service_account.Credentials.from_service_account_file(
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'],
    scopes=["https://www.googleapis.com/auth/cloud-platform"],
)
client = bigquery.Client(credentials=credentials, project=credentials.project_id)

# %%

df_clinical = get_clinical_data(client, project_id=disease)


# %%
case_barcodes = df_clinical['submitter_id'].tolist()
df_mutation = get_somatic_mutation_data(client, gene, case_barcodes)





# %%
# Prepare the data for survival analysis
df_clinical['time'] = df_clinical['days_to_death'].fillna(df_clinical['days_to_last_follow_up'])
df_clinical['event'] = (df_clinical['vital_status'] == 'Dead').astype(bool)

# Merge clinical data with mutation data
df_merged = df_clinical.merge(df_mutation[['case_barcode']], left_on='submitter_id', right_on='case_barcode', how='left', indicator=True)
df_merged['mutated'] = (df_merged['_merge'] == 'both').astype(bool)

# remove rows with NA in 'time' - i.e., there is no death and no follow-up time
df_merged = df_merged.dropna(subset=['time'])

# Ensure 'time' is numeric and positive
df_merged['time'] = pd.to_numeric(df_merged['time'], errors='coerce')
assert (df_merged['time'] > 0).all(), "Found rows with negative or zero 'time' values"

# Create structured array for scikit-survival
y = Surv.from_arrays(event=df_merged['event'], time=df_merged['time'])
X = df_merged[['mutated']].astype(int)  # Convert boolean to int for Cox model

# Compute Kaplan-Meier estimates with confidence intervals
def compute_km_with_ci(event, time):
    time_points, survival_prob, conf_int = kaplan_meier_estimator(event, time, conf_type="log-log")
    return time_points, survival_prob, conf_int

time_mutated, survival_prob_mutated, conf_int_mutated = compute_km_with_ci(y['event'][X['mutated'] == 1], y['time'][X['mutated'] == 1])
time_wildtype, survival_prob_wildtype, conf_int_wildtype = compute_km_with_ci(y['event'][X['mutated'] == 0], y['time'][X['mutated'] == 0])

# Plot
plt.figure(figsize=(10, 6))
plt.step(time_mutated, survival_prob_mutated, where="post", label=f'{gene} Mutated')
plt.fill_between(time_mutated, conf_int_mutated[0], conf_int_mutated[1], alpha=0.3)
plt.step(time_wildtype, survival_prob_wildtype, where="post", label=f'{gene} Wild-type')
plt.fill_between(time_wildtype, conf_int_wildtype[0], conf_int_wildtype[1], alpha=0.3)

plt.ylabel("Survival probability")
plt.xlabel("Time (days)")
plt.legend()
plt.title(f'Survival plot for {gene} mutations in {disease}')

# %%
