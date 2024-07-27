# %%
import os
import pandas as pd
from google.cloud import bigquery
from google.oauth2 import service_account
from dotenv import load_dotenv

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


