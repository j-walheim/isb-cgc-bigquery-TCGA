# %%
import os
import pandas as pd
from google.cloud import bigquery
from google.oauth2 import service_account
from dotenv import load_dotenv

# Load the environment variables from the .env file
load_dotenv()

gene = "PIK3CA"
disease = 'LUAD'


credentials = service_account.Credentials.from_service_account_file(
    os.environ['GOOGLE_APPLICATION_CREDENTIALS'],
    scopes=["https://www.googleapis.com/auth/cloud-platform"],
)
client = bigquery.Client(credentials=credentials, project=credentials.project_id)


query = f"""
    SELECT *
    FROM `isb-cgc-bq.TCGA.masked_somatic_mutation_hg38_gdc_current`
    WHERE Hugo_Symbol = '{gene}' LIMIT 5
    """



# %%

test = client.query(query).to_dataframe()


