import requests
import pubchempy as pcp
def download_sdf(cid: int, save = True, path = "temp/") -> str:
    """Download the 3D structure of a compound in SDF format from PubChem and return the SDF text."""
    try:
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF'
        response = requests.get(url)
        response.raise_for_status()
        sdf =  response.text
    except (requests.exceptions.RequestException, pcp.PubChemHTTPError) as e:
        print(f"Error downloading SDF for CID {cid}: {e}")
        return False, None
    if sdf is None:
        print(f"empty SDF for CID {cid}")
        return False, None 
    #get Smiles of the compound
    download_path = f"{path}{cid}_pubchem.sdf"
    if save:
        with open(download_path, "w") as f:
            f.write(sdf)
    return True, download_path


if __name__ == "__main__":
    id = 5330790 
    download_sdf(id)
