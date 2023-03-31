## Efetch The Python 

<img src="https://huddle.eurostarsoftwaretesting.com/wp-content/uploads/2014/05/python-programming.jpg"  width="20%" height="20%">

This python script, which can be found [here](EfetchThePython.py), facilitates automated GenBank searches and returns a tab-delimited list of hits and optionally a FASTA file with the hits if the flag `--FASTA` is set. The parameter `--term` takes as arguments the typical search code of GenBank searches, e.g. `"Merops[Organism], COI"`. The parameter `--output` defines the prefix for the output files.

Since efetch requires an authentification, you have to provide your NCBI email ID via `--email` and your API key via `--api_key`. See [here](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) for more details how to obtain an API key.

The script uses the Biopython module which can be installed as follows:

```bash
pip install biopython
```

A typical commandline looks like this:

```bash
python EfetchThePython.py \
    --email <your mail address> \
    --api_key <your API key> \
    --Term "Merops[Organism], COI" \
    --Output data/Merops \
    --FASTA
```