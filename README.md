# bio-api-tester

A minimal Streamlit app that takes a **FASTA protein sequence** and queries three biological APIs live:

- **NCBI** (BLAST URL API + E-utilities esearch, protein database)
- **UniProt** REST API (UniParc sequence lookup, then UniProtKB record)
- **InterPro** API

The returned JSON is displayed directly in the app.

All data is fetched live from the APIs. No local biological data (CSV, Excel, JSON, FASTA, etc.) is used, and the pasted sequence is never written to disk.

## Install

```bash
pip install -r requirements.txt
```

## Run

```bash
streamlit run app.py
```

## Input

Paste either:

- a full FASTA record, e.g.

  ```
  >my_protein
  MKTIIALSYIFCLVFADYKDDDDK
  ```

- or just the amino-acid sequence, with no header.

The app removes whitespace and line breaks, drops the header (if any), upper-cases the sequence and checks that it contains only valid amino-acid letters. If several FASTA records are pasted, only the first is used. The parsed sequence length is shown before you query. Use **Load example FASTA** to fill in human hemoglobin subunit alpha.

## How each button works

| Button | Workflow |
| --- | --- |
| **Query UniProt** | Computes the sequence's MD5 checksum and searches UniParc (UniProt's sequence archive) for an exact match, picks the matching UniProtKB accession (reviewed Swiss-Prot entries first), then fetches that record from UniProtKB. The accession and full JSON are shown. |
| **Query InterPro** | Finds the UniProt accession the same way, then queries InterPro with that accession. |
| **Query NCBI** | Runs `blastp` against NCBI's `nr` database through the BLAST URL API, takes the top-hit accession, then looks it up in the protein database with esearch. The top hits (accession, identity, coverage, E-value) and the esearch JSON are shown. |

Notes:

- UniProt and InterPro lookups need an **exact** sequence match. Variants, fragments or novel sequences give a "no match" message.
- NCBI BLAST is asynchronous and can take from several seconds to a few minutes. To speed it up, change `BLAST_DB` at the top of `app.py` (for example to `swissprot` or `refseq_protein`).
- If the same sequence exists in several UniProt entries, the first one is used and the others are listed.

## NCBI email and API key

Both are optional; enter your own in the app. The email is sent to BLAST and E-utilities, and the API key to E-utilities. No API key is stored in the code.

## Example

Click **Load example FASTA** (human hemoglobin subunit alpha, UniProt `P69905`) and press any of the query buttons.
