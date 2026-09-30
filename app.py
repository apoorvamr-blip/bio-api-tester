import hashlib
import re
import time

import requests
import streamlit as st

# ---------------------------------------------------------------- settings --
TIMEOUT = 30  # seconds, for every ordinary HTTP request

UNIPARC_SEARCH = "https://rest.uniprot.org/uniparc/search"
UNIPROTKB_ENTRY = "https://rest.uniprot.org/uniprotkb/{}"
INTERPRO_ENTRY = "https://www.ebi.ac.uk/interpro/api/entry/InterPro/protein/uniprot/{}/"
NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_BLAST = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

BLAST_DB = "nr"          # use "swissprot" or "refseq_protein" for faster (narrower) searches
BLAST_HITS = 5
BLAST_POLL_SECONDS = 15  # gap between BLAST status checks
BLAST_MAX_WAIT = 300     # give up after 5 minutes
TOOL_NAME = "bio-api-tester"

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY" "BXZUOJ")  # 20 standard + ambiguous/rare codes

EXAMPLE_FASTA = (
    ">sp|P69905|HBA_HUMAN Hemoglobin subunit alpha\n"
    "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG\n"
    "KKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTP\n"
    "AVHASLDKFLASVSTVLTSKYR"
)


class LookupFailed(Exception):
    """Raised when an API works but gives no usable match."""


# ----------------------------------------------------------------- parsing --
def parse_fasta(text):
    """Return (header, sequence, error, note). Accepts FASTA or a bare sequence."""
    lines = [line.strip() for line in text.strip().splitlines() if line.strip()]
    if not lines:
        return None, "", "The sequence is empty. Paste a FASTA record or a raw amino-acid sequence.", None

    header, seq_lines, records = None, [], 0
    for line in lines:
        if line.startswith(">"):
            records += 1
            if records > 1:
                break  # only the first record is used
            header = line[1:].strip()
        else:
            seq_lines.append(line)

    note = "Multiple FASTA records found; only the first one is used." if records > 1 else None

    sequence = re.sub(r"\s+", "", "".join(seq_lines)).upper()
    sequence = sequence.rstrip("*")  # some tools end protein FASTA with a stop symbol

    if not sequence:
        return header, "", "The FASTA header has no sequence after it.", note

    invalid = sorted(set(sequence) - VALID_AA)
    if invalid:
        return header, sequence, (
            "Invalid characters in the sequence: " + ", ".join(invalid)
            + ". Only amino-acid letters are allowed (ACDEFGHIKLMNPQRSTVWY, plus B, X, Z, U, O, J); "
            "remove digits, gaps and other symbols."
        ), note

    return header, sequence, None, note


# ------------------------------------------------------------ error output --
def show_request_error(service, e):
    if isinstance(e, requests.Timeout):
        st.error(f"{service} request timed out. Please try again in a moment.")
    elif isinstance(e, requests.ConnectionError):
        st.error(f"Could not connect to {service}. Check your internet connection. ({e})")
    elif isinstance(e, requests.HTTPError):
        st.error(f"{service} returned HTTP {e.response.status_code}: {e}")
    else:
        st.error(f"{service} request failed: {e}")


def guarded(service, func):
    """Run func(); turn expected failures into friendly Streamlit messages."""
    try:
        func()
    except LookupFailed as e:
        st.warning(str(e))
    except (requests.RequestException, ValueError) as e:
        show_request_error(service, e)


# ------------------------------------------------------------------ UniProt --
def find_uniprot_accessions(sequence):
    """Exact-sequence lookup via UniParc (UniProt's sequence archive) using the MD5 checksum.

    Returns a list of UniProtKB accessions, best (reviewed Swiss-Prot) first.
    """
    checksum = hashlib.md5(sequence.encode()).hexdigest().upper()
    r = requests.get(
        UNIPARC_SEARCH,
        params={"query": f"checksum:{checksum}", "format": "json", "size": 10},
        headers={"Accept": "application/json"},
        timeout=TIMEOUT,
    )
    r.raise_for_status()

    rank = {"UniProtKB/Swiss-Prot": 0, "UniProtKB/TrEMBL": 1}
    found = {}
    for entry in r.json().get("results", []):
        for xref in entry.get("uniParcCrossReferences", []):
            db = xref.get("database", "")
            if xref.get("active") and db.startswith("UniProtKB") and xref.get("id"):
                found[xref["id"]] = rank.get(db, 2)  # 2 = isoforms etc.

    if not found:
        raise LookupFailed(
            "No UniProt protein has exactly this sequence. "
            "The lookup needs an exact match; a variant, fragment or novel sequence will not be found."
        )
    return sorted(found, key=lambda acc: (found[acc], acc))


def query_uniprot(sequence):
    accessions = find_uniprot_accessions(sequence)
    accession = accessions[0]
    st.success(f"UniProt accession found: **{accession}**")
    if len(accessions) > 1:
        st.info("Other UniProt entries with the same sequence: " + ", ".join(accessions[1:]))

    r = requests.get(UNIPROTKB_ENTRY.format(accession), headers={"Accept": "application/json"}, timeout=TIMEOUT)
    r.raise_for_status()
    st.json(r.json())


# ---------------------------------------------------------------- InterPro --
def query_interpro(sequence):
    accessions = find_uniprot_accessions(sequence)
    accession = accessions[0].split("-")[0]  # InterPro wants the base accession, not an isoform
    st.success(f"UniProt accession used for InterPro: **{accession}**")

    r = requests.get(INTERPRO_ENTRY.format(accession), timeout=TIMEOUT)
    if r.status_code in (204, 404):
        raise LookupFailed(f"InterPro has no entries for {accession}.")
    r.raise_for_status()
    st.json(r.json())


# --------------------------------------------------------------------- NCBI --
def blast_top_hits(sequence, email):
    """Run NCBI BLAST (blastp) via the URL API and return a short list of top hits."""
    common = {"tool": TOOL_NAME}
    if email:
        common["email"] = email

    put = requests.post(
        NCBI_BLAST,
        data={"CMD": "Put", "PROGRAM": "blastp", "DATABASE": BLAST_DB,
              "QUERY": sequence, "HITLIST_SIZE": BLAST_HITS, **common},
        timeout=TIMEOUT,
    )
    put.raise_for_status()
    rid = re.search(r"RID = (\S+)", put.text)
    rtoe = re.search(r"RTOE = (\d+)", put.text)
    if not rid:
        raise LookupFailed("NCBI BLAST did not accept the sequence (no request ID returned).")
    rid = rid.group(1)

    start = time.time()
    time.sleep(min(int(rtoe.group(1)) if rtoe else BLAST_POLL_SECONDS, 30))
    while True:
        info = requests.get(NCBI_BLAST, params={"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid, **common},
                            timeout=TIMEOUT)
        info.raise_for_status()
        status = re.search(r"Status=(\w+)", info.text)
        status = status.group(1) if status else "UNKNOWN"
        if status == "READY":
            if "ThereAreHits=yes" not in info.text:
                raise LookupFailed("NCBI BLAST finished but found no matching proteins.")
            break
        if status in ("FAILED", "UNKNOWN"):
            raise LookupFailed(f"NCBI BLAST search {status.lower()} (request ID {rid}).")
        if time.time() - start > BLAST_MAX_WAIT:
            raise LookupFailed(f"NCBI BLAST is still running after {BLAST_MAX_WAIT} s (request ID {rid}). Try again later.")
        time.sleep(BLAST_POLL_SECONDS)

    res = requests.get(NCBI_BLAST, params={"CMD": "Get", "FORMAT_TYPE": "JSON2_S", "RID": rid,
                                           "HITLIST_SIZE": BLAST_HITS, **common}, timeout=TIMEOUT)
    res.raise_for_status()
    search = res.json()["BlastOutput2"][0]["report"]["results"]["search"]
    query_len = search.get("query_len") or len(sequence)

    hits = []
    for hit in search.get("hits", []):
        desc, hsp = hit["description"][0], hit["hsps"][0]
        hits.append({
            "accession": desc.get("accession"),
            "title": desc.get("title"),
            "organism": desc.get("sciname"),
            "identity_percent": round(100 * hsp["identity"] / hsp["align_len"], 1),
            "query_coverage_percent": round(100 * (hsp["query_to"] - hsp["query_from"] + 1) / query_len, 1),
            "evalue": hsp["evalue"],
        })
    if not hits:
        raise LookupFailed("NCBI BLAST returned no hits.")
    return hits


def query_ncbi(sequence, email, api_key):
    hits = blast_top_hits(sequence, email)
    accession = hits[0]["accession"]
    st.success(f"Top NCBI protein hit: **{accession}** "
               f"({hits[0]['identity_percent']}% identity, {hits[0]['query_coverage_percent']}% query coverage)")

    # Look the accession up in the NCBI protein database (email / API key are optional)
    params = {"db": "protein", "term": f"{accession}[Accession]", "retmode": "json", "tool": TOOL_NAME}
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    r = requests.get(NCBI_ESEARCH, params=params, timeout=TIMEOUT)
    r.raise_for_status()
    st.json({"blast_top_hits": hits, "esearch": r.json()})


# ---------------------------------------------------------------------- UI --
st.title("Bio API Tester")

if "fasta_input" not in st.session_state:
    st.session_state["fasta_input"] = ""


def load_example():
    st.session_state["fasta_input"] = EXAMPLE_FASTA


fasta_text = st.text_area(
    "FASTA Protein Sequence",
    key="fasta_input",
    height=200,
    placeholder=">my_protein\nMKTIIALSYIFCLVFADYKDDDDK...",
)
st.button("Load example FASTA", on_click=load_example)

header, sequence, error, note = parse_fasta(fasta_text) if fasta_text.strip() else (None, "", None, None)
if fasta_text.strip() and not error:
    label = f" ({header})" if header else ""
    st.caption(f"Parsed sequence length: {len(sequence)} amino acids{label}")
    if note:
        st.caption(note)

email = st.text_input("NCBI email")
api_key = st.text_input("NCBI API key", type="password")


def sequence_ready():
    """Parse again at click time and show a clear error if the input is unusable."""
    _, _, err, _ = parse_fasta(fasta_text)
    if err:
        st.error(err)
        return False
    return True


if st.button("Query NCBI") and sequence_ready():
    with st.spinner("Running NCBI BLAST - this can take a minute or more..."):
        guarded("NCBI", lambda: query_ncbi(sequence, email, api_key))

if st.button("Query UniProt") and sequence_ready():
    with st.spinner("Searching UniProt..."):
        guarded("UniProt", lambda: query_uniprot(sequence))

if st.button("Query InterPro") and sequence_ready():
    with st.spinner("Finding UniProt accession, then querying InterPro..."):
        guarded("InterPro", lambda: query_interpro(sequence))
