import hashlib
import json
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
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
NCBI_ELINK = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
UNIPROTKB_SEARCH = "https://rest.uniprot.org/uniprotkb/search"

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


# Input modes shown in the "Input type" selector
FASTA_MODE = "FASTA Protein Sequence"
LOCUS_MODE = "Locus ID"

# Locus IDs are organism-specific, so the organism must be known for a lookup.
ARABIDOPSIS = "Arabidopsis thaliana"
ARABIDOPSIS_TAXID = 3702
CUSTOM_ORGANISM = "Other / custom organism"
ORGANISMS = [ARABIDOPSIS, CUSTOM_ORGANISM]

# Arabidopsis (TAIR/Araport) locus format, e.g. AT1G01010 (chromosomes 1-5, C = chloroplast, M = mitochondrion)
ARABIDOPSIS_LOCUS_RE = re.compile(r"^AT[1-5CM]G\d{5}$")
# Very permissive format check for other organisms (letters, digits and _ . - : only)
GENERIC_LOCUS_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.\-:]{1,63}$")


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
def fetch_and_show_interpro(accession):
    """Query InterPro with a (base) UniProt accession and display the JSON."""
    r = requests.get(INTERPRO_ENTRY.format(accession), timeout=TIMEOUT)
    if r.status_code in (204, 404):
        raise LookupFailed(f"InterPro has no entries for {accession}.")
    r.raise_for_status()
    st.json(r.json())


def query_interpro(sequence):
    accessions = find_uniprot_accessions(sequence)
    accession = accessions[0].split("-")[0]  # InterPro wants the base accession, not an isoform
    st.success(f"UniProt accession used for InterPro: **{accession}**")

    fetch_and_show_interpro(accession)


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


# ------------------------------------------------------------ Locus ID mode --
# Workflow:  Locus ID -> NCBI Gene record -> RefSeq protein accession(s)
#                     -> UniProt accession (via GeneID cross-reference, else gene-name search)
#                     -> InterPro (using the UniProt accession)
# A Locus ID is NEVER sent directly to UniProt/InterPro as if it were a protein accession.
def normalize_locus(locus, organism):
    """Return (clean_locus, error, note). Validates the format for the chosen organism."""
    locus = re.sub(r"\s+", "", locus or "")
    if not locus:
        return "", "Locus ID is empty. Enter a Locus ID, for example At1g01010.", None
    if not organism:
        return locus, "Organism is empty. Enter the scientific name of the organism (e.g. Oryza sativa).", None

    if organism == ARABIDOPSIS:
        note = None
        # "AT1G01010.1" is a transcript/protein model; the locus is the part before the dot
        if re.match(r"^AT[1-5CM]G\d{5}\.\d+$", locus, re.IGNORECASE):
            locus, note = locus.split(".")[0], "Isoform suffix removed: a Locus ID identifies the gene, not a single transcript."
        locus = locus.upper()
        if not ARABIDOPSIS_LOCUS_RE.match(locus):
            return locus, (
                f"'{locus}' is not a valid Arabidopsis thaliana Locus ID. "
                "Expected the form AT#G##### (for example At1g01010), where # after AT is 1-5, C or M."
            ), None
        return locus, None, note

    if not GENERIC_LOCUS_RE.match(locus):
        return locus, (
            f"'{locus}' does not look like a valid Locus ID. Use only letters, digits and _ . - : "
            "(2-64 characters)."
        ), None
    return locus, None, None


def ncbi_get(url, params, email, api_key):
    """GET an NCBI E-utilities endpoint and return parsed JSON (email / API key are optional)."""
    full = {"tool": TOOL_NAME, "retmode": "json", **params}
    if email:
        full["email"] = email
    if api_key:
        full["api_key"] = api_key
    # NCBI allows ~3 requests/s without an API key and ~10/s with one
    time.sleep(0.12 if api_key else 0.35)
    r = requests.get(url, params=full, timeout=TIMEOUT)
    r.raise_for_status()
    return r.json()


def resolve_ncbi_gene(locus, organism, email, api_key):
    """Locus ID + organism -> NCBI Gene record. Returns a dict with the chosen gene and all candidates."""
    org = organism.replace('"', "")
    term = f'"{locus}"[All Fields] AND "{org}"[Organism] AND alive[prop]'
    found = ncbi_get(NCBI_ESEARCH, {"db": "gene", "term": term, "retmax": 20}, email, api_key)
    ids = found.get("esearchresult", {}).get("idlist", [])
    if not ids:
        raise LookupFailed(
            f"No NCBI Gene record found for Locus ID '{locus}' in {organism}. "
            "Check the spelling and that the selected organism is correct."
        )

    summary = ncbi_get(NCBI_ESUMMARY, {"db": "gene", "id": ",".join(ids)}, email, api_key).get("result", {})
    docs = [summary[uid] for uid in summary.get("uids", []) if uid in summary]
    if not docs:
        raise LookupFailed(f"NCBI found gene IDs {', '.join(ids)} but returned no summary for them.")

    # Prefer records that literally mention the locus (name, alias, locus tag ...)
    exact = [d for d in docs if locus.lower() in json.dumps(d).lower()]
    pool = exact or docs
    return {
        "gene_id": str(pool[0]["uid"]),
        "summary": pool[0],
        "candidates": [
            {"gene_id": str(d["uid"]), "name": d.get("name"), "description": d.get("description")} for d in pool
        ],
    }


def resolve_ncbi_proteins(gene_id, email, api_key):
    """NCBI Gene ID -> linked protein records (RefSeq first, NP_ before XP_)."""
    prot_ids = []
    for linkname in ("gene_protein_refseq", "gene_protein"):  # fall back to all linked proteins
        data = ncbi_get(NCBI_ELINK, {"dbfrom": "gene", "db": "protein", "id": gene_id, "linkname": linkname},
                        email, api_key)
        for linkset in data.get("linksets", []):
            for db in linkset.get("linksetdbs", []):
                prot_ids += [str(x) for x in db.get("links", [])]
        if prot_ids:
            break
    if not prot_ids:
        raise LookupFailed(
            f"NCBI Gene {gene_id} has no linked protein record (it may be a non-coding gene or pseudogene)."
        )

    summary = ncbi_get(NCBI_ESUMMARY, {"db": "protein", "id": ",".join(prot_ids[:20])}, email, api_key).get("result", {})
    proteins = [
        {"uid": uid, "accession": summary[uid].get("accessionversion") or summary[uid].get("caption"),
         "title": summary[uid].get("title")}
        for uid in summary.get("uids", []) if uid in summary
    ]
    proteins = [p for p in proteins if p["accession"]]
    if not proteins:
        raise LookupFailed(f"NCBI Gene {gene_id} links to protein records, but no accessions could be read.")
    proteins.sort(key=lambda p: (not p["accession"].startswith("NP_"), p["accession"]))
    return proteins


def uniprot_search(query):
    """Search UniProtKB; return matching entries with reviewed (Swiss-Prot) entries first."""
    r = requests.get(
        UNIPROTKB_SEARCH,
        params={"query": query, "format": "json", "size": 25,
                "fields": "accession,reviewed,protein_name,gene_names,organism_name"},
        headers={"Accept": "application/json"},
        timeout=TIMEOUT,
    )
    r.raise_for_status()
    results = r.json().get("results", [])
    return sorted(results, key=lambda e: (not e.get("entryType", "").startswith("UniProtKB reviewed"),
                                          e.get("primaryAccession", "")))


def find_uniprot_for_locus(locus, organism, gene_id):
    """Find UniProtKB entries for a locus. Returns (accessions, method_description)."""
    org = organism.replace('"', "")
    org_filter = f"organism_id:{ARABIDOPSIS_TAXID}" if organism == ARABIDOPSIS else f'organism_name:"{org}"'
    strategies = []
    if gene_id:  # most reliable: UniProt cross-references NCBI GeneID
        strategies.append((f"NCBI GeneID cross-reference ({gene_id})", f"xref:geneid-{gene_id}"))
    strategies.append((f"gene name / ordered locus name '{locus}' in {organism}",
                       f'(gene:"{locus}") AND ({org_filter})'))
    for label, query in strategies:
        results = uniprot_search(query)
        if results:
            return [e["primaryAccession"] for e in results], label
    raise LookupFailed(
        f"No UniProt entry found for Locus ID '{locus}' in {organism}. "
        "The gene may have no UniProt record, or the organism selection may be wrong."
    )


def resolve_locus_to_uniprot(locus, organism, email, api_key):
    """Locus ID -> UniProt accession (via NCBI Gene when possible). Shows what was resolved."""
    gene = None
    try:
        gene = resolve_ncbi_gene(locus, organism, email, api_key)
        warn_if_ambiguous_gene(gene)
    except LookupFailed as e:
        st.info(f"{e} Trying UniProt directly with the gene name instead.")

    accessions, method = find_uniprot_for_locus(locus, organism, gene["gene_id"] if gene else None)
    accession = accessions[0]
    show_locus_identifiers(locus, organism, gene=gene, uniprot=accession)
    st.caption(f"UniProt entry found via {method}.")
    if len(accessions) > 1:
        st.warning("Multiple UniProt entries match this Locus ID; using the first "
                   f"(reviewed entries are ranked first). Others: {', '.join(accessions[1:])}")
    return accession


def warn_if_ambiguous_gene(gene):
    if len(gene["candidates"]) > 1:
        others = "; ".join(f"{c['gene_id']} ({c['name']}: {c['description']})" for c in gene["candidates"][1:])
        st.warning(f"Multiple NCBI Gene records match this Locus ID. Using the first (Gene ID {gene['gene_id']}). "
                   f"Other matches: {others}")


def show_locus_identifiers(locus, organism, gene=None, proteins=None, uniprot=None):
    """Display the identifiers discovered during the lookup."""
    lines = ["**Input type:** Locus ID", f"**Locus ID:** {locus}", f"**Organism:** {organism}"]
    if gene:
        s = gene["summary"]
        lines.append(f"**Resolved NCBI Gene ID:** {gene['gene_id']} ({s.get('name')} - {s.get('description')})")
    if proteins:
        lines.append("**Resolved NCBI protein accession(s):** " + ", ".join(p["accession"] for p in proteins))
    if uniprot:
        lines.append(f"**Resolved UniProt accession:** {uniprot}")
    st.success("\n\n".join(lines))


def query_ncbi_locus(locus, organism, email, api_key):
    gene = resolve_ncbi_gene(locus, organism, email, api_key)
    warn_if_ambiguous_gene(gene)
    proteins = resolve_ncbi_proteins(gene["gene_id"], email, api_key)
    show_locus_identifiers(locus, organism, gene=gene, proteins=proteins)
    if len(proteins) > 1:
        st.info(f"{len(proteins)} protein records are linked to this gene (isoforms); "
                f"using {proteins[0]['accession']} for the protein lookup.")

    # Look the primary protein accession up in the NCBI protein database, as in FASTA mode
    protein_search = ncbi_get(NCBI_ESEARCH, {"db": "protein", "term": f"{proteins[0]['accession']}[Accession]"},
                              email, api_key)
    st.json({"gene": gene["summary"], "linked_proteins": proteins, "esearch": protein_search})


def query_uniprot_locus(locus, organism, email, api_key):
    accession = resolve_locus_to_uniprot(locus, organism, email, api_key)
    r = requests.get(UNIPROTKB_ENTRY.format(accession), headers={"Accept": "application/json"}, timeout=TIMEOUT)
    r.raise_for_status()
    st.json(r.json())


def query_interpro_locus(locus, organism, email, api_key):
    accession = resolve_locus_to_uniprot(locus, organism, email, api_key).split("-")[0]  # base accession
    st.info(f"Querying InterPro with UniProt accession {accession}.")
    fetch_and_show_interpro(accession)


# ---------------------------------------------------------------------- UI --
st.title("Bio API Tester")

input_type = st.selectbox("Input type", [FASTA_MODE, LOCUS_MODE])
is_locus = input_type == LOCUS_MODE

# Defaults so the query code below works in either mode
fasta_text, header, sequence, error, note = "", None, "", None, None
locus, organism = "", ""

if not is_locus:
    # ------------------------------------------------ FASTA input (unchanged) --
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
else:
    # ------------------------------------------------------- Locus ID input --
    organism_choice = st.selectbox("Organism", ORGANISMS)
    if organism_choice == CUSTOM_ORGANISM:
        organism = st.text_input("Custom organism (scientific name)", placeholder="e.g. Oryza sativa").strip()
    else:
        organism = organism_choice
    locus_raw = st.text_input("Locus ID", placeholder="At1g01010")
    locus_clean, locus_error, locus_note = normalize_locus(locus_raw, organism) if locus_raw.strip() else ("", None, None)
    if locus_raw.strip() and not locus_error:
        st.caption(f"Locus ID: {locus_clean}" + (f" - {locus_note}" if locus_note else ""))
    st.caption("Locus IDs are organism-specific: the ID is first resolved through NCBI Gene, "
               "then mapped to a UniProt accession. It is never used directly as a protein accession.")

email = st.text_input("NCBI email")
api_key = st.text_input("NCBI API key", type="password")


def inputs_ready():
    """Re-validate at click time and show a clear error if the input is unusable."""
    if is_locus:
        _, err, _ = normalize_locus(locus_raw, organism)
        if err:
            st.error(err)
            return False
        return True
    _, _, err, _ = parse_fasta(fasta_text)
    if err:
        st.error(err)
        return False
    return True


if st.button("Query NCBI") and inputs_ready():
    if is_locus:
        with st.spinner("Resolving Locus ID through NCBI Gene and Protein..."):
            guarded("NCBI", lambda: query_ncbi_locus(locus_clean, organism, email, api_key))
    else:
        with st.spinner("Running NCBI BLAST - this can take a minute or more..."):
            guarded("NCBI", lambda: query_ncbi(sequence, email, api_key))

if st.button("Query UniProt") and inputs_ready():
    if is_locus:
        with st.spinner("Resolving Locus ID, then searching UniProt..."):
            guarded("UniProt", lambda: query_uniprot_locus(locus_clean, organism, email, api_key))
    else:
        with st.spinner("Searching UniProt..."):
            guarded("UniProt", lambda: query_uniprot(sequence))

if st.button("Query InterPro") and inputs_ready():
    if is_locus:
        with st.spinner("Resolving Locus ID to a UniProt accession, then querying InterPro..."):
            guarded("InterPro", lambda: query_interpro_locus(locus_clean, organism, email, api_key))
    else:
        with st.spinner("Finding UniProt accession, then querying InterPro..."):
            guarded("InterPro", lambda: query_interpro(sequence))
