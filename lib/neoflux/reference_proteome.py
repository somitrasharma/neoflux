"""
Reference proteome provider for peptide generation.

Maps a gene symbol to its canonical human (Swiss-Prot, reviewed) protein
sequence so that real variant calls produce real peptide windows instead of
falling back to a synthetic context. Sequences are fetched once from the
UniProt REST API and cached to disk; subsequent lookups are offline.

If a gene cannot be resolved (no network and no cache, or no reviewed human
entry), ``get_protein`` returns ``None`` and callers fall back to their
synthetic-context behavior, so the pipeline still runs without a network.
"""

from __future__ import annotations

import logging
import ssl
import subprocess
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Optional

logger = logging.getLogger("neoflux")

# Cache lives beside the bundled assets so it survives between runs.
_CACHE_DIR = Path(__file__).resolve().parent.parent.parent / "assets" / "proteome_cache"
_UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"

# In-process memo: gene symbol -> sequence (or None if a lookup already failed).
_MEMO: dict[str, Optional[str]] = {}


def _ssl_context() -> Optional[ssl.SSLContext]:
    try:
        import certifi

        return ssl.create_default_context(cafile=certifi.where())
    except Exception:
        return None


def _parse_fasta(text: str) -> Optional[str]:
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines or not lines[0].startswith(">"):
        return None
    seq = "".join(ln.strip() for ln in lines[1:] if not ln.startswith(">"))
    seq = "".join(c for c in seq if c.isalpha())
    return seq or None


def _fetch_uniprot(gene: str) -> Optional[str]:
    query = f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true"
    url = _UNIPROT_SEARCH + "?" + urllib.parse.urlencode(
        {"query": query, "format": "fasta", "size": "1"}
    )

    # Preferred path: urllib with a verified CA bundle.
    ctx = _ssl_context()
    if ctx is not None:
        try:
            with urllib.request.urlopen(url, timeout=20, context=ctx) as resp:
                return _parse_fasta(resp.read().decode())
        except Exception as exc:  # network/SSL/HTTP — fall through to curl
            logger.debug("urllib UniProt fetch failed for %s: %s", gene, exc)

    # Fallback: curl (handles environments where Python lacks a CA bundle).
    try:
        out = subprocess.run(
            ["curl", "-sS", "-m", "20", url],
            capture_output=True, text=True, timeout=25,
        )
        if out.returncode == 0:
            return _parse_fasta(out.stdout)
    except Exception as exc:
        logger.debug("curl UniProt fetch failed for %s: %s", gene, exc)
    return None


def get_protein(gene_symbol: Optional[str]) -> Optional[str]:
    """Return the canonical human protein sequence for a gene symbol, or None.

    Resolution order: in-process memo -> on-disk cache -> UniProt fetch.
    A successful fetch is written to the cache for future offline use.
    """
    if not gene_symbol:
        return None
    gene = gene_symbol.strip().upper()
    if gene in _MEMO:
        return _MEMO[gene]

    cache_file = _CACHE_DIR / f"{gene}.fasta"
    if cache_file.exists():
        seq = _parse_fasta(cache_file.read_text())
        _MEMO[gene] = seq
        return seq

    seq = _fetch_uniprot(gene)
    if seq:
        try:
            _CACHE_DIR.mkdir(parents=True, exist_ok=True)
            cache_file.write_text(f">{gene}|canonical|9606\n{seq}\n")
        except OSError:
            pass
    _MEMO[gene] = seq
    return seq
