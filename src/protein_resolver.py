"""
Mini-NeoScorer: Automated Human Reference Protein Resolver.
Dynamically resolves, validates, and caches canonical full-length human protein sequences
from UniProtKB for any gene symbol, supporting automated 21-mer context window generation.
Includes multi-threaded parallel prefetching for whole-exome scale.
"""

import concurrent.futures
import json
import re
import ssl
import threading
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

BASE_DIR = Path(__file__).resolve().parent.parent
CACHE_FILE = BASE_DIR / "data" / "protein_cache.json"

AA_3_TO_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}


class ProteinResolver:
    def __init__(self, cache_file: Path = CACHE_FILE):
        self.cache_file = cache_file
        self._lock = threading.Lock()
        self.cache: Dict[str, str] = self._load_cache()
        self._ssl_ctx = ssl.create_default_context()
        self._ssl_ctx.check_hostname = False
        self._ssl_ctx.verify_mode = ssl.CERT_NONE

    def _load_cache(self) -> Dict[str, str]:
        if self.cache_file.exists():
            try:
                with open(self.cache_file, "r", encoding="utf-8") as f:
                    return json.load(f)
            except Exception:
                return {}
        return {}

    def _save_cache(self) -> None:
        with self._lock:
            self.cache_file.parent.mkdir(parents=True, exist_ok=True)
            with open(self.cache_file, "w", encoding="utf-8") as f:
                json.dump(self.cache, f, indent=2)

    def fetch_protein_sequence(self, gene_symbol: str) -> Optional[str]:
        """
        Query UniProtKB for the canonical human reviewed (Swiss-Prot) sequence.
        """
        with self._lock:
            if gene_symbol in self.cache:
                return self.cache[gene_symbol]

        query = f"gene_exact:{gene_symbol} AND organism_id:9606 AND reviewed:true"
        url = f"https://rest.uniprot.org/uniprotkb/search?query={urllib.parse.quote(query)}&format=fasta&size=1"

        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "MiniNeoScorer/2.0"}
            )
            with urllib.request.urlopen(req, context=self._ssl_ctx, timeout=8) as resp:
                data = resp.read().decode("utf-8").strip()
                if not data or not data.startswith(">"):
                    return None
                lines = data.split("\n")
                sequence = "".join(lines[1:]).replace(" ", "").strip()
                if sequence:
                    with self._lock:
                        self.cache[gene_symbol] = sequence
                    self._save_cache()
                    return sequence
        except Exception:
            return None

        return None

    def prefetch_proteins(self, gene_symbols: List[str], max_workers: int = 16) -> int:
        """
        High-throughput multi-threaded prefetching of uncached protein sequences.
        """
        needed = [g for g in set(gene_symbols) if g and g not in self.cache]
        if not needed:
            return 0

        print(f"  -> Prefetching {len(needed)} uncached human protein sequences across {max_workers} threads...")
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            list(executor.map(self.fetch_protein_sequence, needed))

        self._save_cache()
        return len(needed)

    @staticmethod
    def parse_protein_change(change_str: str) -> Optional[Tuple[str, int, str]]:
        if not change_str or change_str in (".", "UNKNOWN"):
            return None

        clean = change_str.replace("p.", "")

        match_1 = re.match(r"^([A-Z])(\d+)([A-Z])$", clean)
        if match_1:
            wt, pos, mt = match_1.groups()
            return wt, int(pos), mt

        match_3 = re.match(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", clean)
        if match_3:
            wt3, pos, mt3 = match_3.groups()
            if wt3 in AA_3_TO_1 and mt3 in AA_3_TO_1:
                return AA_3_TO_1[wt3], int(pos), AA_3_TO_1[mt3]

        return None

    def extract_context_window(
        self, gene_symbol: str, protein_change: str, flank: int = 10
    ) -> Optional[Dict]:
        parsed = self.parse_protein_change(protein_change)
        if not parsed:
            return None

        wt_aa, pos_1based, mt_aa = parsed
        protein_seq = self.fetch_protein_sequence(gene_symbol)
        if not protein_seq:
            return None

        pos_0based = pos_1based - 1
        if pos_0based < 0 or pos_0based >= len(protein_seq):
            return None

        if protein_seq[pos_0based] != wt_aa:
            return None

        start = max(0, pos_0based - flank)
        end = min(len(protein_seq), pos_0based + flank + 1)

        upstream = protein_seq[start:pos_0based]
        downstream = protein_seq[pos_0based + 1:end]

        return {
            "gene": gene_symbol,
            "protein_change": protein_change,
            "pos_1based": pos_1based,
            "pos_0based": pos_0based,
            "wt_aa": wt_aa,
            "mt_aa": mt_aa,
            "wt_window": upstream + wt_aa + downstream,
            "mt_window": upstream + mt_aa + downstream,
            "protein_length": len(protein_seq),
        }
