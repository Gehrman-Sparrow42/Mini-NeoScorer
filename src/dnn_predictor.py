"""
Mini-NeoScorer: Deep Neural Network (DNN) HLA Predictor Module.
Provides modular, production-grade peptide-MHC binding affinity and antigen presentation
scoring using:
1. MHCflurry 2.0 (Local feedforward deep neural network with presentation/processing models)
2. NetMHCpan-4.1 (Official NIH IEDB REST API with persistent local caching)
3. Ensemble (Consensus model combining MHCflurry and NetMHCpan for high-fidelity clinical selection)
4. PWM (Ultra-fast position weight matrix baseline)
"""

from abc import ABC, abstractmethod
import gc
import json
import math
import os
from pathlib import Path
import time
from typing import Dict, List, Optional, Tuple
import requests


# Default persistence path for NetMHCpan API cache
CACHE_DIR = Path(__file__).resolve().parent.parent / "data"
NETMHCPAN_CACHE_FILE = CACHE_DIR / "netmhcpan_cache.json"


class BaseHLAPredictor(ABC):
    """Abstract base class for HLA binding affinity predictors."""

    @abstractmethod
    def predict_peptides(
        self, peptides: List[str], allele: str
    ) -> Dict[str, Dict[str, float]]:
        """
        Predict binding metrics for a list of peptides against an HLA allele.
        Returns mapping: peptide -> {
            'ic50': float (nM),
            'percentile_rank': float (0-100),
            'presentation_score': float (0-1),
            'source': str
        }
        """
        pass

    def release_memory(self) -> None:
        """Hook for explicit garbage collection and memory release."""
        pass


class PWMPredictor(BaseHLAPredictor):
    """Local Position Weight Matrix (PWM) scoring baseline for HLA-A*02:01."""

    def __init__(self):
        from src.hla_scorer import score_peptide_a0201, raw_score_to_ic50
        self.score_fn = score_peptide_a0201
        self.ic50_fn = raw_score_to_ic50

    def predict_peptides(
        self, peptides: List[str], allele: str
    ) -> Dict[str, Dict[str, float]]:
        results = {}
        for pep in peptides:
            raw = self.score_fn(pep)
            ic50 = self.ic50_fn(raw)
            # Estimate pseudo-percentile rank from IC50: <50nM ~ 0.5%, 500nM ~ 2.0%
            pseudo_rank = min(100.0, max(0.01, round((ic50 / 50000.0) * 100.0, 2)))
            results[pep] = {
                "ic50": ic50,
                "percentile_rank": pseudo_rank,
                "presentation_score": round(max(0.0, 1.0 - (ic50 / 1000.0)), 3),
                "source": "PWM_Baseline",
            }
        return results


class MHCflurryPredictor(BaseHLAPredictor):
    """
    Local deep neural network predictor using MHCflurry 2.0 (PyTorch backend).
    Predicts both groove binding affinity (IC50) and antigen presentation score
    (integrating proteasomal cleavage and TAP transport probability).
    """

    def __init__(self):
        self.predictor = None

    def _ensure_loaded(self):
        if self.predictor is None:
            from mhcflurry import Class1PresentationPredictor
            self.predictor = Class1PresentationPredictor.load()

    def predict_peptides(
        self, peptides: List[str], allele: str
    ) -> Dict[str, Dict[str, float]]:
        self._ensure_loaded()
        unique_peps = list(dict.fromkeys(peptides))
        if not unique_peps:
            return {}

        # MHCflurry format: HLA-A*02:01 -> HLA-A*02:01
        df = self.predictor.predict(
            peptides=unique_peps,
            alleles=[allele],
            verbose=0,
        )

        results = {}
        for _, row in df.iterrows():
            pep = str(row["peptide"])
            ic50 = float(row.get("affinity", 50000.0))
            pres_score = float(row.get("presentation_score", 0.0))
            pres_rank = float(row.get("presentation_percentile", 100.0))
            results[pep] = {
                "ic50": round(ic50, 2),
                "percentile_rank": round(pres_rank, 2),
                "presentation_score": round(pres_score, 4),
                "source": "MHCflurry_2.0",
            }
        return results

    def release_memory(self) -> None:
        """Purge model weights from RAM and force garbage collection."""
        if self.predictor is not None:
            del self.predictor
            self.predictor = None
            gc.collect()


class NetMHCpanIEDBPredictor(BaseHLAPredictor):
    """
    Global gold-standard neural network predictor via official NIH IEDB REST API.
    Uses NetMHCpan-4.1 binding affinity (BA) algorithm with persistent local caching.
    """

    API_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"

    def __init__(self, cache_file: Path = NETMHCPAN_CACHE_FILE):
        self.cache_file = cache_file
        self.cache: Dict[str, Dict[str, float]] = {}
        self._load_cache()

    def _load_cache(self) -> None:
        if self.cache_file.exists():
            try:
                with open(self.cache_file, "r", encoding="utf-8") as f:
                    self.cache = json.load(f)
            except Exception:
                self.cache = {}

    def _save_cache(self) -> None:
        try:
            self.cache_file.parent.mkdir(parents=True, exist_ok=True)
            with open(self.cache_file, "w", encoding="utf-8") as f:
                json.dump(self.cache, f, indent=2)
        except Exception:
            pass

    def predict_peptides(
        self, peptides: List[str], allele: str
    ) -> Dict[str, Dict[str, float]]:
        unique_peps = list(dict.fromkeys(peptides))
        results: Dict[str, Dict[str, float]] = {}
        missing_peps: List[str] = []

        # Check local cache first
        for pep in unique_peps:
            cache_key = f"{allele}_{pep}"
            if cache_key in self.cache:
                results[pep] = self.cache[cache_key]
            else:
                missing_peps.append(pep)

        # Query NIH IEDB grouped by peptide length
        if missing_peps:
            from collections import defaultdict
            by_len = defaultdict(list)
            for p in missing_peps:
                by_len[len(p)].append(p)

            chunk_size = 150
            for kmer_len, peps_of_len in by_len.items():
                for i in range(0, len(peps_of_len), chunk_size):
                    chunk = peps_of_len[i : i + chunk_size]
                    fasta_payload = "\n".join(
                        [f">p{idx}\n{p}" for idx, p in enumerate(chunk)]
                    )
                    try:
                        response = requests.post(
                            self.API_URL,
                            data={
                                "method": "netmhcpan_ba",
                                "sequence_text": fasta_payload,
                                "allele": allele,
                                "length": str(kmer_len),
                            },
                            timeout=45,
                        )
                        if response.status_code == 200:
                            lines = response.text.strip().splitlines()
                            header = lines[0].split("\t") if lines else []
                            if "peptide" in header and "ic50" in header:
                                p_idx = header.index("peptide")
                                ic_idx = header.index("ic50")
                                rank_idx = (
                                    header.index("percentile_rank")
                                    if "percentile_rank" in header
                                    else -1
                                )

                                for row in lines[1:]:
                                    cols = row.split("\t")
                                    if len(cols) > max(p_idx, ic_idx):
                                        pep_str = cols[p_idx]
                                        ic50_val = float(cols[ic_idx])
                                        rank_val = (
                                            float(cols[rank_idx])
                                            if rank_idx != -1 and len(cols) > rank_idx
                                            else 50.0
                                        )
                                        item = {
                                            "ic50": round(ic50_val, 2),
                                            "percentile_rank": round(rank_val, 2),
                                            "presentation_score": round(
                                                max(0.0, 1.0 - (ic50_val / 1000.0)), 4
                                            ),
                                            "source": "NetMHCpan_4.1_IEDB",
                                        }
                                        results[pep_str] = item
                                        self.cache[f"{allele}_{pep_str}"] = item
                    except Exception as e:
                        # Fallback to conservative estimate if remote network fails
                        for p in chunk:
                            if p not in results:
                                results[p] = {
                                    "ic50": 50000.0,
                                    "percentile_rank": 50.0,
                                    "presentation_score": 0.0,
                                    "source": "NetMHCpan_Unavailable_Fallback",
                                }

            self._save_cache()

        return results


class EnsemblePredictor(BaseHLAPredictor):
    """
    Consensus Ensemble Predictor.
    Combines MHCflurry 2.0 and NetMHCpan-4.1:
    - Consensus IC50: Geometric mean = sqrt(IC50_mhcflurry * IC50_netmhcpan)
    - Consensus Percentile Rank: Arithmetic mean = (Rank_mhcflurry + Rank_netmhcpan) / 2
    - Presentation Score: MHCflurry presentation probability
    Suppresses single-model false positives for clinical trial safety.
    """

    def __init__(self):
        self.mhcflurry = MHCflurryPredictor()
        self.netmhcpan = NetMHCpanIEDBPredictor()

    def predict_peptides(
        self, peptides: List[str], allele: str
    ) -> Dict[str, Dict[str, float]]:
        # Run predictions across both models
        mf_res = self.mhcflurry.predict_peptides(peptides, allele)
        net_res = self.netmhcpan.predict_peptides(peptides, allele)

        ensemble_results = {}
        for pep in peptides:
            mf = mf_res.get(
                pep,
                {"ic50": 50000.0, "percentile_rank": 50.0, "presentation_score": 0.0},
            )
            net = net_res.get(
                pep,
                {"ic50": 50000.0, "percentile_rank": 50.0, "presentation_score": 0.0},
            )

            # Geometric mean of binding affinities
            geom_ic50 = math.sqrt(max(0.01, mf["ic50"]) * max(0.01, net["ic50"]))
            # Arithmetic mean of percentile ranks
            mean_rank = (mf["percentile_rank"] + net["percentile_rank"]) / 2.0
            # Antigen processing & presentation score from MHCflurry
            pres_score = mf.get("presentation_score", 0.0)

            ensemble_results[pep] = {
                "ic50": round(geom_ic50, 2),
                "percentile_rank": round(mean_rank, 2),
                "presentation_score": round(pres_score, 4),
                "mhcflurry_ic50": mf["ic50"],
                "mhcflurry_rank": mf["percentile_rank"],
                "netmhcpan_ic50": net["ic50"],
                "netmhcpan_rank": net["percentile_rank"],
                "source": "Consensus_Ensemble(MHCflurry+NetMHCpan)",
            }

        return ensemble_results

    def release_memory(self) -> None:
        self.mhcflurry.release_memory()
        self.netmhcpan.release_memory()


def get_predictor(name: str = "ensemble") -> BaseHLAPredictor:
    """Factory function to instantiate the selected HLA predictor."""
    lowered = name.lower()
    if lowered == "mhcflurry":
        return MHCflurryPredictor()
    elif lowered in ("netmhcpan", "iedb"):
        return NetMHCpanIEDBPredictor()
    elif lowered == "pwm":
        return PWMPredictor()
    elif lowered in ("ensemble", "consensus"):
        return EnsemblePredictor()
    else:
        raise ValueError(
            f"Unknown predictor '{name}'. Choose from: 'ensemble', 'mhcflurry', 'netmhcpan', 'pwm'."
        )
