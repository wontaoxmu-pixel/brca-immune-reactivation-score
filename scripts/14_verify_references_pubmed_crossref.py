#!/usr/bin/env python3
"""Verify manuscript PMID and DOI metadata using public APIs.

This script parses the References section of the working manuscript, checks
PMIDs against NCBI E-utilities, checks DOI records against CrossRef, and writes
CSV/Markdown audit outputs. It uses only the Python standard library.
"""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = PROJECT_ROOT / "manuscript" / "manuscript_draft_v0.1.md"
OUT_CSV = PROJECT_ROOT / "results" / "tables" / "reference_verification_report.csv"
OUT_MD = PROJECT_ROOT / "manuscript" / "reference_verification_2026-06-14.md"

USER_AGENT = "brca-immune-scrna-paper-reference-check/0.1 (local reproducibility audit)"


def normalize_text(value: str) -> str:
    value = re.sub(r"[^a-z0-9]+", " ", value.lower())
    return re.sub(r"\s+", " ", value).strip()


def normalize_doi(value: str) -> str:
    value = value.strip().rstrip(".")
    value = re.sub(r"^https?://(dx\.)?doi\.org/", "", value, flags=re.I)
    return value.lower()


def request_json(url: str, timeout: int = 20) -> tuple[str, object | None]:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as response:
            return "ok", json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        return f"http_{exc.code}", None
    except Exception as exc:  # noqa: BLE001 - report exact API failure in audit output
        return f"error:{exc.__class__.__name__}", None


def parse_references(text: str) -> list[dict[str, str]]:
    in_refs = False
    refs: list[dict[str, str]] = []
    for line in text.splitlines():
        if line.strip() == "## References":
            in_refs = True
            continue
        if in_refs and line.startswith("## "):
            break
        if not in_refs:
            continue
        match = re.match(r"^(\d+)\.\s+(.*?)\s+PMID:\s*(\d+)\.\s+DOI:\s*(\S+)\.?\s*$", line)
        if not match:
            continue
        ref_no, citation_text, pmid, doi = match.groups()
        title_match = re.search(r"et al\.\s+(.*?)\.\s+[A-Za-z]", citation_text)
        refs.append(
            {
                "ref_no": ref_no,
                "citation_text": citation_text,
                "manuscript_title_fragment": title_match.group(1) if title_match else citation_text,
                "pmid": pmid,
                "doi": normalize_doi(doi),
            }
        )
    return refs


def verify_pubmed(pmid: str) -> dict[str, str]:
    params = urllib.parse.urlencode({"db": "pubmed", "id": pmid, "retmode": "json"})
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?{params}"
    status, payload = request_json(url)
    result = {
        "pubmed_api_status": status,
        "pubmed_title": "",
        "pubmed_doi": "",
        "pubmed_pmid_found": "FALSE",
    }
    if status != "ok" or not isinstance(payload, dict):
        return result
    records = payload.get("result", {})
    record = records.get(pmid, {}) if isinstance(records, dict) else {}
    if isinstance(record, dict) and record:
        result["pubmed_pmid_found"] = "TRUE"
        result["pubmed_title"] = str(record.get("title", "")).rstrip(".")
        articleids = record.get("articleids", [])
        if isinstance(articleids, list):
            for item in articleids:
                if isinstance(item, dict) and str(item.get("idtype", "")).lower() == "doi":
                    result["pubmed_doi"] = normalize_doi(str(item.get("value", "")))
                    break
    return result


def verify_crossref(doi: str) -> dict[str, str]:
    encoded = urllib.parse.quote(doi, safe="")
    url = f"https://api.crossref.org/works/{encoded}"
    status, payload = request_json(url)
    result = {
        "crossref_api_status": status,
        "crossref_title": "",
        "crossref_doi": "",
        "crossref_doi_found": "FALSE",
    }
    if status != "ok" or not isinstance(payload, dict):
        return result
    message = payload.get("message", {})
    if isinstance(message, dict):
        result["crossref_doi_found"] = "TRUE"
        titles = message.get("title", [])
        if isinstance(titles, list) and titles:
            result["crossref_title"] = str(titles[0])
        result["crossref_doi"] = normalize_doi(str(message.get("DOI", "")))
    return result


def main() -> None:
    refs = parse_references(MANUSCRIPT.read_text(encoding="utf-8"))
    rows = []
    for ref in refs:
        pubmed = verify_pubmed(ref["pmid"])
        time.sleep(0.15)
        crossref = verify_crossref(ref["doi"])
        time.sleep(0.15)

        manuscript_title_norm = normalize_text(ref["manuscript_title_fragment"])
        pubmed_title_norm = normalize_text(pubmed["pubmed_title"])
        crossref_title_norm = normalize_text(crossref["crossref_title"])

        rows.append(
            {
                **ref,
                **pubmed,
                **crossref,
                "pmid_doi_match": str(
                    bool(pubmed["pubmed_doi"])
                    and normalize_doi(pubmed["pubmed_doi"]) == normalize_doi(ref["doi"])
                ).upper(),
                "crossref_doi_match": str(
                    bool(crossref["crossref_doi"])
                    and normalize_doi(crossref["crossref_doi"]) == normalize_doi(ref["doi"])
                ).upper(),
                "pubmed_title_contains_manuscript_fragment": str(
                    bool(manuscript_title_norm)
                    and (
                        manuscript_title_norm in pubmed_title_norm
                        or pubmed_title_norm in manuscript_title_norm
                    )
                ).upper(),
                "crossref_title_contains_manuscript_fragment": str(
                    bool(manuscript_title_norm)
                    and (
                        manuscript_title_norm in crossref_title_norm
                        or crossref_title_norm in manuscript_title_norm
                    )
                ).upper(),
            }
        )

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    crossref_resolved = [row for row in rows if row["crossref_doi_found"] == "TRUE"]
    failures = [
        row
        for row in rows
        if row["pubmed_pmid_found"] != "TRUE"
        or row["pmid_doi_match"] != "TRUE"
        or (
            row["crossref_doi_found"] == "TRUE"
            and row["crossref_doi_match"] != "TRUE"
        )
    ]

    md_lines = [
        "# Reference Verification Audit",
        "",
        "Date: 2026-06-14",
        "",
        "## Scope",
        "",
        f"- Manuscript checked: `{MANUSCRIPT.relative_to(PROJECT_ROOT)}`",
        f"- References parsed: {len(rows)}",
        "- APIs used: NCBI E-utilities PubMed ESummary and CrossRef Works API.",
        "",
        "## Summary",
        "",
        f"- PubMed PMIDs found: {sum(row['pubmed_pmid_found'] == 'TRUE' for row in rows)}/{len(rows)}",
        f"- PubMed DOI matches manuscript DOI: {sum(row['pmid_doi_match'] == 'TRUE' for row in rows)}/{len(rows)}",
        f"- CrossRef DOI records found: {len(crossref_resolved)}/{len(rows)}",
        (
            "- CrossRef DOI matches manuscript DOI among resolved records: "
            f"{sum(row['crossref_doi_match'] == 'TRUE' for row in crossref_resolved)}/{len(crossref_resolved)}"
        ),
        "",
        "## Blocking Findings",
        "",
    ]
    if failures:
        for row in failures:
            md_lines.append(
                f"- Ref {row['ref_no']}: PMID found={row['pubmed_pmid_found']}; "
                f"PubMed DOI match={row['pmid_doi_match']}; "
                f"CrossRef found={row['crossref_doi_found']}; "
                f"CrossRef DOI match={row['crossref_doi_match']}."
            )
    else:
        md_lines.append("- None detected by PMID/DOI machine verification. CrossRef network failures, if present, are treated as unresolved checks rather than identifier mismatches.")

    md_lines.extend(
        [
            "",
            "## Output",
            "",
            f"- Full CSV report: `{OUT_CSV.relative_to(PROJECT_ROOT)}`",
            "",
            "## Limitations",
            "",
            "- This audit verifies identifier existence and DOI agreement. It does not replace target-journal reference formatting.",
            "- Title matching is conservative because manuscript references use shortened titles and `et al.` formatting.",
            "- CrossRef API failures are reported as unresolved network checks unless CrossRef returns a DOI record that conflicts with the manuscript DOI.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(md_lines), encoding="utf-8")

    print(f"references parsed: {len(rows)}")
    print(f"blocking verification failures: {len(failures)}")
    print(f"wrote: {OUT_CSV.relative_to(PROJECT_ROOT)}")
    print(f"wrote: {OUT_MD.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
