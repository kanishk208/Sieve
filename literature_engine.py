"""
literature_engine.py — Phase 3: Literature Extraction Engine
Queries PubMed and uses Gemini 1.5 Pro to extract structured evidence.
"""
import os
import json
import requests
from google import genai
from google.genai import types

def fetch_pubmed_abstracts(gene: str, mutation: str, max_results: int = 3) -> str:
    """Queries Entrez E-utilities for the top relevant papers."""
    search_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "pubmed",
        "term": f"{gene} {mutation}",
        "retmode": "json",
        "retmax": max_results
    }
    try:
        res = requests.get(search_url, params=params, timeout=10).json()
        id_list = res.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return ""
        
        # Fetch abstracts for the gathered IDs
        fetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        fetch_params = {
            "db": "pubmed",
            "id": ",".join(id_list),
            "rettype": "abstract",
            "retmode": "text"
        }
        abstracts = requests.get(fetch_url, params=fetch_params, timeout=15).text
        return abstracts
    except Exception as e:
        print(f"[!] PubMed retrieval failed: {e}")
        return ""

def analyze_literature_with_gemini(gene: str, mutation: str, text_data: str) -> dict:
    """Processes abstracts using Gemini Pro to extract structured evidence."""
    if not text_data or not text_data.strip():
        return {
            "concordant_count": 0,
            "discordant_count": 0,
            "experimental_methods": [],
            "summary": "No literature found for this variant."
        }
        
    try:
        # Standard SDK initialization (looks for GEMINI_API_KEY env variable)
        client = genai.Client()
    except Exception as e:
        print(f"[!] Failed to initialize Gemini Client: {e}")
        return {
            "concordant_count": 0,
            "discordant_count": 0,
            "experimental_methods": [],
            "summary": "Literature found, but AI extraction failed (API key missing or error)."
        }
    
    prompt = f"""
    You are an expert biocurator reviewing literature for the variant {gene} {mutation}.
    Analyze the following scientific text and extract data into a valid JSON object with these keys:
    - concordant_count: integer (papers matching a pathogenic/disruptive finding)
    - discordant_count: integer (papers asserting benign/neutral or conflicting findings)
    - experimental_methods: list of strings (e.g., 'In vitro ATPase assay', 'Western blot')
    - summary: a short 2-sentence summary of the consensus.

    Return ONLY a valid JSON object, without markdown formatting blocks.

    Text Data:
    {text_data}
    """
    
    try:
        response = client.models.generate_content(
            model='gemini-1.5-pro',
            contents=prompt,
        )
        
        # Try parsing JSON, stripping any potential markdown formatting
        text = response.text.strip()
        if text.startswith("```json"):
            text = text[7:]
        if text.startswith("```"):
            text = text[3:]
        if text.endswith("```"):
            text = text[:-3]
        
        return json.loads(text.strip())
    except Exception as e:
        print(f"[!] Gemini AI processing failed: {e}")
        return {
            "concordant_count": 0,
            "discordant_count": 0,
            "experimental_methods": [],
            "summary": f"Failed to parse literature: {e}"
        }

def get_literature_evidence(gene: str, mutation: str) -> dict:
    """High-level function to run the full literature pipeline."""
    abstracts = fetch_pubmed_abstracts(gene, mutation)
    evidence = analyze_literature_with_gemini(gene, mutation, abstracts)
    return evidence

if __name__ == "__main__":
    # Test
    print("Testing literature engine for VCP R155H...")
    res = get_literature_evidence("VCP", "R155H")
    print(json.dumps(res, indent=2))
