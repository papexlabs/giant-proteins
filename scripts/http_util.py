"""Small HTTP helper: generic User-Agent, no contact data, simple per-host rate limit.

NCBI E-utilities without an API key allow at most 3 requests/second; we stay at
<= 2.5/s for NCBI hosts. No email/tool parameters are ever sent.
"""
import time
import urllib.parse
import urllib.error
import urllib.request

UA = "PapexResearchBot/1.0 (+https://papex.org)"
_last = {}
_MIN_INTERVAL = {"eutils.ncbi.nlm.nih.gov": 0.4, "www.ncbi.nlm.nih.gov": 0.4,
                 "api.ncbi.nlm.nih.gov": 0.4, "trace.ncbi.nlm.nih.gov": 0.4,
                 "rest.uniprot.org": 0.3, "www.ebi.ac.uk": 0.3}


def get(url, params=None, data=None, retries=6, timeout=300):
    if params:
        url = url + ("&" if "?" in url else "?") + urllib.parse.urlencode(params, doseq=True)
    host = urllib.parse.urlparse(url).netloc
    for attempt in range(retries):
        gap = _MIN_INTERVAL.get(host, 0.1)
        wait = _last.get(host, 0) + gap - time.time()
        if wait > 0:
            time.sleep(wait)
        _last[host] = time.time()
        body = urllib.parse.urlencode(data).encode() if data else None
        req = urllib.request.Request(url, data=body, headers={"User-Agent": UA})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if attempt == retries - 1 or e.code in (400, 404):
                raise
            ra = e.headers.get("Retry-After")
            time.sleep(float(ra) if ra and ra.isdigit() else 10 * (attempt + 1))
        except Exception:  # noqa: BLE001
            if attempt == retries - 1:
                raise
            time.sleep(5 * (attempt + 1))
    return None
