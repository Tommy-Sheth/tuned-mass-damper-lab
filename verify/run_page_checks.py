"""Open index.html in a headless browser and print the page's own self-checks
(the 'Checks against Python' table) plus the default-storm numbers.
Needs:  pip install playwright && playwright install chromium
Run from the repository root:  python verify/run_page_checks.py
"""
import pathlib, sys
from playwright.sync_api import sync_playwright

page_url = (pathlib.Path(__file__).resolve().parent.parent / "index.html").as_uri()
with sync_playwright() as p:
    b = p.chromium.launch()
    pg = b.new_page(viewport={"width": 1300, "height": 900})
    errors = []
    pg.on("pageerror", lambda e: errors.append(str(e)))
    pg.goto(page_url)
    pg.wait_for_selector("#check-table tbody tr", state="attached", timeout=60000)
    pg.wait_for_timeout(1500)
    rows = pg.evaluate("""() => [...document.querySelectorAll('#check-table tbody tr')]
        .map(r => [...r.children].map(c => c.textContent.trim()))""")
    passed = sum(1 for r in rows if r and r[-1] and ('pass' in r[-1].lower() or '✓' in r[-1] or 'ok' in r[-1].lower()))
    for r in rows:
        print(" | ".join(r))
    print("\nrows:", len(rows), " marked pass:", passed, " page errors:", errors or "none")
    print("summary text:", pg.evaluate("(document.querySelector('#check-summary')||document.querySelector('details summary')||{}).textContent"))
    b.close()
    sys.exit(0 if not errors else 1)
