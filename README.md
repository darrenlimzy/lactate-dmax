# Lactate step test: LT1 and LT2 (Log-Poly-ModDmax)

A browser app for the Lactate Analyzer V2 method. Paste the stage data from a
lactate step test and get LT1 (Log-log), LT2 (Log-Poly-ModDmax), Wmax, the five
training zones, the curve plot and CSV/JSON/PNG downloads.

It runs entirely in the browser through [Shinylive](https://posit-dev.github.io/r-shinylive/)
(R compiled to WebAssembly), so it can be hosted as static files on GitHub Pages
and the data entered never leave the user's computer.

## Folders

- `app/analysis.R`: the calculation, copied from the canonical script
  (`cycling-lab-test` skill, `assets/lactate_analyzer_v2.R`). Change the method
  there first, then copy it here.
- `app/app.R`: the Shiny interface.
- `.github/workflows/deploy.yml`: on every push to `main`, exports `app/` with
  Shinylive and publishes it to GitHub Pages.

The app's output was checked against the canonical script on a reference test
(identical LT1, LT2, Wmax and zones).

## Run locally

```bash
Rscript -e 'shinylive::export("app", "site")'
python3 -m http.server 8765 --directory site
```

Then open http://localhost:8765 in Chrome (it does not work from `file://`).

## Privacy

No athlete data are stored in this repository, and the app sends nothing to a
server: all calculation happens in the viewer's browser.
