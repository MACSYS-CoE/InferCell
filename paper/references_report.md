# references.bib verification report (2026-10-23)

Every entry in `references.bib` was checked against the live source below on
2026-10-23. Title, authors, year, venue, volume, issue, pages and DOI are as that
source gives them. Only formatting was changed: capitals protected with braces,
HTML and Unicode turned into LaTeX, and initials punctuated. Every key that
`paper/main.tex` and `paper/sections/` cite at present (20 keys) is in the file; bibtex with icml2026.bst resolves all 20 with no errors or warnings.

## The 4D minimal cell paper: first author

**The first author is Zane R. Thornburg, not Zaida Luthey-Schulten.** Luthey-Schulten
is the last (senior) author. The key is `thornburg2026bringing`.

Evidence:
- CrossRef `https://api.crossref.org/works/10.1016/j.cell.2026.02.009` gives the
  title "Bringing the genetically minimal cell to life on a computer in 4D", Cell
  189(9), 2582–2597.e27, April 2026, with 18 authors. The first is
  `{"given": "Zane R.", "family": "Thornburg", "sequence": "first"}` and the 18th is
  Zaida Luthey-Schulten.
- The same record has `alternative-id: S0092867426001741`. That is the PII in the
  cell.com URL the README gives, `https://www.cell.com/cell/fulltext/S0092-8674(26)00174-1`,
  so the DOI and the URL are the same paper. (cell.com itself returned 403 to curl.)
- The companion bib (`syn3a-ensemble/paper/references.bib`, `thornburg2026`) agrees.
- The README §References ("Luthey-Schulten et al. (2026)") and the outline's
  Literature section are wrong and should say "Thornburg et al. (2026)". I did not
  edit them, as the brief asked.

## Key → verified by

| Key | Verified by |
|---|---|
| karr2012whole | https://api.crossref.org/works/10.1016/j.cell.2012.05.044 |
| karr2012wholecellkb | https://api.crossref.org/works/10.1093/nar/gks1108 (CrossRef gives only published-online 2012-11-21, no print date, so year = 2012; the vol. 41 D1 issue is dated January 2013) |
| macklin2020simultaneous | https://api.crossref.org/works/10.1126/science.aav3751 |
| thornburg2019kinetic | https://api.crossref.org/works/10.3389/fmolb.2019.00130 |
| thornburg2022fundamental | https://api.crossref.org/works/10.1016/j.cell.2021.12.025 |
| thornburg2026bringing | https://api.crossref.org/works/10.1016/j.cell.2026.02.009 |
| thornburg2021minimal | https://api.datacite.org/dois/10.5281/zenodo.5780120 (Zenodo v1.0, 2021; creators "Zanert2" [Thornburg's GitHub handle], "Luthey-Schulten-Lab") |
| hutchison2016design | https://api.crossref.org/works/10.1126/science.aad6253 |
| breuer2019essential | https://api.crossref.org/works/10.7554/eLife.36842 |
| agmon2022vivarium | https://api.crossref.org/works/10.1093/bioinformatics/btac049 |
| kimpson2026whole | `/fred/oz022/tkimpson/syn3a-ensemble/paper/main.tex` (title, author). **Not on arXiv**: `export.arxiv.org/api/query?search_query=ti:"minimal cell in a minute"` and `all:syn3A AND au:Kimpson` both return 0 results. Entered as `@misc`, "In preparation". **FLAG:** the main.tex co-author block is still a TODO, so the author list may change. Update the entry when it is posted. |
| babtie2017how | https://api.crossref.org/works/10.1098/rsif.2017.0237 |
| lubitz2010parameter | https://api.crossref.org/works/10.1021/jp108764b |
| gutenkunst2007universally | https://api.crossref.org/works/10.1371/journal.pcbi.0030189 |
| lindensantangeli2025increasing | https://api.crossref.org/works/10.1038/s41467-025-62415-4 (the outline's "authors to check": Linden-Santangeli, Zhang, Kramer, Rangamani) |
| liu2009modularization | https://api.semanticscholar.org/graph/v1/paper/DOI:10.1214/09-BA404 (author order Liu, Bayarri, Berger). **Note:** CrossRef and OpenAlex list Bayarri, Berger, Liu, which looks like an alphabetical deposit. I followed Semantic Scholar and the usual citation. Projecteuclid.org could not be fetched to settle it. Check against the PDF before submission. |
| plummer2015cuts | https://api.crossref.org/works/10.1007/s11222-014-9503-z |
| jacob2017better | https://export.arxiv.org/api/query?id_list=1708.08719 (no journal_ref; arXiv preprint) |
| carmona2020semimodular | https://proceedings.mlr.press/v108/carmona20a.html (PMLR 108:4226–4235). PMLR gives "Christian Carmona"; arXiv 2003.06804 gives "Chris U. Carmona". I used PMLR. |
| yu2023variational | https://api.crossref.org/works/10.1214/23-sts886 |
| chakraborty2023modularized | https://api.crossref.org/works/10.1007/s11222-023-10207-5 (article no. 33) |
| liu2025general | https://api.crossref.org/works/10.1093/jrsssb/qkaf012 |
| frazier2025cutting | https://api.crossref.org/works/10.1214/24-ba1448 |
| goudie2019joining | https://api.crossref.org/works/10.1214/18-BA1104. **Correction:** it is in *Bayesian Analysis* 14(1), not *JASA* as the outline says, and it has a fifth author, Wernisch. |
| manderson2023combining | https://api.crossref.org/works/10.1214/22-BA1327 |
| minka2001expectation | https://export.arxiv.org/abs/1301.2294 (comment: "Appears in Proceedings of the Seventeenth Conference on Uncertainty in Artificial Intelligence (UAI2001)"). No page numbers are given there, so none are recorded. |
| barthelme2014expectation | https://api.crossref.org/works/10.1080/01621459.2013.864178 |
| vehtari2020expectation | https://jmlr.org/papers/v21/18-817.html |
| andrieu2010particle | https://api.crossref.org/works/10.1111/j.1467-9868.2009.00736.x |
| lindsten2014particle | https://jmlr.org/papers/v15/lindsten14a.html (and arXiv 1401.0604 journal_ref) |
| andrieu2009pseudomarginal | https://api.crossref.org/works/10.1214/07-AOS574 |
| doucet2015efficient | https://api.crossref.org/works/10.1093/biomet/asu075 (CrossRef gives initials only; recorded as given) |
| golightly2011bayesian | https://api.crossref.org/works/10.1098/rsfs.2011.0047 (the outline's "Golightly and Wilkinson. Particle MCMC for stochastic kinetic models") |
| hobolth2009simulation | https://api.crossref.org/works/10.1214/09-AOAS247 |
| gillespie1977exact | https://api.crossref.org/works/10.1021/j100540a008 |
| wilkinson2018stochastic | https://api.crossref.org/works/10.1201/9781351000918 (3rd ed., Chapman and Hall/CRC, 2018) |
| vankampen2007stochastic | https://api.crossref.org/works/10.1016/B978-0-444-52965-7.X5000-4 (title, Elsevier, 2007, ISBN 9780444529657; CrossRef has **no author field**). The author and "Third Edition" come from https://openlibrary.org/isbn/9780444529657.json and https://openlibrary.org/authors/OL2708276A.json (N.G. Van Kampen). |
| schnoerr2017approximation | https://api.crossref.org/works/10.1088/1751-8121/aa54d9 |
| wood2010statistical | https://api.crossref.org/works/10.1038/nature09319 |
| price2018bayesian | https://api.crossref.org/works/10.1080/10618600.2017.1302882 |
| cranmer2020frontier | https://api.crossref.org/works/10.1073/pnas.1912789117 |
| geffner2023compositional | https://proceedings.mlr.press/v202/geffner23a.html (PMLR 202:11098–11116) |
| kennedy2001bayesian | https://api.crossref.org/works/10.1111/1467-9868.00294 |
| kyzyurova2018coupling | https://api.crossref.org/works/10.1137/17m1157702 |
| ming2021linked | https://api.crossref.org/works/10.1137/20m1323771 |
| hoffman2014nouturn | https://jmlr.org/papers/v15/hoffman14a.html |
| vehtari2021rank | https://api.crossref.org/works/10.1214/20-BA1221 (CrossRef's garbled "Rˆ" written as `$\widehat{R}$`) |
| talts2018validating | https://export.arxiv.org/api/query?id_list=1804.06788 (no journal_ref; arXiv preprint) |
| sailynoja2022graphical | https://api.crossref.org/works/10.1007/s11222-022-10090-6 (article no. 32) |
| saltelli2010variance | https://api.crossref.org/works/10.1016/j.cpc.2009.09.018 |
| sobol2001global | https://api.crossref.org/works/10.1016/S0378-4754(00)00270-6 |
| bezanson2017julia | https://api.crossref.org/works/10.1137/141000671 |
| rackauckas2017differentialequations | https://api.crossref.org/works/10.5334/jors.151 |
| revels2016forwardmode | https://export.arxiv.org/api/query?id_list=1607.07892 (no journal_ref; arXiv preprint) |
| ge2018turing | https://proceedings.mlr.press/v84/ge18b.html (AISTATS 2018, PMLR 84:1682–1690) |

## Unverified or dropped (not in the .bib)

These are named in the candidate lists, mostly in §8 of
`dev/notes/modular-bayesian-inference-heterogeneous-modules.md`. They were not
verified in this pass, so they are not in the .bib. Add any of them after checking
it against a source.

- Carmona & Nicholls (2022), variational meta-posteriors for SMI (§8).
- Pompe & Jacob (2021), cut-posterior asymptotics and posterior bootstrap (§8).
- Poole & Raftery (2000), Bayesian melding (§8).
- Manderson & Goudie (2022), the numerically stable two-stage melding algorithm (§8).
- Liu & Goudie, arXiv 2605.22301, divide-and-conquer SMC for chained melding (§8).
- Heess, Tarlow & Winn (2013); Jitkrittum et al. (2015): learned EP message operators (§8).
- Bagaev & de Vries, reactive message passing / RxInfer.jl (§8).
- Linhart et al. (TMLR); Arruda et al. (arXiv 2505.14429); Gloeckler et al. (ICLR 2025); Gloeckler et al. (Simformer); Radev et al. (2023); Ward et al. (2022): compositional and amortised SBI (§8).
- Arya, Schauer, Schäfer & Rackauckas (2022), StochasticAD; "A differentiable Gillespie algorithm" (eLife 2025 / arXiv 2407.04865) (§8).
- Ming & Williamson (arXiv 2306.01212), linked deep GPs (§8).
- Bissiri, Holmes & Walker (2016); Miller & Dunson (2019) (§8 foundations).
- Lebl et al.: not in any of the candidate lists, and no reference was found to check. Dropped.

Notes:
- `karr2012wholecellkb` (WholeCellKB) is not in the candidate lists either. It was
  added because the brief asked for it, and is verified. Delete it if the paper
  does not cite it.
- Companion-bib entries that are not in InferCell's candidate lists were not
  carried over: sanghvi, elowitz, swain, roberts, petzold, virtanen, davis,
  salmon, holm, clopper, efron, loman, ma2021mtk, hill, montoison, utkarsh,
  kimpson2026lowprec and the `infercell` software entry.
- The companion bib's keys (`thornburg2022`, `karr2012`, …) were renamed to this
  repo's `firstauthorYEARword` convention. The companion's `minimalcell2021` is
  `thornburg2021minimal` here.
