# PŘEDÁVACÍ ZPRÁVA: sjednocení tree_reconstruction (velké PC), 2026-11-25

Pro projekt **Skeny**. Po dokončení se vracím do projektu **Úklid PC → SAG** (viz poslední oddíl).

## Kdo jsem a jak se mnou pracovat

- Výzkumnice a pedagožka na LDF MENDELU, měření stromů, laserové skeny, ANSYS. Komunikuj česky.
- V Terminálu, Gitu a Pythonu jsem začátečnice. Každý příkaz vysvětli krok za krokem: co dělá a proč.
- Postup vždy: **analýza (jen čtení) → návrh → moje schválení → provedení.** Nic nemazat; co je navíc, přesunout do složky `_KE_SMAZANI`.
- Před každou změnou kódu mi ukaž návrh a počkej na souhlas. Zachovat styl stávajícího kódu, komentáře srozumitelné pro začátečníka, nastavitelné parametry nahoře ve skriptu.

## Proč to řešíme

V projektu SAG (`freq_from_dendrometry`) chci přidat objemy z rekonstrukcí **skenů po redukci** (série 2023 olistěné a 2024 bez listí). Rekonstrukce dělám v `tree_reconstruction`. Ten je ale teď rozdělený na dvě verze:

- **Velké PC (Windows):** hlavní verze, zdrojová data, push na GitHub (repo `tree_reconstruction`).
- **Mac:** kopie z 2026-09-11, na které jsem v září pracovala. Změny jsou zabalené v zipu `PRENOS_tree_reconstruction_2026-10-27.zip`. Na velkém PC už zip je, uložila jsem ho pod názvem `PRENOS_tree_reconstruction_2026_10_27`.

Nové rekonstrukce chci dělat až ve **sloučené verzi**, aby byly objemy po redukci srovnatelné s objemy ze skenů 2021 (stejné skripty, stejné parametry adQSM/TreeQSM/AdTree).

## Co je v zipu z Macu (podle `ZMENY_NA_MACU.txt` uvnitř zipu)

- 737 souborů, kořen `tree_reconstruction/`.
- **Commity na Macu** (24 commitů, 2026-09-14 až 2026-09-15, od `d8c6b2c` po `89f8fba`). Hlavně `scripts/python/ansys_eval/`:
  - nové: `check_png_provenance.py`, `check_scenario_dicts_in_sync.py`, `plot_frequency_spectra.py`
  - změněné: `load_static_scenarios.py`, `method_deviation_overview.py`, `plot_modal_harmonic.py`, `plot_stiffness_comparison.py`
- **Neuložené změny** (necommitnuté): `plot_modal_harmonic.py`, `plot_stiffness_comparison.py` a nové `make_boxplot.py`, `make_size_category_plot.py`.
- 268 změněných výstupů (CSV, grafy). Ty se dají znovu vygenerovat, nepřebírat je naslepo.
- **Pozor:** ve skriptech z Macu jsou cesty `/Users/bara/...`. Převzít jen věcné změny, ne cesty.
- Commity z Macu **nejsou na GitHubu**. Push se vždy dělal z velkého PC, takže historie na Macu začíná vlastním „Initial commit“ a nedá se jednoduše `git pull`/`merge`. Změny se přebírají porovnáním souborů.

## Postup na velkém PC

### 0. Příprava
1. Na velkém PC otevřít aplikaci Claude a projekt **Skeny**, tuto zprávu přiložit.
2. Připojit složku `tree_reconstruction` (tlačítko „Add folder“) a složku, kde je rozbalený `PRENOS_tree_reconstruction_2026_10_27`.
3. Git příkazy, které mění repozitář (`pull`, `commit`, `push`), spouštím **sama ve VS Code nebo v Terminálu**, Claude mi je připraví a vysvětlí. (Claude ve svém prostředí nemá povoleno mazat soubory, a git při těchto příkazech maže pomocné soubory `.lock`.)

### 1. Analýza (jen čtení)
1. `git status` a `git log --oneline -10` na velkém PC: jsou tu necommitnuté změny? Shoduje se poslední commit s GitHubem (`git fetch`, pak `git status`)?
2. Porovnat **velké PC** a **zip z Macu** soubor po souboru (hash SHA-256), jen skripty (`*.py`, `*.mac`, `*.json`, `*.md`), bez výstupů:
   - soubory jen na Macu,
   - soubory změněné jen na Macu (na PC od 2026-09-11 beze změny) → bezpečně převzít,
   - soubory změněné **na obou** stranách → konflikt, ukázat mi rozdíly a rozhodnout ručně.
3. Najít pevně zapsané cesty (`C:\...`, `/Users/...`) ve všech skriptech. Na Macu jich bylo cca 45 řádků ve 32 skriptech (29 souborů ve `scripts/python`). Ve `scripts/python/` už existuje `paths.py`: zjistit, co obsahuje a kdo ho používá.
4. Výsledek analýzy shrnout v tabulce a navrhnout postup. **Nic neměnit.**

### 2. Sloučení (po mém schválení)
1. `git pull` z GitHubu (aby PC mělo poslední stav).
2. Převzít změny z Macu (VS Code → Compare, nebo zkopírovat soubory, kde nejde o konflikt). Cesty `/Users/bara/...` nepřebírat.
3. Commit (jasná zpráva, např. „Merge changes made on Mac 2026-09-14..15 (ansys_eval)“) a push.

### 3. Sjednocení cest (po mém schválení, samostatný commit)
1. `scripts/python/project_paths.py` (nebo rozšířit stávající `paths.py`, podle analýzy): kořen projektu si dopočítá sám přes `Path(__file__).resolve()`.
2. Cesta k datům (skenům) v místním souboru `local_paths.json`, zapsaném v `.gitignore` → stejný kód poběží na PC i na Macu.
3. Upravit pevné cesty ve skriptech na `project_paths`. Po úpravě spustit pár klíčových skriptů a ověřit, že výstupy (např. `csv/volume_summary_settled*.csv`) vycházejí stejně jako před úpravou.
4. Commit a push.

### 4. Příprava na rekonstrukce skenů po redukci
1. Navrhnout, kam v `tree_reconstruction` patří výsledky nových sérií, **bez přepsání výsledků B21**.
2. Výstup objemů pro SAG má mít **stejné sloupce** jako dnešní soubor:
   `tree,adqsm_variant,treeqsm_tag,adqsm,adtree_raw,adtree_calibrated,treeqsm_optimal,treeqsm_simplified`
   a označení stromů ve tvaru `B23_S01`, `B24_S07` (série + číslo stromu, jako `B21_S07`).
3. **Každá série do samostatného souboru**, např. `volume_summary_settled_B24.csv` (+ `_trunk`, `_branch`). Nikdy nepřidávat řádky B23/B24 do souborů B21. Skript v SAG (`build_crown_area_volume_mass_table.py`) převádí `B21_S07` → `BK07` podle čísla stromu, takže `B24_S07` ve stejném souboru by vytvořil duplicitní `BK07`.
4. Které skeny zpracovat, rozhodnu sama. K dispozici (vše jen v iCloudu na Macu, ve složce Babice `03_skeny/naskenovane_stromy_vse/`; na PC zjistit, kde jsou):
   - `skeny_2023_olisteny_po redukci`: 13 stromů, vyčištěno 6 (S01, S08, S09, S11, S12, S21) jako `.las/.txt`
   - `skeny_2024_bezlisty_po redukci_Leica`: 13 stromů, jen `no_plate/*.las`, nevyčištěné

### 5. Úklid po sloučení
- Na velkém PC: rozbalený `PRENOS_...` přesunout do `_KE_SMAZANI` (mažu sama).
- Na Macu (v projektu Úklid PC): `01_VYZKUM/ERC_CZ_LL1909/04_skripty_python/tree_reconstruction_MAC_kopie` a zip v `_INBOX` přesunout do `_KE_SMAZANI`.

## Kontrola, že je hotovo

- [ ] Velké PC a GitHub jsou na stejném commitu, `git status` čistý.
- [ ] Změny z Macu převzaté (seznam souborů), konflikty vyřešené.
- [ ] Žádné pevné cesty ve skriptech, `local_paths.json` v `.gitignore`.
- [ ] Výstupy B21 po úpravě cest stejné jako před ní.
- [ ] Je jasné, kam půjdou výsledky B23/B24 a v jakém formátu je dostane SAG.

## Návrat do SAG (projekt Úklid PC)

Až bude hotovo, vracím se do projektu Úklid PC a říkám, v jakém stavu je `tree_reconstruction`. Tam čeká (připraveno, **zatím neprovedeno**):

1. Vytvořit složky v `01_VYZKUM/SAG/04_skripty_python/freq_from_dendrometry/`:
   - `tools/plochy_rekonstrukce/outputs_B23_olisteny_po_redukci/`
   - `tools/plochy_rekonstrukce/outputs_B24_bezlisty_po_redukci/`
   - `data_raw/babice_dynatree/rekonstrukce_po_redukci/`
2. Do `shared/paths.py` přidat složky sérií skenů (`SCANS_B21_DIR`, `SCANS_B23_DIR`, `SCANS_B24_DIR`, odvozené od `EXTERNAL_SCANS_DIR`) a výstupní složky.
3. V `tools/plochy_rekonstrukce/tree_area.py` nahradit pevnou cestu `FILE_PATH` parametry „série + soubor“ a přesměrovat výstupy podle série.
4. Spočítat dendrometrické parametry (plochy, rozměry) ze skenů po redukci.
5. Doplnit objemy z rekonstrukcí po redukci a do analýzy (A03) přidat porovnání buků před a po redukci.
6. Commit dvou opravených README v SAG (`README.md`, `tools/plochy_rekonstrukce/README_SETUP.md`) udělat sama.

Stav úklidu Documents k dnešku: viz `HANDOFF_2026-11-25.md` v projektu Úklid PC. Dnes navíc přesunuto lékařské PDF do `07_OSOBNI` (záznam `_analyza_struktury/zaznam_presunu/2026-11-25_D11_pdf_osobni.csv`).
