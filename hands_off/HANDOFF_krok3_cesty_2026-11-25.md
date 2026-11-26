# PŘEDÁVACÍ ZPRÁVA: tree_reconstruction – krok 3 (sjednocení cest), 2026-11-25

Pro Claude Code ve VS Code na velkém PC. Navazuje na `hands_off/HANDOFF_tree_reconstruction_2026-11-25.md`
(celkový plán, kroky 0–5). Kroky 1 a 2 jsou HOTOVÉ, tato zpráva řeší **krok 3** (a připravuje krok 4).

## Jak se mnou pracovat (důležité)

- Komunikuj **česky**. V Pythonu, Gitu a Terminálu jsem začátečnice: každý příkaz a každou úpravu vysvětli krok za krokem (co dělá a proč).
- Postup vždy: **analýza (jen čtení) → návrh → moje schválení → provedení.**
- **Nic nemazat.** Co je navíc, přesunout do `_KE_SMAZANI/` (je v `.gitignore`).
- Před každou změnou kódu ukaž návrh (diff) a počkej na souhlas. Zachovat styl stávajícího kódu, komentáře srozumitelné pro začátečníka, nastavitelné parametry nahoře ve skriptu.
- `git commit` a `git push` spouštím **sama**, ty mi připravíš a vysvětlíš příkazy.

## Stav k 2026-11-25 17:25

- PC i GitHub na commitu `31d4cda`, `git status` čistý.
  - `8f23ea7` Track all scripts (nový `.gitignore`: sleduje jen kód `*.py *.m *.mac *.md *.bat *.json` + `volume_results.csv` + `csv/volume_summary_settled*.csv`; data, výsledky ANSYS, npz, grafy a `docs/` se ignorují). Původní `.gitignore` je zálohovaný v `_KE_SMAZANI/gitignore_puvodni_2026-11-25.txt`.
  - `31d4cda` Merge changes made on Mac 2026-09-11..15 (ansys_eval): 4 změněné + 6 nových skriptů v `scripts/python/ansys_eval/`, cesty už přepsané na PC.
- Pozor: skripty z Macu mění příponu výstupů `_noacel` → `_no_acel` (výstupy ANSYS vyhodnocení se budou muset přegenerovat).
- `TODO_investigations.md` bod 14: chyby v Poissonových číslech v `scripts/ansys/bk_kp.mac`. **Neopravovat teď**, opravím sama později.
- `runsken.m`, `bk_kp.mac`, `export_geom_ansys.py`: obsahují pracovní nastavení (číslo stromu, `simp_SmallRadii` atd.), která průběžně měním. Nastavení neměnit.

## Úkol: krok 3 – sjednotit cesty

Cíl: stejný kód poběží na PC i na Macu bez přepisování cest.

### 3a. Analýza (jen čtení)
1. `scripts/python/paths.py` obsahuje jen relativní názvy složek (`PLOTS_DIR="plots"`, `CSV_DIR="csv"`, `NPZ_DIR="npz"`, `ANSYS_GEOM_DIR="ansys_geom"`) a funkce `ensure_*()`. Importuje jen `os`. Používá ho 17 skriptů, které předpokládají spuštění se `scripts/python` jako pracovní složkou.
2. Pevné cesty: zhruba 42 řádků ve 29 souborech obsahuje `C:\Users\Spravce\Documents\BARA\01_Skeny_Babice\tree_reconstruction\...` (všechny uvnitř projektu). Najdi je znovu (`*.py`, `*.m`), ať je seznam aktuální.
   Výjimky mimo projekt:
   - `treegraph_export_geom.py`, `treegraph_add_to_comparison.py`: `...\01_Skeny_Babice\01_TreeGraph\treegraph_work\results\...` → přes `local_paths.json`
   - `scratch_order_alignment.py`: `C:\Users\Spravce\AppData\Local\Temp\claude\...` → nechat, jen okomentovat
   - `visualize_pointcloud_region.py`: ukázková cesta `C:\cesta\k\tvemu\souboru\...` → nechat
   - `scripts/matlab/runsken.m`: 2 řádky s cestou → navrhni řešení pro MATLAB zvlášť (např. `fileparts(mfilename('fullpath'))`)
3. Výsledek ukaž v tabulce (soubor, řádek, proměnná, navržená náhrada). **Nic neměnit.**

### 3b. Návrh (schválený v principu, detaily ukázat)
1. **Rozšířit `paths.py`** (nezakládat nový soubor), stávající obsah ponechat:
   ```python
   from pathlib import Path
   import json
   PYTHON_DIR     = Path(__file__).resolve().parent      # scripts/python
   PROJECT_ROOT   = PYTHON_DIR.parent.parent             # tree_reconstruction
   SCRIPTS_DIR    = PROJECT_ROOT / "scripts"
   ANSYS_DIR      = SCRIPTS_DIR / "ansys"
   MATLAB_DIR     = SCRIPTS_DIR / "matlab"
   ANSYS_EVAL_DIR = PYTHON_DIR / "ansys_eval"
   DATA_DIR       = PROJECT_ROOT / "data"
   LOCAL_PATHS_FILE = PROJECT_ROOT / "local_paths.json"
   # + funkce, která local_paths.json načte a při chybějícím klíči vypíše srozumitelnou chybu
   ```
2. **`local_paths.json`** v kořeni projektu (už je v `.gitignore`), klíče např. `scans_dir`, `treegraph_results_dir`. Do Gitu přidat vzor `local_paths.example.json`.
3. **Úprava skriptů:** zachovat názvy proměnných, měnit jen pravou stranu, např.
   `INPUT_FOLDER = str(paths.ANSYS_DIR)`.
   Skripty v `ansys_eval/` nevidí `paths.py` přímo → na začátek přidat `sys.path.insert(0, ...)` s komentářem pro začátečníka.

### 3c. Ověření (před a po úpravě)
1. PŘED úpravou spustit `volume_summary_table.py` (vyrábí `csv/volume_summary_settled*.csv`, které jsou v Gitu) a aspoň jeden skript z `ansys_eval` (např. `load_static_scenarios.py`). Uložit kontrolní součty (SHA-256) výstupů.
2. Upravit cesty.
3. PO úpravě spustit znovu a porovnat. `git status` / `git diff` na CSV musí být prázdný (kromě očekávaných změn).
4. Ověřit, že skripty fungují i při spuštění z jiné pracovní složky (proto `PROJECT_ROOT` z `__file__`).
5. Připravit mi commit C: „Replace hardcoded paths with paths.py + local_paths.json“ a push.

## Krok 4 (připravit, až bude krok 3 hotový)

- Najít na PC, kde jsou skeny po redukci (zatím neznámé; na Macu jsou v iCloudu ve složce Babice `03_skeny/naskenovane_stromy_vse/`):
  - `skeny_2023_olisteny_po redukci`: 13 stromů, vyčištěno 6 (S01, S08, S09, S11, S12, S21) jako `.las/.txt`
  - `skeny_2024_bezlisty_po redukci_Leica`: 13 stromů, jen `no_plate/*.las`, nevyčištěné
- Navrhnout, kam v `tree_reconstruction` patří výsledky sérií B23/B24, **bez přepsání výsledků B21**.
- Výstup pro SAG: **každá série ve vlastním souboru** (`volume_summary_settled_B23.csv`, `_B24.csv` + `_trunk`, `_branch`), stejné sloupce jako dnes:
  `tree,adqsm_variant,treeqsm_tag,adqsm,adtree_raw,adtree_calibrated,treeqsm_optimal,treeqsm_simplified`,
  označení stromů `B23_S01`, `B24_S07`. Nikdy nepřidávat B23/B24 do souborů B21 (SAG převádí `B21_S07` → `BK07` podle čísla stromu, smíchání by vytvořilo duplicitní `BK07`).
- Které skeny zpracovat, rozhodnu sama.
