# -*- coding: utf-8 -*-
# =====================================================================
#  Vizualizace mracna bodu jednoho stromu, se zamerenim na libovolnou
#  vybranou oblast (napr. podezrely shluk bodu, ktery zpusobuje
#  nesmyslne velke polomery pri rekonstrukci v AdTree/AdQSM/TreeQSM).
#
#  POUZITI:
#   1) Nastav INPUT_XYZ nize na cestu ke svemu .xyz souboru.
#   2) Nastav Z_LO/Z_HI/X_LO/X_HI/Y_LO/Y_HI na oblast, kterou chces
#      prohlednout (napr. souradnice podezreleho mista z vystupu
#      kalibrace/kontroly poloměru).
#   3) Spust cely skript (F5 ve VS Code, nebo "python
#      visualize_pointcloud_region.py" v terminalu).
#   4) Tri obrazky se ulozi do OUTPUT_DIR.
#
#  Co skript dela, krok za krokem:
#   1) Nacte cely mrak bodu (sloupce x,y,z, bez hlavicky).
#   2) Obrazek 1: bocni pohled na CELY strom (x-z rovina), s cervenym
#      obdelnikem vyznacujicim, kde presne lezi zvolena oblast - dava
#      kontext, kde v celem stromu ta oblast je.
#   3) Obrazek 2: pohled SHORA (x-y rovina) jen na body ve vyskovem
#      pasu Z_LO-Z_HI - ukaze, jestli tvar pripomina kruhovy prurez
#      vetve (ocekavany tvar), nebo je to neuspoladany/rozdvojeny shluk.
#   4) Obrazek 3: bocni pohled (x-z rovina) na uzsi vyrez kolem oblasti
#      (automaticky o kus vetsi nez zvolena oblast), barevne odlisi
#      body podle hloubky (y) - pomaha rozeznat, jestli jde o jednu
#      vetev, nebo o dve ruzne vetve/objekty, ktere se v tomhle pohledu
#      jen prekryvaji.
#
#  Zavislosti: numpy, matplotlib   (install: pip install numpy matplotlib)
# =====================================================================

import numpy as np
import matplotlib.pyplot as plt

# =====================  PARAMETRY  ====================================
# Cesta ke vstupnimu mracnu bodu (format: 3 sloupce x y z, bez hlavicky,
# oddelene mezerou - presne jako "IND07_083_-_Cloud.xyz").
INPUT_XYZ = r"C:\cesta\k\tvemu\souboru\IND07_083_-_Cloud.xyz"

# Kam ulozit vysledne obrazky.
OUTPUT_DIR = "."

# Kazdy N-ty bod se pouzije pro CELKOVY pohled (obrazek 1) - jen kvuli
# rychlosti kresleni u velkych mraku (stovky tisic az miliony bodu).
# Detailni obrazky (2 a 3) pouzivaji VSECHNY body ve sve oblasti.
SUBSAMPLE_EVERY_N = 5

# Souradnice oblasti, kterou chces prohlednout (v METRECH, ve stejnem
# souradnem systemu jako vstupni .xyz - z je typicky uz vyska nad bazi).
Z_LO, Z_HI = 24.0, 27.0
X_LO, X_HI = 4.0, 10.0
Y_LO, Y_HI = 11.5, 13.0

# O kolik metru rozsirit oblast v obrazku 3 (detail zboku) oproti
# zvolenemu ramecku, aby bylo videt i okoli (napr. kam vetev vede dal).
ZOOM_MARGIN_Z = 3.0   # +/- kolem Z_LO..Z_HI
ZOOM_MARGIN_XY = 2.0  # +/- kolem X_LO..X_HI a Y_LO..Y_HI
# =====================================================================


print("Nacitam mrak bodu z '%s' ..." % INPUT_XYZ)
pts = np.loadtxt(INPUT_XYZ)   # shape (N, 3): sloupce x, y, z
print("  nacteno %d bodu" % len(pts))

# =====================================================================
# OBRAZEK 1: cely strom zboku (x-z), s vyznacenou oblasti
# =====================================================================
sub_full = pts[::SUBSAMPLE_EVERY_N]
fig, ax = plt.subplots(figsize=(6, 10))
ax.scatter(sub_full[:, 0], sub_full[:, 2], s=0.5, c="gray", alpha=0.4)
ax.add_patch(plt.Rectangle((X_LO, Z_LO), X_HI - X_LO, Z_HI - Z_LO,
                            fill=False, edgecolor="red", linewidth=2))
ax.set_xlabel("x [m]")
ax.set_ylabel("z / vyska [m]")
ax.set_title("Cely strom (bocni pohled)\ncerveny ramecek = zvolena oblast")
ax.set_aspect("equal")
fig.tight_layout()
fig.savefig(OUTPUT_DIR + "/1_cely_strom_zboku.png", dpi=150)
plt.close(fig)
print("Ulozeno: 1_cely_strom_zboku.png")

# =====================================================================
# OBRAZEK 2: pohled SHORA na zvoleny vyskovy pas, barva = vyska
# =====================================================================
mask_slice = (pts[:, 2] >= Z_LO) & (pts[:, 2] <= Z_HI)
slice_pts = pts[mask_slice]

fig, ax = plt.subplots(figsize=(8, 7))
if len(slice_pts) > 0:
    sc = ax.scatter(slice_pts[:, 0], slice_pts[:, 1], s=2, c=slice_pts[:, 2],
                     cmap="viridis")
    plt.colorbar(sc, ax=ax, label="vyska z [m]")
else:
    print("  POZOR: v zadanem vyskovem pasu nejsou zadne body.")
ax.add_patch(plt.Rectangle((X_LO, Y_LO), X_HI - X_LO, Y_HI - Y_LO,
                            fill=False, edgecolor="red", linewidth=2))
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
ax.set_title("Pohled SHORA, jen body ve vysce %.1f-%.1f m\ncerveny ramecek = zvolena oblast" % (Z_LO, Z_HI))
ax.set_aspect("equal")
fig.tight_layout()
fig.savefig(OUTPUT_DIR + "/2_pohled_shora.png", dpi=150)
plt.close(fig)
print("Ulozeno: 2_pohled_shora.png")

# =====================================================================
# OBRAZEK 3: detail zboku, sirsi vyrez kolem oblasti, barva podle y
# =====================================================================
mask_zoom = ((pts[:, 2] >= Z_LO - ZOOM_MARGIN_Z) & (pts[:, 2] <= Z_HI + ZOOM_MARGIN_Z) &
             (pts[:, 0] >= X_LO - ZOOM_MARGIN_XY) & (pts[:, 0] <= X_HI + ZOOM_MARGIN_XY) &
             (pts[:, 1] >= Y_LO - ZOOM_MARGIN_XY) & (pts[:, 1] <= Y_HI + ZOOM_MARGIN_XY))
zoom_pts = pts[mask_zoom]

fig, ax = plt.subplots(figsize=(9, 7))
if len(zoom_pts) > 0:
    sc = ax.scatter(zoom_pts[:, 0], zoom_pts[:, 2], s=3, c=zoom_pts[:, 1],
                     cmap="plasma")
    plt.colorbar(sc, ax=ax, label="y [m] (hloubka)")
else:
    print("  POZOR: v rozsirene oblasti nejsou zadne body.")
ax.add_patch(plt.Rectangle((X_LO, Z_LO), X_HI - X_LO, Z_HI - Z_LO,
                            fill=False, edgecolor="red", linewidth=2))
ax.set_xlabel("x [m]")
ax.set_ylabel("z / vyska [m]")
ax.set_title("Detail zboku (x-z), barva = y\ncerveny ramecek = zvolena oblast")
fig.tight_layout()
fig.savefig(OUTPUT_DIR + "/3_detail_zboku.png", dpi=150)
plt.close(fig)
print("Ulozeno: 3_detail_zboku.png")

print("\nHotovo. 3 obrazky ulozeny do '%s'." % OUTPUT_DIR)
