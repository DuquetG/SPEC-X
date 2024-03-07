import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sigfig import round

# Traitement des données brutes sous format .mca
# Prend en argument le fichier .mca et retourne une liste de valeurs filtrées
# Détecte automatiquement le nombre de canaux

def traitement(fichier):

    with open(fichier, 'r') as file:
        data = file.read()

    valeurs = re.search(r'<<DATA>>(.*?)<<END>>', data, re.DOTALL)
    liste_valeurs = list(map(float, valeurs.group(1).strip().split("\n")))

    match_canaux = re.search(r'MCA Channels: (\d+)', data)
    if match_canaux:
        canaux = int(match_canaux.group(1))
    else:
        canaux = None

    liste_canaux = range(canaux)

    return [liste_canaux, liste_valeurs]


# Fonction d'initialisation de courbe gaussienne

def gaussienne(x, a, sigma, mu):
    return a * np.exp(-((x - mu) / sigma) ** 2)


# Courbe d'ajustement

def gauss_fit(x_data, y_data, pos):
    popt, pcov = curve_fit(gaussienne, x_data, y_data, p0=[np.max(y_data), np.sqrt(np.std(y_data)), pos])

    fit = gaussienne(x_data, *popt)
    # Calcul du R^2 du fit
    residuals = fit - y_data
    ss_res = np.sum(residuals ** 2)
    ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
    r_squared = 1 - (ss_res / ss_tot)

    plt.plot(x_data, fit, label=f"Courbe d'ajustement du pic d'énergie à ~{pos} canaux", alpha=0.8)

    return popt


# Mention des valeurs approximatives des pics d'énergie pour chaque source en canaux

sources = ["Am", "Cd", "Co"]
pics = {"Am": [178, 225, 750],
        "Cd": [281, 1107],
        "Co": [184, 1534]}


# Valeurs de référence de l'américium en keV

ref = [13.95, 17.74, 59.54]


# Ajustement d'une relation linéaire ax + b

a, b = np.polyfit(pics["Am"], ref, 1)


# Initialisation des listes des valeurs d'énergie et d'incertitude en keV

kev = []
inc = []


# Transformation de canaux à énergie des valeurs pour chaque source

def etalonnage(a, b, pics):
    for i, pic in pics.items():
        energie_list = []
        for ii in pic:
            energie = ii * a + b
            energie_list.append(energie)
            kev.append(energie)
        print(f"Liste des pics d'énergie du {i} en keV: {energie_list}")


# Identification des incertitudes

def incertitudes(a, b, sources, pics):
    for i in sources:
        données = traitement(f"Données\{i}.mca")
        x = données[0]
        y = données[1]

        incertitude_list = []
        for ii in pics[i]:
            popt = gauss_fit(x, y, ii)
            incertitude = a*popt[1]+b
            incertitude_list.append(incertitude)
            inc.append(incertitude)

        print(f"Liste des incertitudes des pics d'énergie du {i} en keV: {incertitude_list}")


# Affichage des spectres bruts sans étalonnage

def graphiques(sources, pics):
    for i in sources:
        données = traitement(f"Données\{i}.mca")
        x = données[0]
        y = données[1]

        plt.plot(x, y, label="Spectre brut")

        for ii in pics[i]:
            gauss_fit(x, y, ii)

        plt.title(f"Graphique d'identification des pics dans le spectre du {i}")
        plt.xlabel('Nombre de canaux')
        plt.ylabel('Nombre de comptes')
        plt.tick_params("both", direction="in")
        plt.legend()
        plt.show()


# Affichage du grapique de l'efficacité en fonction de l'énergie

def efficacité(kev, inc):
    abs = []
    rel = []
    for i in range(len(inc)):
        abs.append(2.35*inc[i])
        rel.append(100*2.35*inc[i]/kev[i])

    for i in [abs, rel]:
        plt.clf()
        plt.xlabel('Énergie [keV]')
        if i == abs:
            plt.ylabel('Résolution absolue [keV]')
            plt.title("Graphique de la résolution absolue en énergie en fonction de l'énergie des photons X et γ pour le CdTe")
        elif i == rel:
            plt.ylabel('Résolution relative [keV]')
            plt.title("Graphique de la résolution relative en énergie en fonction de l'énergie des photons X et γ pour le CdTe")
        plt.tick_params("both", direction="in")
        plt.scatter(kev, i, marker='D')
        plt.show()


# etalonnage(a, b, pics)
# incertitudes(a, b, sources, pics)
# efficacité(kev, inc)

# À exécuter séparément

graphiques(sources, pics)