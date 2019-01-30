import numpy as np
import matplotlib.pyplot as plt


# LINK: http://blog.bharatbhole.com/creating-boxplots-with-matplotlib/

PAR = { 's': 10 ** 4,  # externer Nährstoff
        'dm': 0.1,  # mRNA-Abbaurate
        'ns': 0.5,  # Nährstoffeffizienz
        'nr': 7459,  # Ribosomenlänge
        'nx': 300,  # Länge nicht-ribosomaler Proteine
        'ymax': 1260,  # max. übersetzen Dehnungsrate
        'Ky': 7,  # Übersetzung Verlängerungsschwelle
        'vt': 726,  # max. Nährstoffimportrate
        'Kt': 1000,  # Nährstoffimportschwelle
        'vm': 5800,  # max. enzymatische Rate
        'Km': 1000,  # enzymatic threshold
        'wr': 930,  # max. Ribosomen-Transkriptionsrate
        'we': 4.14,  # max. Enzymtranskriptionsrate   we = wt = wm
        'wq': 948.93,  # max. q-Transkriptionsrate
        'thetar': 426.87,  # Ribosomentranskriptionsschwelle
        'thetanr': 4.38,  # Nicht-ribosomale Transkriptionsschwelle
        'Kq': 152.219,  # q-Autoinhibitionsschwelle
        'hq': 4,  # q-Autoinhibition Hill-Koeffizient
        'kb': 1,  # mRNA-Ribosomen-Bindungsrate
        'ku': 1,  # mRNA-Ribosomen-Nichtbindungsrate
        'M': 10 ** 8,  # total cell mass
        'kcm': 0.00599,  # Chloramphenicol-Bindungsrate
        'thetax': [4.38, 4.38, 426.87],# 4.38],  #transkriptionswelle #transcriptional energy-thresholds
        'wx': [4.15, 4.15, 930]#, 948.93]  # wt,wm,wr,wq transkriptionsraten #maximal transcription rates
}

## Create data
np.random.seed(10)
ky = [PAR["Ky"]]
kq = [PAR["Kq"]]
thetar = [PAR["thetar"]]
thetanr = [PAR["thetanr"]]
wr = [PAR["wr"]]
wq = [PAR["wq"]]
we = [PAR["we"]]
kcm = [PAR["kcm"]]
y = [310.96192, 432.97502, 0.0414, 471.364, 471.364, 471.364, 49.333, 49.333, 851.4558, 881.1960, 0, 357.898, 348.675, 0] # ribosome-bound mRNA



## combine these different collections into a list
data_to_plot = [ky, kq, PAR["thetax"], wr, wq, we, kcm, y]

# Create a figure instance
fig = plt.figure(1, figsize=(9, 6))

# Create an axes instance
ax = fig.add_subplot(111)

# Create the boxplot
bp = ax.boxplot(data_to_plot)

# Save the figure
#fig.savefig('fig1.png', bbox_inches='tight')

plt.show()