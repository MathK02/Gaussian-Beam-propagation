import numpy as np
import matplotlib.pyplot as plt

def gaussian(wo, x, y):
    X, Y = np.meshgrid(x, y)
    X2 = X**2
    Y2 = Y**2
    U = np.exp(-(X2 + Y2) / wo**2)
    return U

def super_gaussian(wo, x, y, n):
    X, Y = np.meshgrid(x, y)
    X2 = X**2
    Y2 = Y**2
    U = np.exp(-((X2 + Y2) / wo**2)**n)
    return U

def aberration_front(x, y, r0, sph, coma, asti0, asti45, defocus):
    N = len(x)
    X, Y = np.meshgrid(x, y)
    M1 = np.double(X**2 + Y**2 <= r0**2)
    r = (X**2 + Y**2)**0.5
    r = r * M1
    u = r / r0

    phi = np.arctan2(X, Y)
    phi[N//2 + 1, N//2 + 1] = 0
    Delta_aberr = sph * 5**0.5 * (6 * (r/r0)**4 - (6 * (r/r0)**2) + 1) + \
                  coma * 8**0.5 * (3 * u**3 - 2 * u) * np.cos(phi) + \
                  asti0 * (6**0.5 * u**2) * np.cos(2 * phi) + \
                  (asti45 * 6**0.5 * u**2) * np.sin(2 * phi) + \
                  defocus * 3**0.5 * (r/r0)**2
    return Delta_aberr

def PorteDisqueR1(N, S, P):
    if N < 2 or N % 2 != 0:
        raise ValueError("N n'est pas un nombre entier PAIR strictement positif")
    dx = 2 * S / N
    x = dx * np.arange(-N/2, N/2)
    y = x
    X, Y = np.meshgrid(x, y)
    M = np.double(X**2 + Y**2 <= P**2)
    return M, x, y

def fonction_transfert_espacelibre(Vx, Vy, d, lambda_):
    Hd = np.exp(1j * np.pi * lambda_ * d * (Vx**2 + Vy**2))
    return Hd

def VisuIdB(*args):
    nargin = len(args)

    if nargin in [2, 3]:  # cas VisuIdB(ImAC,seuil_dB) ou VisuIdB(ImAC,seuil_dB,tronc_dB)
        AC = args[0]
        seuil_dB = args[1]
        if nargin == 3:
            tronc_dB = args[2]
            if not np.isscalar(tronc_dB):
                raise ValueError('La syntaxe "VisuIdB(x,y,ImAC)" est invalide. Il faut préciser "seuil_dB"')
            tronc = True
        else:
            tronc = False
        m, n = AC.shape
        x = [1, n]
        y = [1, m]
        YDir = 'reverse'

    elif nargin in [4, 5]:  # cas VisuIdB(x,y,ImAC,seuil_dB)  ou   VisuIdB(x,y,ImAC,seuil_dB,tronc_dB)
        x = args[0]
        y = args[1]
        YDir = 'normal'
        AC = args[2]
        seuil_dB = args[3]
        if nargin == 5:
            tronc_dB = args[4]
            tronc = True
        else:
            tronc = False

    else:
        raise ValueError('Syntaxe d\'appel incorrect! Il faut entre 2 et 5 arguments en entrée; voir l\'aide-en-ligne.')

    U = np.real(AC * np.conj(AC))
    seuil_eff = 10**(seuil_dB / 10)

    if not tronc:
        # seuil_dB RELATIF:
        if seuil_dB >= 0:
            raise ValueError('L\'argument "seuil_dB" doit être négatif (valeur relative au max à 0dB)')
        plt.figure()
        plt.imshow(10 * np.log10(np.maximum(seuil_eff, U / np.max(U))), cmap='gray', extent=(x[0], x[1], y[0], y[1]))
        plt.colorbar(label='dB_{relatif}')
    else:
        plt.figure()
        plt.imshow(np.minimum(tronc_dB, 10 * np.log10(np.maximum(seuil_eff, U))), cmap='gray', extent=(x[0], x[1], y[0], y[1]))
        plt.colorbar(label='dB_{absolu}')

    plt.gca().set_aspect('auto')
    plt.gca().invert_yaxis() if YDir == 'reverse' else None
    plt.gca().tick_params(direction='out')
    plt.show()
    
    

N = 2**8  # Nombre d'échantillonnage, de 2.^8 pour lambda/4 à 2.^10 pour lambda/1
# /!\ ne pas aller au delà de 2.^9 pour l'affichage 3D
# /!\ ne pas aller au delà de 2.^10 dans tous les cas
S = 0.05  # demi-côté du support (carré) de la pupille (circulaire) [en mètre]
r0 = 0.020  # rayon pupille [en mètre]
f = 0.15  # focale de la lentille [en metre]
lambda_ = 600 * 10**-9  # longueur d'onde de l'onde incidente [en mètre]
l = lambda_  # abrégé

dx = 2 * S / N  # pas d'échantillonnage spatial, abcisse [en mètre]
df = 1 / (dx * N)  # pas d'échantillonnage fréquentiel [en mètre^-1]
x = np.arange(-N / 2, N / 2) * dx  # tableau (1,N) des pas d'échantillonnage spatiaux, centré en 0, de pas dx
y = x  # pas d'échantillonnage spatial, ordonnée [en mètre]
X, Y = np.meshgrid(x, y)  # matrice de taille (N,N).

eps0 = (10**-4) / 5  # pas d'échantillonnage de la défocalisation [en metre]
eps = eps0 * np.arange(-100, 101)  # liste de défocalisation


# onde incidente, onde plane d'amplitude uniforme, sans incidence
wo = 1e-3
w = lambda_ * f / (wo * np.pi)   #waist attendu à l'arrivé pour l'onde gaussienne.
Uinc = gaussian(wo, x, y)  # onde incidente gaussienne
Uinc = super_gaussian(wo,x,y,5)

#Uinc = 1  # onde plane pour étudier les aberrations
Delta_aberr = aberration_front(x, y, r0, 0, 0, lambda_/2, 0, 0)  # coeffs sigma de Zernike (sph, coma, ast0 , ast45, defocus,trefoil) en metre (de l'ordre de lambda/10)

disq, x, y = PorteDisqueR1(N, S, r0)  # génération de la pupille
U0 = disq * Uinc * np.exp(2j * np.pi * Delta_aberr / lambda_)  # onde juste après la lentille avec ajout d'aberrations

M2 = np.fft.ifftshift(U0)  # shiftage dans le plan objet pour annuler la variabilité, éviter le saut de phase.
U_focal = np.fft.fftshift(np.fft.fft2(M2))  # amplitude de l'onde incidente au foyer paraxial, espace réel.

dx1 = df * (lambda_ * f)  # pas d'échantillonnage spatial du plan focal, abcisse
dy1 = dx1  # pas d'échantillonnage spatial du plan focal, ordonnée
x1 = dx1 * np.arange(-N/2, N/2)  # tableau des pas d'échantillonnage spatiaux du plan focal centré en 0, abscisse
y1 = x1  # tableau des pas d'échantillonnage spatiaux du plan focal centré en 0, ordonnée
X1, Y1 = np.meshgrid(x1, y1)
S1 = x1[int(N/2) - 1]  # demi-côté du plan étudié en mètre

U_focal_tild = U0  # Raccourci pour faire moins de calcul. numériquement on a bien U_focal_tild = U0, car double TF.

df1 = 1 / (dx1 * N)  # pas d'échantillonnage fréquentiel dans l'espace de Fourier
vx1 = df1 * np.arange(-N/2, N/2)  # tableau des pas d'échantillonnage fréquentiels dans l'espace de Fourier centrés en 0, abscisse
vy1 = vx1  # tableau des pas d'échantillonnage fréquentiels dans l'espace de Fourier centrés en 0, ordonnée
Vx1, Vy1 = np.meshgrid(vx1, vy1)

# Fonction de déplacement. Permettra de déplacer U_focal_tild d'une quantité "eps0" (troisième argument de Hd) autour du foyer paraxial dans l'espace de Fourier.
Hd = fonction_transfert_espacelibre(Vx1, Vy1, 40e-6, lambda_)  # Faire attention à bien mettre les variables [Vx, Vy] (= meshgrid(vx,vy)) dans les deux premiers arguments de Hd /!\

# Amplitude de l'onde dans l'espace de Fourier, déplacée de "eps0" autour du foyer paraxial.
U_defocal_tild = Hd * U_focal_tild

M2 = np.fft.ifftshift(U_defocal_tild)  # shiftage dans le plan objet pour annuler la variabilité
U_defocal1 = np.fft.fftshift(np.fft.ifft2(M2))  # TF^-1 de U_defocal_tild, retour dans l'espace réel.

Plan_coupe_vertical = np.zeros((N, len(eps)))  # rappel: eps = eps0 * np.arange(-100, 101)
Rayon_Tot = np.zeros((N, N, len(eps)))  # on incrémentera chaque plan calculé pour faire au final un volume d'épaisseur 201

# initialisation des critères d'optimisations selon l'axe Z (epsilon)
rE = abs(S1 / 6)  # definition du cercle d'énergie. Ici, diametre = 1/6 du plan X,Y d'étude
E_encerclee = []
RMS_spot = []

for i in range(len(eps)):  # Remplissage de la matrice Plan_coupe_vertical et de la matrice 3D Rayon_Tot
    Hd = fonction_transfert_espacelibre(Vx1, Vy1, eps[i], lambda_)

    U_defocal_tild = Hd * U_focal_tild

    M2 = np.fft.ifftshift(U_defocal_tild)
    U_defocal = np.fft.fftshift(np.fft.ifft2(M2))

    I_defocal = np.abs(U_defocal)**2

    Rayon_Tot[:, :, i] = U_defocal  # c'est une matrice en 3D regroupant en chaque point l'amplitude du faisceau incident.
    disq1, _, _ = PorteDisqueR1(N, S1, rE)  # Matrice masque disque

    E_encerclee.append(np.sum(I_defocal * disq1) / np.sum(I_defocal))

    RMS_spot.append(np.sqrt(np.sum(((X1**2 + Y1**2) * I_defocal)) / np.sum(I_defocal)))

    for j in range(N):
        Plan_coupe_vertical[j, i] = U_defocal[j, N//2]  # collage pour réaliser une vue en coupe du faisceau le long de l'axe

# -----posbesrms et E encerle----
Best_RMS = min(RMS_spot)
ind_Best_RMS = RMS_spot.index(Best_RMS)
Best_eps_RMS = eps[ind_Best_RMS]
print('Le meilleur foyer RMS est situé en f +', Best_eps_RMS * 1000, 'mm et vaut :', Best_RMS * 10**6, 'micro mètre')

Best_Eenc = max(E_encerclee)
ind_Best_Eenc = E_encerclee.index(Best_Eenc)
Best_eps_Eenc = eps[ind_Best_Eenc]
print('Le meilleur foyer E_{encerclé} pour un rayon de', rE * 10**6, 'micro mètre est situé en f +',
      Best_eps_Eenc * 1000, 'mm et vaut :', Best_Eenc)
# ----------------------------------
# E_encerclé au meilleur foyer RMS en fonction du rayon considéré
I_defocal_meilleur_foyer = np.abs(Rayon_Tot[:, :, ind_Best_RMS])**2

E_encerclee_meilleur_foyer = []
r_E = np.linspace(0, np.abs(S1/4), N//2 + 1)

for i in range(1, N//2 + 2):
    disq2, _, _ = PorteDisqueR1(N, S1, r_E[i - 1])
    E_encerclee_meilleur_foyer.append(np.sum(I_defocal_meilleur_foyer * disq2) / np.sum(I_defocal_meilleur_foyer))

# ----------------------------------------------------------------

eps1 = np.array(eps)*100
RMS_spot1 = np.array(RMS_spot)*10**6

eps2 = np.array(eps)*1000
E_encerclee1 = np.array(E_encerclee)*100

r_E1 = np.array(r_E)*10**6

x1 = np.array(x1)*10**6
y1 = np.array(y1)*10**6

VisuIdB(x1,y1,U_focal,-15)

plt.figure()
plt.plot(eps1, RMS_spot1)
plt.title("rayon RMS selon z")
plt.xlabel("epsilon [mm]")
plt.ylabel("rayon RMS[μm]")
plt.figure()
plt.plot(eps2, E_encerclee1)
plt.title("E_{encerclée} selon z")
plt.xlabel("epsilon [mm]" )
plt.ylabel("E_{encerclée} [%]")
plt.figure()
plt.plot(r_E1, E_encerclee_meilleur_foyer)
plt.title("E_{encerclée} en fonction du rayon")
plt.xlabel("r [\mum]" )
plt.ylabel("E_{encerclée} [%]")


