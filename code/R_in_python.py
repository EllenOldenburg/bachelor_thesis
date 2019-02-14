"""

imports

"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate


def Integration(x0):
    Integ = scipy.integrate.ode(changeValues).set_integrator('lsoda').set_initial_value(x0, 0)
    Integ.set_f_params(PAR)
    return Integ


def timeCourse(t, x0):
    #x0 = INPUT_VALUES
    integrator = Integration(x0)
    array = [x0]
    cnt = 1

    while cnt < len(t):
        array.append(integrator.integrate(t[cnt]))
        cnt += 1

    return array


"""

Constants

"""

PAR = { #'s': 10 ** 4,  # externer Nährstoff
        'kin': 0, # substrate supply rate
        'dn': 0,  # 1/min, culture dilution rate
        'dm': 0.1,  # mRNA-Abbaurate
        'ns': 0.5,  # Nährstoffeffizienz
        'nr': 7459,  # Ribosomenlänge
        'nt': 300,  # aa/molecs, non-ribo proteins
        'nm': 300,  # "
        'nq': 300,  # "
        'nx': 300,  # Länge nicht-ribosomaler Proteine
        'gmax': 1260,  # max. übersetzen Dehnungsrate
        'Kg': 7,  # Übersetzung Verlängerungsschwelle
        'vt': 726,  # max. Nährstoffimportrate
        'Kt': 1000,  # Nährstoffimportschwelle
        'vm': 5800,  # max. enzymatische Rate
        'Km': 1000,  # enzymatic threshold
        'wr': 930,  # max. Ribosomen-Transkriptionsrate
        'wt': 4.14,
        'wm': 4.14,
        'wq': 948.93,  # max. q-Transkriptionsrate
        'Tr': 426.87,  # molecs/cell, ribo transcr threshold
        'Tt': 4.38,  # molecs/cell, non-ribo transcr threshold
        'Tm': 4.38,  # "
        'Tq':4.38,  # "
        'Kr': float('inf'),  # molecs/cell, non-q autoinhib threshold
        'Kt': float('inf'),  # "
        'Km': float('inf'),  # "
        'Kq': 152219,  # q-Autoinhibitionsschwelle
        'hr': 1,  # , non-q autoinhib Hill coeff
        'ht': 1,  # "
        'hm': 1,  # "
        'hq': 4,  # , q autoinhib Hill coeff
        'kb': 1,  # cell/(min molecs),  mRNA-ribo binding rate
        'ku': 1,  # 1/min, mRNA-ribo unbinding rate
        'M': 10**8,  # aa, total cell mass
        'kcm': .00599,  # min/uM, chloramphenicol binding rate
        'hq': 4,  # q-Autoinhibition Hill-Koeffizient
        'kb': 1,  # mRNA-Ribosomen-Bindungsrate
        'ku': 1,  # mRNA-Ribosomen-Nichtbindungsrate
        'M': 10 ** 8,  # total cell mass
        'kcm': 0.00599,  # Chloramphenicol-Bindungsrate
        'thetax': [4.38, 4.38, 426.87, 4.38],  # transkriptionswelle #transcriptional energy-thresholds
        'wx': [4.15, 4.15, 930, 948.93],  # wt,wm,wr,wq transkriptionsraten #maximal transcription rates
        'l': 1
}




INPUT_VALUES = [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, .001, 10**4,]

"""

Functions

"""


def vimp(et, s, par):
    return et * ((par['vt'] * s) / (par['Kt'] + s))


def vcat(em, si, par):
    return em * ((par['vm'] * si) / (par['Km'] + si))


# leitet die nettorate der Translation eines Proteins x ab
# 4x
def vx(a, cx, par):
    return cx * (gamma(a, par) / par['nx'])


def vr(a, par, cr):
    return cr * (gamma(a, par) / par['nx'])


def vq(a, cq, par):
    return cq * (gamma(a, par) / par['nx'])


# γ ist die Geschwindigkeit der Translationsdehnung mit maximaler Rate γ max = k 2 und Schwelle K γ = k 2 / K p für halbmaximale Dehnung.
def gamma(a, par):
    return (par['gmax'] * a) / (par['Kg'] + a)


# Wenn wir davon ausgehen, dass der Energieverbrauch in jedem Dehnungsschritt konstant ist, folgt daraus die effektive Transkriptionsrate
# 4x
def omegax(a, wx, thetaX):
    omegaResult = []
    for i in range(0, 3):
        omegaResult.append(wx[i] * (a / (thetaX[i] + a)))

    return omegaResult


def I(q, par):
    return 1 / (1 + ((q / par['Kq']) ** par['hq']))


def omegaq(a, par, q):
    return par['wq'] * (a / (par['Tq'] + a)) * I(q, par)


# Die Wachstumsrate λ ist entscheidend, um die zellulären Prozesse mit Wachstum zu verbinden, da alle intrazellulären Spezies durch Umverteilung des Zellinhalts zwischen Mutter- und Tochterzellen verdünnt werden
# nur im steady state
def lamda(a, par, cValues):
    return (gamma(a, par) / par['M']) * sum(cValues)


# die Gesamtmasse der Zelle als Gesamtproteinmasse (einschließlich gebundener Ribosomen)
# def M(par):
# return(sum(par['nx']*x)+pa['nr']*sum(cx))

# def M(cValue, par):  # die Gesamtmasse der Zelle als Gesamtproteinmasse (einschließlich gebundener Ribosomen)
# return (sum(par['nx'] * x) + par['nr'] * sum(cValue))

def M(cValues, r, et, q, em, par):
    return par['nr'] * sum(cValues) + par['nr'] * r + par['nx'] * et + par['nx'] * q + par['nx'] * em


def nr_r(cValues, et, q, em, par):
    return - (par['nr'] * sum(cValues) + par['nx'] * et + par['nx'] * q + par['nx'] * em + par['M'])


"""

Differentialgleichungen


"""


# Differentialgleichung für den inneren Nährstoff
# 1x
def dsi_dt(si, s, par, et, em, lamdaResult):
    return vimp(et, s, par) - vcat(em, si, par) - lamdaResult * si


# Gleichung für die zellulare Energie
# 1x
def da_dt(a, em, cx, si, par, lamdaResult):
    sum_all_protein_in_cell = []
    for i in range(0, 4):
        sum_all_protein_in_cell.append(par['nx'] * vx(a, cx[i], par))
    return par['ns'] * vcat(em, si, par) - sum(sum_all_protein_in_cell) - (lamdaResult * a)


# Gleichung für freie Ribosomen  für vr und cr
# 1x
def dr_dt(a, cx, mx, r, par, lamdaResult, cr):
    sum_proteins_ribosomes = []
    for i in range(0, 4):
        sum_proteins_ribosomes.append(vx(a, cx[i], par) - par['kb'] * r * mx + par['ku'] * cx[i])
    return vr(a, par, cr) - lamdaResult * r + sum(sum_proteins_ribosomes)


# intrazelluläre Ribosomen Sei r die Anzahl der freien Ribosomen. Bezeichne k b und k u die Bindungs- und Unbindungsraten eines Ribosoms an mRNA (für alle mRNAs als identisch angenommen) und lasse die mRNA für ein Protein x dann m x sein
# 4x
def dmx_dt(a, cx, mx, r, par, lamdaResult, omegaResult):
    return omegaResult - (lamdaResult + par['dm']) * mx + vx(a, cx, par) - par['kb'] * r * mx + par['ku'] * cx


def dmq_dt(a, mq, r, par, lamdaResult, q, cq):
    return omegaq(a, par, q) - (lamdaResult + par['dm']) * mq + vq(a, cq, par) - par['kb'] * r * mq + par['ku'] * cq


# Mit cx wird der Komplex zwischen einem Ribosom bezeichnet und die mRNA für Protein x
# 4x
def dcxt_dt(a, cx, mx, r, par, lamdaResult):
    return (-lamdaResult * cx) + par['kb'] * r * mx - par['ku'] * cx - vx(a, cx, par)


# Transporterenzyme für vt und ct
# 1x
def det_dt(a, ct, et, par, lamdaResult):
    return vx(a, ct, par) - lamdaResult * et


# Enzym,das si in a umwandelt (metabolische Enzyme)
# 1x
def dem_dt(a, cm, em, par, cValues):
    return vx(a, cm, par) - lamda(a, par, cValues) * em


# house-keeping Proteine
# 1x
def dq_dt(a, cq, q, par, cValues):
    return vx(a, cq, par) - lamda(a, cValues, par) * q


def dn_dt(lamdaResult, N, par):
    return lamdaResult * N - par['dn'] * N


def ds_dt(par, et, s, N):
    return par['kin'] - vimp(et, s, par) * N - par['dn'] * s


"""

Change INPUT_VALUES

"""

"""

Change INPUT_VALUES

"""


def changeValues(time, i, par):
    si = i[0]
    a = i[1]
    et = i[3]
    q = i[5]
    em = i[4]
    r = par['M'] / par['nr']
    # r = i[2]
    mt = i[6]
    mm = i[7]
    mr = i[8]
    mq = i[9]
    ct = i[10]
    cm = i[11]
    cr = i[12]
    cq = i[13]
    s = i[14]
    N = i[15]
    cx = [ct, cm, cr, cq]
    # r = nr_r(cx, et, q, em, par)

    omegaResult = omegax(a, par["wx"], par["thetax"])
    # lamdaResult = lamda(a, par, cx)
    lamdaResult = par['l']

    detResult = det_dt(a, ct, et, par, lamdaResult)
    demResult = det_dt(a, cm, em, par, lamdaResult)
    dqResult = det_dt(a, cq, q, par, lamdaResult)
    dsiResult = dsi_dt(si, s, par, et, em, lamdaResult)
    daResult = da_dt(a, em, cx, si, par, lamdaResult)
    drResult = dr_dt(a, cx, mt, r, par, lamdaResult, cr)
    dmtResult = dmx_dt(a, ct, mt, r, par, lamdaResult, omegaResult[0])
    dmmResult = dmx_dt(a, cm, mm, r, par, lamdaResult, omegaResult[1])
    dmrResult = dmx_dt(a, cr, mr, r, par, lamdaResult, omegaResult[2])
    dmqResult = dmq_dt(a, mq, r, par, lamdaResult, q, cq)
    dctResult = dcxt_dt(a, ct, mt, r, par, lamdaResult)
    dcmResult = dcxt_dt(a, cm, mm, r, par, lamdaResult)
    dcrResult = dcxt_dt(a, cr, mr, r, par, lamdaResult)
    dcqResult = dcxt_dt(a, cq, mq, r, par, lamdaResult)
    dsResult = ds_dt(par, et, s, N)
    dnResult = dn_dt(lamdaResult, N, par)

    return [
        dsResult,
        dnResult,
        dsiResult,
        daResult,
        drResult,
        detResult,
        demResult,
        dqResult,
        dmtResult,
        dmmResult,
        dmrResult,
        dmqResult,
        dcmResult,
        dcrResult,
        dctResult,
        dcqResult
    ]

'''

ode_weisse < - function(time,
                        init=c(N=0.001, s=1e4, si=1, a=0,
                               et=0, em=0, eq=0, er=0,
                               ct=0, cm=0, cq=0, cr=0,
                               mt=0, mm=0, mq=0, mr=0),
                        parms=c(kin=0,  # substrate supply rate
                                dn=0,  # 1/min, culture dilution rate
                                dm=0.1,  # 1/min, mRNA degradation rate
                                ns=0.5,  # , nutrient efficiency
                                nr=7459,  # aa/molecs, ribo protein length
                                nt=300,  # aa/molecs, non-ribo proteins
                                nm=300,  # "
                                nq=300,  # "
                                gmax=1260,  # aa/(min molecs), max transl.elong
                                Kg=7,  # molecs/cell, transl.elong threshold
                                vt=726,  # 1/min, max nutrient import rate
                                Kt=1000,  # molecs, nutrient import threshold
                                vm=5800,  # 1/min, max enzymatic rate
                                Km=1000,  # molecs/cell, enzymatic threshold
                                wr=930,  # molecs/(min cell), max ribo transcr
                                wt=4.14,  # molecs/(min cell), max non-ribo transc
                                wm=4.14,  # "
                                wq=948.93,  # ", but proteins q
                                Tr=426.87,  # molecs/cell, ribo transcr threshold
                                Tt=4.38,  # molecs/cell, non-ribo transcr threshold
                                Tm=4.38,  # "
                                Tq=4.38,  # "
                                Kr=Inf,  # molecs/cell, non-q autoinhib threshold
                                Kt=Inf,  # "
                                Km=Inf,  # "
                                Kq=152219,  # molecs/cell, q autoinhib threshold
                                hr=1,  # , non-q autoinhib Hill coeff
                                ht=1,  # "
                                hm=1,  # "
                                hq=4,  # , q autoinhib Hill coeff
                                kb=1,  # cell/(min molecs),  mRNA-ribo binding rate
                                ku=1,  # 1/min, mRNA-ribo unbinding rate
                                M=1e8,  # aa, total cell mass
                                kcm=.00599  # min/uM, chloramphenicol binding rate
                                ), ...)
{

with(as.list(c(parms, init)), {

## NOTE: equation numbers in parentheses (x)
## refer to supplemental material text, while
## numbers in square brackets [x] refer to
## equations in the main text

## transcription rate for proteins x
## NOTE: Kx=Inf, except for proteins q
ox < - function(wx, Tx, a, mx, Kx, hx) {
wx * a / (Tx+a) * 1 / (1+ (mx / Kx) ^ hx)}  # (10) * (11), with Kx!=Inf for proteins q
## translation elongation rate
## NOTE: in SBML: Kgamma = gmax/Kp = 1260/180.137803092828 = 6.994645
## see text for eq. (9) for derivation
g < - function(a) {gmax * a / (Kg+a)}  # (9)
## translation rate for proteins x
vx < - function(cx, nx, a) {g(a) * cx / nx}  # (9)

## substrate import and catabolic rates
vimp < - function(et, s) {et * vt * s / (Kt+s )}  # (7)
vcat < - function(em, si) {em * vm * si / (Km+si)}  # (7)

## protein classes
x < - c("t", "m", "q", "r")

## ODEs

## CULTURE LEVEL
## "We define the total mass of the cell as the total protein mass
##  (including bound ribosomes):
##  M = sum_x(nx*cx) + n_r * sum_x(cx) # (12) "
## "M is the mass of a mid-log phase cell"
## "for simulations we assume that eq. (14), with M=Ms,
##  also holds away from steady state"
## GROWTH RATE
Rt < - (cr + ct + cm + cq)  # total translating ribosomes
mu < - g(a) * Rt / M  # (14) "at steady state"
dN = mu * N                - dn * N  # [14]
ds = kin - vimp(et, s) * N - dn * s  # [15]

## CELL LEVEL
## internal substrate
dsi < - vimp(et, s) - vcat(em, si) - mu * si  # (1)

## energy equivalent (eg. a=ATP)
da < - ns * vcat(em, si) - mu * a
for ( i in 1:length(x) ) {  # all translation costs!
    cx < - get(paste0("c", x[i]), mode="numeric")
nx < - get(paste0("n", x[i]), mode="numeric")
da < - da - nx * vx(cx, nx, a)
}  # (2)

## free ribosomes
der < - vx(cr, nr, a) - mu * er
for (i in 1:length(x) ) {
cx < - get(paste0("c", x[i]), mode="numeric")
nx < - get(paste0("n", x[i]), mode="numeric")
mx < - get(paste0("m", x[i]), mode="numeric")
der < - der + vx(cx, nx, a) - kb * mx * er + ku * cx
}  # (3)

## protein X ODE
det < - vx(ct, nt, a) - mu * et  # (4)
dem < - vx(cm, nm, a) - mu * em  # (4)
deq < - vx(cq, nq, a) - mu * eq  # (4)

## mRNA X ODE
dmx < - function(mx, cx, nx, a, wx, Tx, Kx, hx) {
ox(wx, Tx, a, mx, Kx, hx) + ku * cx + vx(cx, nx, a) -
kb * mx * er - dm * mx - mu * mx
}  # (5)
dmt < - dmx(mt, ct, nt, a, wt, Tt, Kt, ht)
dmm < - dmx(mm, cm, nm, a, wm, Tm, Km, hm)
dmq < - dmx(mq, cq, nq, a, wq, Tq, Kq, hq)
dmr < - dmx(mr, cr, nr, a, wr, Tr, Kr, hr)

## mRNA:ribosome X ODE
dcx < - function(mx, cx, nx)  {
kb * mx * er - ku * cx - vx(cx, nx, a) - mu * cx
}  # (6)
dct < - dcx(mt, ct, nt)
dcm < - dcx(mm, cm, nm)
dcq < - dcx(mq, cq, nq)
dcr < - dcx(mr, cr, nr)


list(c(dN, ds, dsi, da,
det, dem, deq, der,
dct, dcm, dcq, dcr,
dmt, dmm, dmq, dmr))  # , log_y=unname(log(y))))
})

}



# ' Simple weisse growth model
# ' @family growth models
# ' @rdname grow_weisse
# ' @export
grow_weisse < - function(time, parms, ...)
{
    out < - deSolve:: ode(init, time, growthmodels::ode_weisse, parms = parms)
}
## attach names of parameters as attributes
attr(grow_weisse, "fname") < - c("grow_weisse")
attr(grow_weisse, "pnames") < - c("s", "y", "phi", "sin", "mumax", "K", "Y")


class(grow_weisse) < - c("growthmodel", "function")

'''

time = np.linspace(0, 30, 100)
#INPUT_VALUESdsi_dt = np.array(changeValues(time, INPUT, PAR))
# INPUT_VALUESdsi_dt = changeValues(INPUT_VALUES, PAR)
#print("INPUT DSI: ", INPUT_VALUESdsi_dt)
#print("LEN INPUT DSI: " , len(INPUT_VALUESdsi_dt))

"""

plot

"""

results = timeCourse(time, INPUT_VALUES)

names = ['s','N','si', 'a','et', 'em', 'eq','er', 'ct', 'cm', 'cq', 'cr','mt', 'mm', 'mq', 'mr']
#phi, sin, mumax, K, Y
#parms = [0, 0, 1, 1, 0.5]
#s=1, y=0.01
#y0 = [1, 0.01]
#out = timeCourse(time, INPUT_VALUES)

plt.title('Cell model with parameters from the paper', size = 20)
plt.xlabel('Time', size = 20)
plt.ylabel('Concentration', size = 20)
plt.xticks(size = 15)
plt.yticks(size = 15)
lines = plt.plot(time, results)
plt.legend(lines[:16], names, prop = {'size': 12}, loc = 'upper left', frameon=True, ncol=2)
plt.show()




