


import modelbase
import modelbase.ratelaws as rl


import numpy as np


defaultParameters = { 's': 1e4,  # external nutrient [ molecs ]
                      'dm': 0.1,  # mRNA degradation rate [1/min ]
                      'ns': 0.5,  # nutrient efficiency [ none ]
                      'nr': 7459,  # ribosome length [ aa / molecs ]
                      'nt': 300,  # length of non-ribosomal proteins [ aa / molecs ]
                      'nm': 300,  # length of non-ribosomal proteins [ aa / molecs ]
                      'nq': 300,  # length of non-ribosomal proteins [ aa / molecs ]
                      'gammamax': 1260,  # max. transl. elongation rate [ aa / min molecs ]
                      'Kgamma': 7,  # transl. elongation threshold [ molecs / cell ]
                      'vt': 726,  # max. nutrient import rate [ 1/min ]
                      'Kt': 1000,  # nutrient import threshold [ molecs ]
                      'vm': 5800,  # max. enzymatic rate [ 1/min ]
                      'Km': 1000,  # enzymatic threshold
                      'wr': 930,  # max. ribosome transcription rate [ molecs / min cell ]
                      'wt': 4.14,  # max. enzyme transcription rate [ molecs / min cell ]
                      'wm': 4.14,  # max. enzyme transcription rate [ molecs / min cell ]
                      'wq': 948.93,  # max. q-transcription rate
                      'thetar': 426.87,  # ribosome transcription threshold [ molecs / cell ]
                      'thetat': 4.38,  # non-ribosomal transcription threshold [ molecs / cell ]
                      'thetam': 4.38,  # non-ribosomal transcription threshold [ molecs / cell ]
                      'thetaq': 4.38,  # non-ribosomal transcription threshold [ molecs / cell ]
                      'Kq': 152219,  # q-autoinhibition threshold [ molecs / cell ]
                      'hq': 4,  # q-autoinhibition Hill coeff. [ none ]
                      'kb': 1,  # mRNA-ribosome binding rate
                      'ku': 1,  # mRNA-ribosome unbinding rate
                      'M': 1e8,  # total cell mass
                      'kcm': 0.00599, # chloramphenicol-binding rate
}


indexx = ['r','t','m','q']


class WeisseBaseModel(modelbase.Model):
    '''
    class to define the basic model proposed by Weiße et al, 2015, PNAS
    '''

    def __init__(self, pars = {}):
        '''
        constructor method for basic Weisse model
        '''

        defPars = defaultParameters.copy()

        super(WeisseBaseModel,self).__init__(defPars)

        self.par.update(pars)

        self.add_cpds(['si','a']) # si: internal substrate, a: energy
        self.add_cpds(['m'+x for x in indexx]) # mRNA
        self.add_cpds(['c'+x for x in indexx]) # RNA/protein complexes
        self.add_cpds(['p'+x for x in indexx]) # proteins

        def gamma(par, a):
            return np.array([par.gammamax * a[0] / (par.Kgamma + a[0])])
        self.add_algebraicModule(gamma, 'gamma', ['a'], ['gamma'])
        
        def vimp(par, et):
            return et * par.vt * par.s / (par.Kt + par.s)
        self.add_reaction('vimp',vimp,{'si':1},'pt')

        def vcat(par, em, si):
            return em * par.vm * si / (par.Km + si)
        self.add_reaction('vcat',vcat,{'si':-1,'a':self.par.ns},'pm','si')

        def makevx(x):
            def vx(par, cx, gamma):
                return cx * gamma / getattr(par,'n'+x)
            return vx
        
        for x in indexx: # all vx's
            rname = 'v'+x
            cx = 'c'+x
            nx = 'n'+x
            mx = 'm'+x
            px = 'p'+x
            stDict = {'a':-getattr(self.par,nx), 'pr': 1, px: 1, mx: 1, cx: -1}
            if x == 'r':
                stDict['pr'] = 2
            vx = makevx(x)
            self.add_reaction(rname, vx, stDict, cx, 'gamma')
            
        def complexbu(par, r, mx, cx):
            return par.kb * r * mx - par.ku * cx

        for x in indexx: # all complex binding/unbinding reactions
            rname = 'cbu'+x
            mx = 'm'+x
            cx = 'c'+x
            stDict = {'pr':-1, mx: -1, cx: 1}
            self.add_reaction(rname, complexbu, stDict, 'pr', mx, cx)

        def inhibitionq(par, q):
            return np.array([1 / (1 + (q[0]/par.Kq) ** par.hq)])
        self.add_algebraicModule(inhibitionq, 'Iq', ['pq'], ['Iq'])
                
        def makeomegax(x):
            def omegax(par, a):
                return getattr(par, 'w'+x) * a / (getattr(par, 'theta'+x) + a)
            return omegax

        def degmx(par, mx):
            return par.dm * mx
            
        for x in indexx: # all transcription rates omegax
            rname = 'omega'+x
            mx = 'm'+x
            stDict = {mx: 1}
            if x == 'q':
                def omegax(par, a, Iq):
                    return Iq * par.wq * a / (par.thetaq + a)
                self.add_reaction(rname, omegax, stDict, 'a', 'Iq')
            else:
                omegax = makeomegax(x)
                self.add_reaction(rname, omegax, stDict, 'a')
            
            stDictDeg = {mx: -1}
            self.add_reaction('degm'+x, degmx, stDictDeg, mx)

        def stGrowth(par, y):
            return np.array([y[0] * y[1:].sum() / par.M])
        self.add_algebraicModule(stGrowth, 'lambda', ['gamma', 'cr', 'ct', 'cm', 'cq'], ['lambda'])

        def dilution(par, lam, x):
            return lam * x

        for cpd in self.cpdNames:
            rname = 'dilution_'+cpd
            stDict = {cpd: -1}
            self.add_reaction(rname, dilution, stDict, 'lambda', cpd)

        
