# -*- coding: utf-8 -*-

from __future__ import print_function

import numpy as np
from astroquery.lamda import Lamda
from astroquery.hitran import read_hitran_file, download_hitran, cache_location
from scipy import constants
import itertools
import os

def omega(theta):
    """Solid angle subtended by the Sun

    Parameters
    ----------
    theta : float
        Angular diameter of the Sun

    Returns
    -------
    omega : float
        The resulting solid angle for the Sun
        https://en.wikipedia.org/wiki/Solid_angle
    """
    omega = 2*np.pi*(1-np.cos(theta/2))
    return omega


def zlamda(levels, temp):
    """Partition function from LAMDA data

    Parameters
    ----------
    levels : pandas DataFrame
        LAMDA levels table

    temp : float
        Gas temperature

    Returns
    -------
    zlamda : float
        Partition function
    """
    pop = levels['Weight']*np.exp(-constants.h*constants.c*1e2*levels['Energy']/constants.k/temp)
    return pop.sum()


def glu(sigma, gu, gl, El, Z):
    """Pumping rate from lower to upper level glu
    divided by the Einstein coefficient Aul
    (equation A7 from Crovisier & Encrenaz 1983)

    Parameters
    ----------
    sigma : float
        wavenumber in units of cm-1
    gu : float
        upper level degeneracy
    gl : float
        lower level degeneracy
    El : float
        energy of lower level in cm-1
    temp : float
        gas temperature
    Z : float
        partition function

    Returns
    -------
    glu : float
        This is actually glu/Aul
    """
    # solid angle subtended by the Sun at 1 AU
    om = omega(9.35e-3)
    # black-body Sun temperature in K
    Tbb = 5778
    # relative population of the lower level
    pl = gl * np.exp(-constants.h*constants.c*El*1e2/constants.k/50)/Z
    glu = om/4./np.pi*gu/gl*pl/\
        (np.exp(constants.h*constants.c*sigma*1e2/constants.k/Tbb)-1)
    return glu


def hitran_bands(df):
    """Print all bands in HITRAN transitions

    Parameters
    ----------
    df : pandas DataFrame
        HITRAN transitions returned by read_hitran_file
    """
    bands = df.global_upper_quanta.unique()
    for b in bands:
        nu = df[(df.global_upper_quanta.str.contains(b))
                & (df.global_lower_quanta.str.contains(b'GROUND'))
                ].nu.mean()
        print(b, nu)


def lamda_to_hitran(level, mol):
    """Convert levels from Lamda to HITRAN notation

    Parameters
    ----------
    level : str
        LAMDA level notation
    mol : str
        molecule name
    """
    import re
    if mol == "aCH3OH":
        m = re.match('(\d+)_(-?)(\d+)', level)
        if m.group(2) == '-':
            quanta = m.group(1).rjust(3)+m.group(3).rjust(3)+'  A-     '
        else:
            quanta = m.group(1).rjust(3)+m.group(3).rjust(3)+'  A+     '
    elif mol == "eCH3OH":
        m = re.match('(\d+)_(-?)(\d+)', level)
        if m.group(2) == '-':
            quanta = m.group(1).rjust(3)+m.group(3).rjust(3)+'  E2     '
        else:
            quanta = m.group(1).rjust(3)+m.group(3).rjust(3)+'  E1     '
    elif mol == "H2O":
        m = re.match('(\d+)_(\d+)_(\d+)', level)
        quanta = m.group(1).rjust(3)+m.group(2).rjust(3)+m.group(3).rjust(3)+'      '
    return quanta.encode('utf-8')


# Read path to HITRAN and LAMDA data from environment variable
HITRAN_DATA = os.environ.get(
    "HITRAN_DATA",
    cache_location
)

LAMDA_DATA = os.environ.get(
    "LAMDA_DATA",
    Lamda.cache_location
)

# names in LAMDA database
lamda = {'H2O': 'oh2o@daniel',
        'HDO': 'hdo',
        'aCH3OH': 'a-ch3oh',
        'eCH3OH': 'e-ch3oh'}
# identifiers in HITRAN database
hitran_ids = {'H2O': (1,1),
             'HDO': (1,4),
             'aCH3OH': (1,1),
             'eCH3OH': (1,1)}
# name of ground vibrational state in LAMDA database
ground_state= {'H2O': b'0 0 0',
               'HDO': b'0 0 0',
               'aCH3OH': b'GROUND',
               'eCH3OH': b'GROUND'}
# subset of bands with transitions to ground state
excitation_band = {'H2O': ['0 0 1', '0 1 0', '1 0 0', '1 0 1', '0 1 1'],
                   'HDO': b'.*',
                   'aCH3OH': ['4V12', 'V12', 'V7', 'V8'],
                   'eCH3OH': ['3V12', 'V12', 'V7', 'V8']}
# species identifier
species = {'H2O': b'',
           'HDO': b'',
           'aCH3OH': b'A',
           'eCH3OH': b'E'}

# print(trans[~trans.global_lower_quanta.str.contains(b'0 0 0')])

class gfactor(object):
    """
    Calculate effective infrared pumping rates excited by solar radiation as a
    black body radiation field
    """

    def __init__(self, mol, nlev=False):
        """
        pumping rates excited by solar radiation ignoring hot band cascades
        """
        self.mol = mol
        # download all transitions of the molecule between the wavenumbers of 0 and 20000 cm-1
        download_hitran(*hitran_ids[mol], 0, 20000)
        self.tbl = read_hitran_file(os.path.join(HITRAN_DATA,
                                            '{0}.data'.format(''.join(x for x in mol if not x.islower()))),
                                            # formatfile=os.path.join(HITRAN_DATA, 'readme.txt')
                                            ).to_pandas()

        collrates,radtransitions,enlevels = Lamda.query(mol=lamda[mol])

        levels = enlevels.to_pandas()

        # create new column for quantum numbers using HITRAN notation
        levels.loc[:,'local_quanta'] = levels['J'].apply(lamda_to_hitran, args=(mol,))

        # create new column for relative population

        lamda_levels = levels.Energy.values
        if mol == 'aCH3OH':
            # Fix a bug in energy levels for aCH3OH in HITRAN 2012
            lamda_levels += 128.1069

        self.trans = self.tbl[
                    # select species
                    self.tbl["local_lower_quanta"].str.contains(species[mol]) &
                    # select and bands with vibrational transitions to ground state
                    self.tbl['global_upper_quanta'].isin([i.rjust(15).encode('utf-8') for i in excitation_band[mol]])
                    # mol['global_lower_quanta'].str.contains(ground_state[mol]) &
                    # ~mol['global_upper_quanta'].str.contains(ground_state[mol])
                    # & np.isclose(mol["elower"].values[:, None], lamda_levels, atol=.01).any(axis=1)
                    # & mol["local_lower_quanta"].isin(levels["local_quanta"])
                    ]

        # hitran_bands(trans, mol)

        # hitran_levels = self.tbl[self.tbl['global_lower_quanta'].str.contains(ground_state[mol])]["local_lower_quanta"].unique()
        # nlev = len(hitran_levels)

        if not nlev:
            nlev = len(levels)
        self.gcube = np.zeros((nlev, nlev))

        # group transitions by common upper level
        grouped = self.trans.groupby(['global_upper_quanta', 'local_upper_quanta'])

        # select upper level from any vibrational state different from ground
        # nu    u----
        # --------||\----------
        #         || \
        #         ||  \
        #      Au1||Au2\Au3
        #         ||    \
        # nu'     ||   3---
        # --------||-----------
        #         |-- 2 (up)
        # nu''   --- 1 (lo)
        # ---------------------
        Z = zlamda(levels[:nlev], 50)
        for _, group in grouped:
            # transitions that go to the lamda levels in the ground vibrational state
            ground = group[
                group['global_lower_quanta'].str.contains(ground_state[mol]) &
                group['local_lower_quanta'].isin(levels["local_quanta"][:nlev])
                ]
            if len(ground) >= 2:
                # transitions from k to ground level
                Asum = group['a'].sum()
                # add a new column with glu values (equation A7 from Crovisier & Encrenaz 1983)
                ground.loc[:,'glu'] = glu(ground['nu'], ground['gp'], ground['gpp'], ground['elower'], Z)
                # combination of pairs of possible levels in ground vibrational state
                for lo, up in itertools.combinations(ground['local_lower_quanta'], 2):
                    trans_lo = ground[ground['local_lower_quanta'] == lo]
                    trans_up = ground[ground['local_lower_quanta'] == up]
                    # multiply Aeins for the pair
                    Aprod = ground[ground['local_lower_quanta'].isin([lo, up])]['a'].product()
                    # g12 += g1u * Au2/(Au1 + Au2 + Au3)
                    if lo != up:
                        i = levels[levels["local_quanta"] == lo]["Level"].values[0]
                        j = levels[levels["local_quanta"] == up]["Level"].values[0]
                        self.gcube[i-1, j-1] += trans_lo['glu'].values[0]*Aprod/Asum
                        self.gcube[j-1, i-1] += trans_up['glu'].values[0]*Aprod/Asum


    def hotbands(self):
        """
        Add contribution from hotbands to effective pumping rates
        """
        # selection for hot bands
        trans = self.tbl[
                    self.tbl['global_lower_quanta'].str.contains(ground_state[self.mol])
                    # & ~self.tbl['global_upper_quanta'].str.contains(ground_state[self.mol])
                    & self.tbl['global_upper_quanta'].isin([i.rjust(15).encode('utf-8') for i in excitation_band[self.mol]])
                    & self.tbl["local_lower_quanta"].isin(levels["local_quanta"])
                    ]

        trans_hb = self.tbl[self.tbl['global_upper_quanta'].isin(trans['global_upper_quanta'])
                    & self.tbl['local_upper_quanta'].isin(trans['local_upper_quanta'])]
        grouped2 = trans.groupby(['global_upper_quanta', 'local_upper_quanta'])

        # hot bands
        for k, hb in trans.iterrows():
            i = levels[levels["local_quanta"] == hb['local_lower_quanta']]["Level"].values[0]
            # all transitions from upper level in vibrational state
            ground = self.tbl[self.tbl['global_upper_quanta'].str.contains(hb['global_upper_quanta'], regex=False)
                        & self.tbl['local_upper_quanta'].str.contains(hb['local_upper_quanta'], regex=False)]
            # transitions to level in excited vibrational state
            mid = ground[~ground['global_lower_quanta'].str.contains(ground_state[self.mol])]
            if len(mid) == 1:
                ground2 = self.tbl[self.tbl['global_upper_quanta'].str.contains(mid['global_upper_quanta'].values[0])
                        & self.tbl['local_upper_quanta'].str.contains(mid['local_upper_quanta'].values[0], regex=False)]
                # hot band transitions that go to ground state
                ground3 = ground2[ground2['global_lower_quanta'].str.contains(ground_state[self.mol])
                                & ground2['local_lower_quanta'].isin(levels['local_quanta'])
                                ]
                for kk, g in ground3.iterrows():
                    j = levels[levels["local_quanta"] == g['local_lower_quanta']]["Level"].values[0]
                    if i != j:
                        Asum = ground['a'].sum()
                        Aprod = hb['a']*mid['a'].values[0]
                        Asum2 = ground2['a'].sum()
                        g_ik = glu(hb['nu'], hb['gp'], hb['gpp'])
                        gcube[i-1, j-1] += g_ik*Aprod/Asum*g['a']/Asum2

        for name, group in grouped:
            # groups that have transitions to the ground state
            ground = group[group['global_lower_quanta'].str.contains(ground_state[self.mol])
                    & group['local_lower_quanta'].isin(levels["local_quanta"])]
            if len(ground) >=1:
                for _, group2 in grouped:
                    upper1 = group2[group2['global_lower_quanta'].str.contains(group['global_upper_quanta'].values[0])
                        & group2['local_lower_quanta'].str.contains(group['local_upper_quanta'].values[0])]
                    upper2 = group2[group2['global_lower_quanta'].str.contains(ground_state[self.mol]) &
                        group2['local_lower_quanta'].isin(levels["local_quanta"])]
                    if len(upper1) & len(upper2):
                        Asum = group2['a'].sum()
                        Asum2 = ground['a'].sum()
                        for k, hb in upper2.iterrows():
                            i = levels[levels["local_quanta"] == hb['local_lower_quanta']]["Level"].values[0]
                            g_ik = glu(hb['nu'], hb['gp'], hb['gpp'])
                            for k, hbj in ground.iterrows():
                                j = levels[levels["local_quanta"] == hbj['local_lower_quanta']]["Level"].values[0]
                                Aprod = hb['a']*upper1['a'].values[0]
                                gcube[i-1, j-1] += g_ik*Aprod/Asum*hbj['a']/Asum2


    def scale(self, rh):
        """
        Scale infrared pumping rates by heliocentric distance
        """
        self.gcube /= rh**2
