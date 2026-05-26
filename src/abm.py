#!/usr/bin/python3
# -*- coding: utf-8 -*-
import random
import numpy as np


def draw():
    """Draw a uniform random number in [0, 1)."""
    return np.random.random_sample()


class HealthCareWorker:
    """
    Attributes:
        strain_set (set): Carried pathogens as (species, phenotype) tuples.
        record (dict):    History of status changes.
    """

    def __init__(self):
        self.strain_set = set()
        self.record = {}

    def recordStatus(self, reason, ts):
        idx = 0 if not self.record else max(self.record.keys()) + 1
        self.record[idx] = {
            'strain_set': sorted(list(self.strain_set)),
            'ts': ts,
            'reason': reason,
        }


class Patient:
    """
    Attributes:
        status (str):       'C' (colonized/carrier) or 'I' (infected).
        species (str):      ESKAPE species key (see config.PATHOGENS) or 'none'.
        phenotype (str):    Resistance phenotype string or 'none'.
        drug_agent (str):   Treating agent name (see config.AGENT_LEVEL), or None if not on treatment.
        treat_time (int):   Hours since treatment started; None if not infected.
        convt_time (int):   Hour at which dominant pathogen changed; None if unchanged.
        lab_result (tuple): (species, phenotype) confirmed by lab, or None.
        time_inICU (int):   Total hours in ICU.
        super_infe (bool):  True if patient has experienced super-infection.
        infct_flag (bool):  True if patient is destined to develop infection.
        infct_time (int):   Simulation day on which infection will manifest.
    """

    def __init__(self, args, name, current_day, species='none', phenotype='none'):
        self.name = name
        self.status = 'C'
        self.species = species
        self.phenotype = phenotype
        self.drug_agent = None
        self.treat_time = None
        self.total_treat_time = 0
        self.convt_time = None
        self.lab_result = None
        self.time_inICU = 0
        self.super_infe = False
        self.infct_flag = False
        self.infct_time = None
        self.record = {}

        self._set_infection_time(args, current_day)

    # ── Infection-time helpers ────────────────────────────────────────────────

    def _set_infection_time(self, args, current_day, time_interval=5):
        """Stochastically decide if/when colonization will progress to infection."""
        if self.species == 'none':
            return
        sigma = self._get_sigma(args)
        if draw() < sigma:
            self.infct_flag = True
            self.infct_time = current_day + random.randint(0, time_interval)

    def reset_infection_time(self, args, current_day, time_interval=5):
        """Re-determine infection time after the dominant pathogen changes."""
        # Reset existing flag
        self.infct_flag = False
        self.infct_time = None
        if self.species == 'none':
            return
        sigma = self._get_sigma(args)
        if draw() < sigma:
            self.infct_flag = True
            self.infct_time = current_day + random.randint(0, time_interval)

    def _get_sigma(self, args):
        """Look up sigma for current species/phenotype from config."""
        from config import PATHOGENS
        cfg = PATHOGENS.get(self.species)
        if cfg:
            return cfg['sigma'].get(self.phenotype, 0.10)
        return 0.0  # uncolonized (species='none') never progresses to infection

    # ── Record-keeping ────────────────────────────────────────────────────────

    def recordStatus(self, reason, currentT):
        idx = 0 if not self.record else max(self.record.keys()) + 1
        self.record[idx] = {
            'reason': reason,
            'currentTS': currentT,
            'attributes': {
                'name':       self.name,
                'status':     self.status,
                'species':    self.species,
                'phenotype':  self.phenotype,
                'drug_agent':        self.drug_agent,
                'total_treat_time':  self.total_treat_time,
                'lab_result':        self.lab_result,
                'time_inICU': self.time_inICU,
                'infct_time': self.infct_time,
                'treat_time': self.treat_time,
                'convt_time': self.convt_time,
                'super_infe': self.super_infe,
            },
        }
