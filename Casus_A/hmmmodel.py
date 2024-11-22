#!/usr/bin/python3
import random
from logging import logMultiprocessing

import numpy as np
from fontTools.misc.cython import returns
from orca.orca import start
from orca.settings import startingProfile


class HiddenMarkovModel:

    def __init__(self, n_components, n_features):
        self.n_components = n_components
        self.n_features = n_features
        self.startprob_ = None
        self.transmat_ = None
        self.emissionprob_ = None

    def __str__(self):
        return f"HiddenMarkovModel(n_components={self.n_components}, n_features={self.n_features})"

    def __repr__(self):
        return (f"HiddenMarkovModel(startprob_={self.startprob_}, "
                f"transmat_={self.transmat_}, emissionprob_={self.emissionprob_})")

    def sample(self, size):

        assert size > 0
        self.not_empty()

        my_components = range(self.n_components)
        my_features = range(self.n_features)
        current_state = random.choices(my_components, self.startprob_)[0]
        states = [current_state]
        current_emission = random.choices(my_features, self.emissionprob_[current_state])[0]
        emissions = [current_emission]

        for i in range(size - 1):
            current_state = random.choices(my_components, self.transmat_[current_state])[0]
            states.append(current_state)
            current_emission = random.choices(my_features, self.emissionprob_[current_state])[0]
            emissions.append(current_emission)

        return states, emissions

    def not_empty(self):
        assert self.startprob_ is not None
        assert self.transmat_ is not None
        assert self.emissionprob_ is not None
        return True

    def score(self, emissions, states):
        self.not_empty()
        logprob = np.log(self.startprob_[states[0]] * self.emissionprob_[self.startprob_[states[0]]][emissions[0]]) # Start
        for i in range(len(emissions) - 1):
            # transmat kans
            self.transmat_[states[i]][0]
            # emission prob
            self.emissionprob_[states[i+1]][emissions[i+1]]
            logprob += np.log()

        return logprob


