#!/usr/bin/python3
import math
import random

import numpy as np


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


    def score(self, emissions, states):
        self.not_empty()

        logprob = np.log(self.startprob_[states[0]] * self.emissionprob_[states[0]][emissions[0]])

        for i in range(len(emissions) - 1):
            transition_prob = self.transmat_[states[i]][states[i + 1]]
            emission_prob = self.emissionprob_[states[i + 1]][emissions[i + 1]]

            logprob += np.log(transition_prob) + np.log(emission_prob)

        return logprob

    def not_empty(self):
        assert self.startprob_ is not None
        assert self.transmat_ is not None
        assert self.emissionprob_ is not None
        return True

    def predict(self, emissions):
        self.not_empty()
        if emissions is None or emissions == []:
            return []

        n_rounds = len(emissions)

        # opslaan voor viterbi
        viterbi = [[-math.inf] * self.n_components for _ in range(n_rounds)]
        previous_max = [[None] * self.n_components for _ in range(n_rounds)]

        # beurt 1
        curr_emission = emissions[0]
        for i in range(self.n_components):
            viterbi[0][i] = math.log(self.startprob_[i]) + math.log(self.emissionprob_[i][curr_emission]) # fix me

        # alle andere beurten
        for i in range(1, n_rounds):
            curr_emission = emissions[i]
            for t in range(self.n_components):
                max_prob = -math.inf
                max_prev_state = 0
                for p in range(self.n_components):
                    previous = viterbi[i-1][p]
                    trans_prob = math.log(max(self.transmat_[p][t], 1e-10))
                    emission_prob = math.log(max(self.emissionprob_[t][curr_emission], 1e-10))
                    prob = previous + trans_prob + emission_prob
                    if prob > max_prob:
                        max_prob = prob
                        max_prev_state = p
                viterbi[i][t] = max_prob
                previous_max[i][t] = max_prev_state

        # bepaal laatste tafel
        final_state = max(range(self.n_components), key=lambda t: viterbi[-1][t])

        # achterstevoren de states terug gaan
        best_path = [final_state + 1]
        for i in range(n_rounds - 1, 0, -1):
            final_state = previous_max[i][final_state]
            best_path.insert(0, final_state + 1)

        return best_path
