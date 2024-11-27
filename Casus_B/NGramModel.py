#!/usr/bin/python3
import numpy as np
from scipy.constants import value


class NGramModel:

    def __init__(self, key):
        self.skip = key
        self.model_ = {}

    def fit(self, text):
        if self.skip <1:
            pass
        # loop once
        for i in range(len(text)):
            if (i + self.skip) < len(text):
                key = text[i:i+self.skip]
                value = text[i + self.skip]
                self.__add_to_model__(key, value)

    def __add_to_model__(self, key, value):
        if key in self.model_:
            if value in self.model_[key]:
                self.model_[key][value] += 1
            else:
                self.model_[key][value] = 1
        else:
            self.model_[key] = {value: 1}

    def predict_proba(self, key):
        return self.model_[key] if key in self.model_ else []

    def predict(self, seed, length):
        result = seed

        for _ in range(length - len(seed)):
            key = result[-self.skip:]
            cases = self.predict_proba(key)
            possabilities= [k for k in cases.keys()]
            chances = [c / sum(cases.values()) for c in cases.values()]
            result += np.random.choice(possabilities, p=chances)

        return result