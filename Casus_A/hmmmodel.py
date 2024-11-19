#!/usr/bin/python3

class HiddenMarkovModel:

    startprob_ = 0
    transmat_ = 0
    emmisionprob_ = 0

    def __init__(self, n_components, n_features):
        print("init")

    def __str__(self):
        print("To string")

    def __repr__(self):
        print("To repr")

    def sample(self):
        print("To sample")
