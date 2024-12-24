# ! /User/bin/env python3
from random import random
import numpy as np


class BaseDifferential():
    """
    Main class containing the ODE solvers, fit, and __init__ functions for 
    the differential (tumor) equations.
    """
    def __init__(self, solver="runge-kutta", fitter="direct", solution="BIC"):
        self.params = {'c':0.10, "y0":0.50}
        self.solver = solver
        self.fitter = fitter
        self.solution = solution
        self.options = ["runge-kutta", "heun", "euler"]

    def __repr__(self):
        pass

    def __str__(self):
        return f"parameters = {self.params}"

    def formula(self, y, params=None):
        """
        This is the linear differential equation, used to calculate the
        growth. It is used by the ODE solvers, fit, and BIC or AICc calculations
        : param y: The y-value used to calculate the growth
        : param params: The parameters for the equation used to calculate the growth
        : return: The calculated growth (from the equation)
        """
        if y < 0:
            return 0
        if params is None:
            params = self.params
        c = params['c']
        return c

    def heun_method(self, t, Nsteps=None, params=None):
        """
        The Heun method uses the formula of a differential equation to determine the growth
        between two data points, with t/N_steps on the x-axis.  
        It calculates the average of the growth at the beginning and at the end.  
        : param t: The x value (int) being targeted.  
        : param Nsteps: The number of steps used to reach t (precision of the growth).  
        : return y: The y value corresponding to t
        """
        # Makes the Nsteps 100 times the t value + 1
        if Nsteps is None:
            Nsteps = int(abs(t)/0.01) + 1
        # uses the params from the module if non is given
        if params is None:
            params = self.params
        dt = t/Nsteps
        y = params["y0"]
        for _ in range(Nsteps):
            growth = self.formula(y, params)
            y2 = y + growth * dt
            # Calculate second growth with previous growth
            growth2 = self.formula(y2, params)
            # Use the average growth
            avg_growth = (growth+growth2) / 2.0
            y += avg_growth * dt

        return y

    def euler_method(self, t, Nsteps=None, params=None):
        """
        The Euler method uses the formula of a differential equation to determine the growth
        between two data points, with t/N_steps on the x-axis.  
        : param t: The x value (int) being targeted.  
        : param Nsteps: The number of steps used to reach t (precision of the growth).  
        : return y: The y value corresponding to t.
        """
        # Makes the Nsteps 100 times the t value + 1
        if Nsteps is None:
            Nsteps = int(abs(t)/0.01) + 1
        # uses the params from the module if non is given
        if params is None:
            params = self.params
        dt = t/Nsteps
        y = params["y0"]
        for _ in range(Nsteps):
            # Calculate the growth
            growth = self.formula(y, params) * dt
            y = y + growth

        return y

    def runge_kutta(self, t, Nsteps=None, params=None):
        """
        The runge_kutta method uses the formula of a differential equation to determine the growth
        between two data points, with t/N_steps on the x-axis.  
        It calculates the growth using calculated growths from different datapoint 
        : param t: The x value (int) being targeted.  
        : param Nsteps: The number of steps used to reach t (precision of the growth).  
        : return y: The y value corresponding to t
        """
        # Makes the Nsteps 100 times the t value + 1
        if Nsteps is None:
            Nsteps = int(abs(t)/0.01) + 1
        # uses the params from the module if non is given
        if params is None:
            params = self.params
        dt = t/Nsteps
        y = params["y0"]
        for _ in range(Nsteps):
            # Calculate multiple growths (using the previous growth)
            growth1 = self.formula(y, params)
            y2 = y + growth1 * dt * 0.5
            growth2 = self.formula(y2, params)
            y3 = y + growth2  * dt * 0.5
            growth3 = self.formula(y3, params)
            y4 = y + growth3* dt
            growth4 = self.formula(y4, params)

            # Calculate the growth with growth from certain points
            growth = (growth1 + growth2 * 2 + growth3 * 2 + growth4) / 6.0
            y += growth * dt
        return y
    
    def mean_squared_error(self, ts, ys, params):
        """
        mean_squared_error calculates how much the y-values from the calculated graph
        differ from the actual y-values. This suggets the accuracy of the calculated graph
        : param ts: list with all the x-values connected to the y-values (ys)
        : param ys: list with all the y-values connected to the x-values (ts)
        : param params: all the parameters used to made the calculated graph
        : return: The mean squared error value of the calculated graph
        """
        N = len(ts)
        total = 0.0
        for datapoint in range(N):
            if self.solver == "runge-kutta":
                prediction = self.runge_kutta(ts[datapoint], Nsteps=None, params=params)
            elif self.solver == "heun":
                prediction = self.heun_method(ts[datapoint], Nsteps=None, params=params)
            else:
                prediction = self.euler_method(ts[datapoint], Nsteps=None, params=params)
            error = ys[datapoint] - prediction
            total += error * error
        return total / N

    def direct(self, ts, ys):
        """
        Calculates the parameters that fit the y-values the best on the selected
        differential equation.
        The Hookes and jeeves (direct search) to find the best parameters
        : param ts: list with all the x-values connected to the y-values (ys)
        : param ys: list with all the y-values connected to the x-values (ts)
        : return: It saves the best parameters, so other functions can be used on
        these parameters
        """
        # Value error when solver argument is not correct
        if self.solver not in self.options:
            raise ValueError("Solver must be: 'euler', 'heun' or 'runge-kutta'")

        params = self.params
        # The step size is set at 0.3 for all the params
        delta = {key:0.3 for key in self.params}
        mse = self.mean_squared_error(ts, ys, self.params)

        # Only when the highets step size is close to zero
        while max(abs(delta) for delta in delta.values()) > 1e-6:
            # For all the parameters of the differential equation
            for key in params:
                new_param = params.copy()
                new_param[key] = params[key] + delta[key]
                # MSE is calculated for positive delta + param
                new_mse = self.mean_squared_error(ts, ys, new_param)
                if new_mse < mse:
                    # The step size is made bigger in the direction
                    mse = new_mse
                    params = new_param
                    delta[key] *= 1.2
                    continue
                new_param[key] = params[key] - delta[key]
                new_mse = self.mean_squared_error(ts, ys, new_param)
                if new_mse < mse:
                    # The step size is made negative for the opposite direction
                    mse = new_mse
                    params = new_param
                    delta[key] *= -1.0
                    continue
                # The step size is made smaller when both direction have no improvement
                delta[key] *= 0.5

        # The new parameters are saved in the model
        self.params = params

    def random_fit(self, ts, ys):
        """
        Calculates the parameters that fit the y-values the best on the selected
        differential equation.
        The Monte Carlo (random search) to find the best parameters
        : param ts: list with all the x-values connected to the y-values (ys)
        : param ys: list with all the y-values connected to the x-values (ts)
        : return: It saves the best parameters, so other functions can be used on
        these parameters
        """
        tries = 0
        params = self.params
        mse = self.mean_squared_error(ts, ys, params)
        # Only when there is no improvement within the amount of tries
        while tries < 1000:
            # The new parameters are chosen at random
            new_params = {key: params[key] + random() - 0.5 for key in params}
            new_mse = self.mean_squared_error(ts, ys, new_params)
            if new_mse < mse:
                # parameters and mse are saved
                mse = new_mse
                params = new_params
                # tries are reset when there's improvement
                tries = 0
            tries += 1

        self.params = params

    def fit(self, ts, ys):
        """
        Uses the chosen fit model by the user. When the user gives a non existing
        fit model, a value error is given to the user
        : param ts: list with all the x-values connected to the y-values (ys)
        : param ys: list with all the y-values connected to the x-values (ts)
        : return: The parameters for the best model are saved after fitting, or
        ValueError when the wrong fitter is chosen
        """
        if self.fitter == "direct":
            self.direct(ts, ys)
        elif self.fitter == "random":
            self.random_fit(ts, ys)
        else:
            raise ValueError("Fitter must be either: 'direct' or 'random'")

    def quality(self, ts, ys):
        """
        Return either a BIC or AICc value for the fitted graph with the x-
        and y-values. The BIC and AICc values can be used to choose the best
        differential equantion for the ts and ys data
        : param ts: list with all the x-values connected to the y-values (ys)
        : param ys: list with all the y-values connected to the x-values (ts)
        : return: either the BIC or the AICc score for the graph. If the chosen
        solution doesn't exist, it will be mentioned to the user
        """
        if self.solution == "BIC":
            return self.bic(ts, ys)
        if self.solution == "AICC":
            return self.aicc(ts, ys)
        else:
            return "Solution must be either: 'BIC' or 'AICC'"

    def bic(self, ts, ys):
        """
        Calculates the BIC value for the fitted graph with the x-
        and y-values. This value can be used to choose the best
        differential equation for the ts and ys data
        : param ts: list with all the x-values connected to the y-values (ys)
        : param ys: list with all the y-values connected to the x-values (ts)
        : return: The BIC value belonging to the fitted graph
        """
        mse = self.mean_squared_error(ts, ys, self.params)
        n = len(ts)
        k = len(self.params)
        return n * np.log(mse) + k * np.log(n)

    def aicc(self, ts, ys):
        """
        Calculates the AICc value for the fitted graph with the x-
        and y-values. This value can be used to choose the best
        differential equation for the ts and ys data
        : param ts: list with all the x-values connected to the y-values (ys)
        : param ys: list with all the y-values connected to the x-values (ts)
        : return: The AICc value belonging to the fitted graph 
        """
        mse = self.mean_squared_error(ts, ys, self.params)
        n = len(ts)
        k = len(self.params)
        return (n * np.log(mse) + 2 * k)  + (2 * k * (k+1)) / (n - k - 1)


class Linear(BaseDifferential):
    """
    This class is used to run the functions of the `BaseDifferential` class with
    the linear differential equation.
    """
    def __str__(self):
        return f"{self.params['c']: .2f}"

    def __repr__(self):
        return f"Linear(solver={self.solver}, fitter={self.fitter}, solution={self.solution})"

    def formula(self, y, params=None):
        """
        This is the linear differential equation, used to calculate the
        growth. It is used by the ODE solvers, fit, and BIC or AICc calculations
        : param y: The y-value used to calculate the growth
        : param params: The parameters for the equation used to calculate the growth
        : return: The calculated growth (from the equation)
        """
        if y < 0:
            return 0
        if params is None:
            params = self.params
        c = params['c']
        return c


class Mendelsohn(BaseDifferential):
    """
    This class is used to run the functions of the `BaseDifferential` class with
    the Mendelsohn differential equation.
    """
    def __init__(self, solver="runge-kutta", fitter="direct", solution="BIC"):
        super().__init__(solver, fitter, solution)
        self.params['d'] = 0.50

    def __str__(self):
        return f"{self.params['c']: .2f} ⋅ V^{self.params['d']: .2f}"

    def __repr__(self):
        return f"Mendelsohn_growth(solver={self.solver}, fitter={self.fitter}, solution={self.solution})"

    def formula(self, y, params=None):
        """
        This is the Mendelsohn differential equation, used to calculate the
        growth. It is used by the ODE solvers, fit, and BIC or AICc calculations
        : param y: The y-value used to calculate the growth
        : param params: The parameters for the equation used to calculate the growth
        : return: The calculated growth (from the equation)
        """
        if y < 0:
            return 0
        if params is None:
            params = self.params
        c = params['c']
        d = params['d']
        return c * np.power(y, d)


class Montroll(BaseDifferential):
    """
    This class is used to run the functions of the `BaseDifferential` class with
    the Montroll differential equation.
    """
    def __init__(self, solver="runge-kutta", fitter="direct", solution="BIC"):
        super().__init__(solver, fitter, solution)
        self.params["V_max^d"] = 0.6
        self.params["d"] = 0.10

    def __str__(self):
        return f"{self.params['c']: .2f} ⋅ V ⋅ {self.params['V_max^d']: .2f} - V^{self.params['d']: .2f}"

    def __repr__(self):
        return f"Montroll(solver={self.solver}, fitter={self.fitter}, solution={self.solution})"

    def formula(self, y, params=None):
        """
        This is the Montroll differential equation, used to calculate the
        growth. It is used by the ODE solvers, fit, and BIC or AICc calculations
        : param y: The y-value used to calculate the growth
        : param params: The parameters for the equation used to calculate the growth
        : return: The calculated growth (from the equation)
        """
        V_max = params["V_max^d"]
        if y < 0:
            return 0
        if params is None:
            params = self.params
        c = params['c']
        d = params['d']
        return c * y * (V_max - np.power(y, d))


class Gompertz(BaseDifferential):
    """
    This class is used to run the functions of the `BaseDifferential` class with
    the Gompertz differential equation.
    """
    def __init__(self, solver="runge-kutta", fitter="direct", solution="BIC"):
        super().__init__(solver, fitter, solution)
        self.params["V_max"] = 0.90

    def __str__(self):
        return f"{self.params['c']: .2f} ⋅ V ⋅ ln({self.params['V_max']: .2f} / V)"

    def __repr__(self):
        return f"Gompertz(solver={self.solver}, fitter={self.fitter}, solution={self.solution})"

    def formula(self, y, params=None):
        """
        This is the Gompertz differential equation, used to calculate the
        growth. It is used by the ODE solvers, fit, and BIC or AICc calculations
        : param y: The y-value used to calculate the growth
        : param params: The parameters for the equation used to calculate the growth
        : return: The calculated growth (from the equation)
        """
        V_max = params["V_max"]
        if y <= 0 or V_max <= 0:
            return 0
        if params is None:
            params = self.params
        c = params['c']
        return c * y * np.log(abs(V_max / y))


class Exponential(BaseDifferential):
    """
    This class is used to run the functions of the `BaseDifferential` class with
    the Exponential differential equation.
    """
    def __str__(self):
        return f"{self.params['c']: .2f} ⋅ V"

    def __repr__(self):
        return f"Exponential(solver={self.solver}, fitter={self.fitter}, solution={self.solution})"

    def formula(self, y, params=None):
        """
        This is the Exponential differential equation, used to calculate the
        growth. It is used by the ODE solvers, fit, and BIC or AICc calculations
        : param y: The y-value used to calculate the growth
        : param params: The parameters for the equation used to calculate the growth
        : return: The calculated growth (from the equation)
        """
        if y < 0:
            return 0
        if params is None:
            params = self.params
        c = params['c']
        return c * y

class Logistic(BaseDifferential):
    """
    This class is used to run the functions of the `BaseDifferential` class with
    the Logistic differential equation.
    """
    def __init__(self, solver="runge-kutta", fitter="direct", solution="BIC"):
        super().__init__(solver, fitter, solution)
        self.params["V_max"] = 1.00

    def __str__(self):
        return f"{self.params['c']: .2f} ⋅ ({self.params['V_max']: .2f} - V)"

    def __repr__(self):
        return f"Logistic(solver={self.solver}, fitter={self.fitter}, solution={self.solution})"

    def formula(self, y, params=None):
        """
        This is the logistic differential equation, used to calculate the
        growth. It is used by the ODE solvers, fit, and BIC or AICc calculations
        : param y: The y-value used to calculate the growth
        : param params: The parameters for the equation used to calculate the growth
        : return: The calculated growth (from the equation)
        """
        if y < 0:
            return 0
        if params is None:
            params = self.params
        c = params['c']
        V_max = params["V_max"]
        return c * y * (V_max - y)


class LinearLimited(BaseDifferential):
    """
    This class is used to run the functions of the `BaseDifferential` class with
    the Linear limited differential equation.
    """
    def __init__(self, solver="runge-kutta", fitter="direct", solution="BIC"):
        super().__init__(solver, fitter, solution)
        self.params['d'] = 0.5

    def __str__(self):
        return f"{self.params['c']: .2f} ⋅ V / (V ⋅ {self.params['d']: .2f})"

    def __repr__(self):
        return f"LinearLimited(solver={self.solver}, fitter={self.fitter}, solution={self.solution})"

    def formula(self, y, params=None):
        """
        This is the linear limited differential equation, used to calculate the
        growth. It is used by the ODE solvers, fit, and BIC or AICc calculations
        : param y: The y-value used to calculate the growth
        : param params: The parameters for the equation used to calculate the growth
        : return: The calculated growth (from the equation)
        """
        if y < 0:
            return 0
        if params is None:
            params = self.params
        c = params['c']
        d = params['d']
        if y + d == 0:
            return 0
        return c * y / (y + np.exp(d))

class ExponentialAddition(BaseDifferential):
    """
    This class is used to run the functions of the `BaseDifferential` class with
    the Exponential addition differential equation.
    """
    def __init__(self, solver="runge-kutta", fitter="direct", solution="BIC"):
        super().__init__(solver, fitter, solution)
        self.params['d'] = 1.00

    def __str__(self):
        return f"{self.params['c']: .2f} ⋅ V + {self.params['d']: .2f}"

    def __repr__(self):
        return f"ExponentialAddition(solver={self.solver}, fitter={self.fitter}, solution={self.solution})"

    def formula(self, y, params=None):
        """
        This is the Exponential addition differential equation, used to calculate the
        growth. It is used by the ODE solvers, fit, and BIC or AICc calculations
        : param y: The y-value used to calculate the growth
        : param params: The parameters for the equation used to calculate the growth
        : return: The calculated growth (from the equation)
        """
        if y < 0:
            return 0
        if params is None:
            params = self.params
        c = params['c']
        d = params['d']
        return c * y  + d
