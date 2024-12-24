# differential


## About

When presented with data of which you need the know the mechanisms it is possible to try and fit a differential model to it. Mostly for the purpose of research. To fit the data the package will use the following standard models:
* Exponential Growth
* Exponential Addition Growth
* Gompertz Growth
* Logistic Growth
* Linear Growth
* Mendelsohn Growth
* Montreoll Growth

Other differential models can be easily added as will be discussed int the `Expandable` chapter. The package will makes it possible for the user to fit each model using one of our three built in solvers. These solvers are the following:

* Euler
* Huen
* Runge-Kutta

A [Mean Sqeured Error (MSE)](https://en.wikipedia.org/wiki/Mean_squared_error) will be calculated by the solver to determine the best fit. The lower this MSE, the better the fit. This package supports [AICc](https://en.wikipedia.org/wiki/Akaike_information_criterion) and [BIC](https://en.wikipedia.org/wiki/Bayesian_information_criterion) as a means of quality assessment, when comparing different models make sure the same quality method.

# Models

The following are the current supported models. Each can be used to fit the data presented. 

### Parameters

Not all parameters are present in each model. Below is a quick rundown of what they represent:

|      Parameter      |   Symbol   | Definition                                                                  |
|:--------------:|:----------:|-----------------------------------------------------------------------------|
|      Time      |   $ t $    | The independent variable, representing time (x-value).                                |
|    Constant    |   $ c $    | Growth rate constant influencing the dynamics.                              |
|     Volume     |   $ V $    | Represents the changing volume (y-value).                                             |
|     Delta      |   $ d $    | An additive component in the growth formula.                                |
| Maximum Volume | $ V\_max $ | The upper limit for the volume in models like Gompertz and Logistic Growth. |



## Exponential Growth

Exponential Growth shows a pattern of data that increases over time. Fastly adding volume, by each step of time passing. It resembles an Exponential function.

### Formula

$$
\frac{\text{d}V}{\text{d}t} = c \cdot V
$$

## Exponential Addition Growth

It is alike to Exponential growth. Exponential Addition Growth shows a patten in data that increases over time. Unlike Exponential Growth it has an additional component.
### Formula

$$
\frac{\text{d}V}{\text{d}t} = c \cdot V + d
$$

## Gompertz Growth

Named after Benjamin Gompertz, this type of growth curve is classified as a sigmoid function. Meaning growth is slowest at the start and the end of a time series.

### Formula
$$\frac{\text{d}V}{\text{d}t} = c \cdot V \cdot \ln \left( \frac{V\_max}{V} \right)$$


## Logistic Growth

A sigmoid curve type of growth. Meaning growth is slowest at the start and the end of a time series.  Similar to Gompertz Growth. 

### Formula
$$
\frac{\text{d}V}{\text{d}t} = c \cdot V \cdot \left( V\_max - V \right)
$$

## Linear Growth

This represents the Growth at its basis. A steady increase of volume over a given time series.
### Formula 
$$
\frac{\text{d}V}{\text{d}t} = c
$$

## Mendelsohn Growth

Built upon the exponential growth it can make an adjustment to the exponent. This wat proposed my Mendolsohn in 1963. 

### Formula
$$
\frac{\text{d}V}{\text{d}t} = c \cdot V^d
$$


## Montroll Growth

Just like logistic the growth decreases over time as the volume approaches the maximum value. The difference lies in the additive exponential component (`d`).

### Formula

The variable $V\_max^d$ is in our model will be considered as one. 
$$
\frac{\text{d}V}{\text{d}t} = c \cdot V \cdot \left( V\_max^d - V^d \right)
$$

# Usage

Choose a model from the `differential_eq.py` to fit. Each model takes the following parameters for initialization:

|   name   | options                        | default       |             definition             |
|:--------:|--------------------------------|---------------|:----------------------------------:|
|  solver  | "runge-kutta", "heun", "euler" | "runge-kutta" | Type of differential solver used.  |
|  fitter  | "direct", "random"             | "direct"      |    Parameter finding strategy.     |
| solution | "AICC", "BIC"                  | "BIC"         | Determines the quality of the fit. |


The data the model requires will be a time series `t` and data points `y`. These will reside in a list with the equal length. Use these as parameters to call the `model.fit(ts, ys)` method. Where `ts` and `ys` are time and volume data respectively.

For results and output consult the `Example` chapter. For the explination of the implementation, examine the well documented code in `differential_eq.py`.

## Scaling data
If the data includes high numbers the program could run into out of memory numbers and give unexpected results. It is recommended to scale the dataset `y` between 0 and 1 and the `t` dataset between -1.0 and 1.0 respectively. 

# Example

For use cases look at the jupyter-notebooks `performing.ipynb` and `results.ipynb`. 

# Expendable

It is possible to add personalized for own usage. This can be achieved by creating child objects from the `BaseDifferential` Class.

## Logistic Example
With the example below is shown how to implement this.
```
class Logistic(BaseDifferential):
    """
    This class is used to run the functions of the `BaseDifferential` class with
    the Logistic differential equation.
    """
    def __init__(self, solver, fitter, solution):
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
```

Create a `__init__` which and add the parameters needed to the `self.params` dictionary. And call the `super().__init__(solver, fitter, solution)` for initialization.
Create a `formula` function copy the `if` statements as shown. Below implement the personal model to your satisfaction.

It is optional to create the `__str__` and `__repr__` methods.
