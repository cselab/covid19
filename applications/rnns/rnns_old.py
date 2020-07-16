#!/usr/bin/env python3
# Author: Pantelis R. Vlachas
# Date:   20/4/2020
# Email:  pvlachas@ethz.ch
from pathlib import Path
import numpy as np

import json
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'


import pickle
import sys

PROJECT_DIR =  Path("/Users/pvlachas/Documents/PhD/project_covid/covid19")


sys.path.append(str(PROJECT_DIR))

from epidemics.utils.io import download_and_save


DATA_DIR = Path("/Users/pvlachas/Documents/PhD/project_covid/covid19/rnns/data")
DATA_DOWNLOADS_DIR = DATA_DIR / 'downloads'

print(DATA_DIR)
print(DATA_DOWNLOADS_DIR)


##########################################################################
## GPYTORCH - Gaussian processes in pytorch
##########################################################################
import math
import torch
import gpytorch
from matplotlib import pyplot as plt


os.makedirs(str(DATA_DIR), exist_ok=True)

def fetch_openzh_covid_data(*, cache_duration=3600):
    """
    Returns a dictionary of lists {canton abbreviation: number of cases per day}.
    """
    url = 'https://raw.githubusercontent.com/daenuprobst/covid19-cases-switzerland/master/covid19_cases_switzerland_openzh.csv'
    # url = 'https://raw.githubusercontent.com/daenuprobst/covid19-fatalities-switzerland/master/covid19_cases_switzerland_openzh.csv'
    path = DATA_DOWNLOADS_DIR / 'covid19_cases_switzerland_openzh.csv'

    raw = download_and_save(url, path, cache_duration=cache_duration)
    rows = raw.decode('utf8').split()
    cantons = rows[0].split(',')[1:-1]  # Skip the "Date" and "CH" cell.

    data = {canton: [] for canton in cantons}
    for day in rows[1:]:  # Skip the header.
        cells = day.split(',')[1:-1]  # Skip "Date" and "CH".
        assert len(cells) == len(cantons), (len(cells), len(cantons))

        for canton, cell in zip(cantons, cells):
            data[canton].append(float(cell or 'nan'))
    return data


cases_per_country = fetch_openzh_covid_data()

# print(cases_per_country)

# for key in cases_per_country:
#     print(key)
#     print(np.shape(cases_per_country[key]))
#     print(cases_per_country[key])

### TODO: Correct Nans

# time_series = [cases_per_country["ZH"], cases_per_country["VS"]]
time_series_train = cases_per_country["ZH"]
time_series_train[0] = 0.0
time_series_train[1] = 0.0
time_series_train = np.array(time_series_train)
print(time_series_train)

time_series_test = cases_per_country["VS"]
time_series_test[0] = 0.0
time_series_test[1] = 0.0
time_series_test[2] = 0.0
time_series_test[11] = 0.5 * (time_series_test[10] + time_series_test[12])
time_series_test = np.array(time_series_test)
print(time_series_test)


test_input  = time_series_test[:-1]
test_target = time_series_test[1:] - time_series_test[:-1]

# Building the training and testing data-sets
test_x = torch.Tensor(test_input)
print(test_x.size())

test_y = torch.Tensor(test_target)
print(test_y.size())


train_input  = time_series_train[:-1]
train_target = time_series_train[1:] - time_series_train[:-1]

# Building the training and testing data-sets
train_x = torch.Tensor(train_input)
print(train_x.size())

train_y = torch.Tensor(train_target)
print(train_y.size())

# # Initialize plot
# f, ax = plt.subplots()
# # Plot training data as black stars
# ax.plot(train_x.numpy(), train_y.numpy(), 'go', label="Train")
# ax.plot(test_x.numpy(), test_y.numpy(), 'bx', label="Test")
# ax.legend()
# plt.show()

# train_x_np = train_x.numpy()
# test_x_np = test_x.numpy()
# # Initialize plot
# f, ax = plt.subplots()
# # Plot training data as black stars
# ax.plot(np.arange(np.shape(train_x_np)[0]), train_x_np, 'g--', label="Train")
# ax.plot(np.arange(np.shape(test_x_np)[0]), test_x_np, 'b--', label="Test")
# ax.legend()
# plt.show()





# We will use the simplest form of GP model, exact inference
class ExactGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood):
        super(ExactGPModel, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel())

    def forward(self, x):
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)

# initialize likelihood and model
likelihood = gpytorch.likelihoods.GaussianLikelihood()
model = ExactGPModel(train_x, train_y, likelihood)


training_iter = 100000


# Find optimal model hyperparameters
model.train()
likelihood.train()

# Use the adam optimizer
optimizer = torch.optim.Adam([
    {'params': model.parameters()},  # Includes GaussianLikelihood parameters
], lr=0.01)

# "Loss" for GPs - the marginal log likelihood
mll = gpytorch.mlls.ExactMarginalLogLikelihood(likelihood, model)

for i in range(training_iter):
    # Zero gradients from previous iteration
    optimizer.zero_grad()
    # Output from model
    output = model(train_x)
    # Calc loss and backprop gradients
    loss = -mll(output, train_y)
    loss.backward()
    print('Iter %d/%d - Loss: %.3f   lengthscale: %.3f   noise: %.3f' % (
        i + 1, training_iter, loss.item(),
        model.covar_module.base_kernel.lengthscale.item(),
        model.likelihood.noise.item()
    ))
    optimizer.step()


# f_preds = model(test_x)
# y_preds = likelihood(model(test_x))

# f_mean = f_preds.mean
# f_var = f_preds.variance
# f_covar = f_preds.covariance_matrix
# f_samples = f_preds.sample(sample_shape=torch.Size(1000,))

# Get into evaluation (predictive posterior) mode
model.eval()
likelihood.eval()

# Test points are regularly spaced along [0,1]
# Make predictions by feeding model through likelihood
with torch.no_grad(), gpytorch.settings.fast_pred_var():
    observed_pred = likelihood(model(test_x))

    # Initialize plot
    f, ax = plt.subplots()
    # Get upper and lower confidence bounds
    lower, upper = observed_pred.confidence_region()
    # Plot training data as black stars
    ax.plot(train_x.numpy(), train_y.numpy(), 'go', label='Observed Data')
    # Plot predictive means as blue line
    ax.plot(test_x.numpy(), observed_pred.mean.numpy(), 'rx', label='Predictions')
    # Shade between the lower and upper confidence bounds
    ax.fill_between(test_x.numpy(), lower.numpy(), upper.numpy(), color='r', alpha=0.5, label='Confidence')
    # Plot GROUNDTRUTH test data as red line
    ax.plot(test_x.numpy(), test_y.numpy(), 'bx', label='Groundthruth')
    lgd = ax.legend(loc="upper left", bbox_to_anchor=(1.05, 1), borderaxespad=0., frameon=False)
    plt.tight_layout()
    plt.show()

    train_y_np = train_y.numpy()
    train_time_series_np = np.cumsum(train_y_np)
    test_y_np = test_y.numpy()
    test_time_series_np = np.cumsum(test_y_np)
    test_y_pred_np = observed_pred.mean.numpy()
    test_time_series_pred_np = np.cumsum(test_y_pred_np)
    f, ax = plt.subplots()
    # Plot training data as black stars
    ax.plot(np.arange(np.shape(train_time_series_np)[0]), train_time_series_np, 'g--', label="Train")
    ax.plot(np.arange(np.shape(test_time_series_np)[0]), test_time_series_np, 'b--', label="Groundtruth")
    ax.plot(np.arange(np.shape(test_time_series_pred_np)[0]), test_time_series_pred_np, 'r--', label="Prediction")
    ax.legend()
    plt.show()


# train_input  = time_series_train[:-1]
# train_target = time_series_train[1:] - time_series_train[:-1]

# # Building the training and testing data-sets
# train_x = torch.Tensor(train_input)
# print(train_x.size())

# train_y = torch.Tensor(train_target)
# print(train_y.size())
# 
# 
# 







