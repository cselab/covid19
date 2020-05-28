import math
import torch
import gpytorch
from matplotlib import pyplot as plt

# We will use the simplest form of GP model, exact inference
class ExactGPModel(gpytorch.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood):
        super(ExactGPModel, self).__init__(train_x, train_y, likelihood)
        self.mean_module = gpytorch.means.ConstantMean()
        self.covar_module = gpytorch.kernels.ScaleKernel(gpytorch.kernels.RBFKernel())
        self.linear1 = torch.nn.Linear(1, 100, bias=True)
        self.linear2 = torch.nn.Linear(100, 1, bias=True)
    def forward(self, x):
        # # print(x.size())
        # # print(x.min())
        # # print(x.max())
        # # print(ark)
        # x = self.linear1(x)
        # # x = self.tanh(x)
        # x = self.linear2(x)
        # # print(x.size())
        mean_x = self.mean_module(x)
        covar_x = self.covar_module(x)
        return gpytorch.distributions.MultivariateNormal(mean_x, covar_x)

def train_model(x_data,y_data):


    # Training data
    train_x = torch.Tensor(x_data)
    train_y = torch.Tensor(y_data)

    # initialize likelihood and model
    likelihood = gpytorch.likelihoods.GaussianLikelihood()
    model = ExactGPModel(train_x, train_y, likelihood)


    training_iter = 20000
    # training_iter = 1000

    # Find optimal model hyperparameters
    model.train()
    likelihood.train()

    # Use the adam optimizer
    optimizer = torch.optim.Adam([
        {'params': model.parameters()},  # Includes GaussianLikelihood parameters
    ], lr=0.1)

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

    model.eval()
    likelihood.eval()

    # Test points are regularly spaced along [0,1]
    # Make predictions by feeding model through likelihood
    with torch.no_grad(), gpytorch.settings.fast_pred_var():
        observed_pred = likelihood(model(train_x))

        # Initialize plot
        f, ax = plt.subplots()
        # Get upper and lower confidence bounds
        lower, upper = observed_pred.confidence_region()
        # Plot training data as black stars
        ax.plot(train_x.numpy(), train_y.numpy(), 'go', label='Observed Data')
        # Plot predictive means as blue line
        ax.plot(train_x.numpy(), observed_pred.mean.numpy(), 'rx', label='Predictions')
        # Shade between the lower and upper confidence bounds
        ax.fill_between(train_x.numpy(), lower.numpy(), upper.numpy(), color='r', alpha=0.5, label='Confidence')
        ax.set_xlabel("Total cases")
        ax.set_ylabel("New cases")
        lgd = ax.legend(loc="upper left", bbox_to_anchor=(1.05, 1), borderaxespad=0., frameon=False)
        plt.tight_layout()
        # plt.show()
        plt.savefig("./plots/prediction_NC_TC.pdf")

        train_y_np = train_y.numpy()
        train_time_series_np = np.cumsum(train_y_np)

        test_y_pred_np = observed_pred.mean.numpy()
        test_time_series_pred_np = np.cumsum(test_y_pred_np)

        test_y_pred_lower_np = lower.numpy()
        test_y_pred_upper_np = upper.numpy()
        # test_y_pred_lower_np = test_y_pred_np[0] + test_y_pred_lower_np
        # test_y_pred_upper_np = test_y_pred_np[0] + test_y_pred_upper_np
        test_y_pred_lower_np = np.cumsum(test_y_pred_lower_np)
        test_y_pred_upper_np = np.cumsum(test_y_pred_upper_np)

        f, ax = plt.subplots()
        # Plot training data as black stars
        ax.plot(np.arange(np.shape(train_time_series_np)[0]), train_time_series_np, 'g--', label="Train")
        ax.plot(np.arange(np.shape(test_time_series_pred_np)[0]), test_time_series_pred_np, 'r--', label="Prediction")
        ax.fill_between(np.arange(np.shape(test_time_series_pred_np)[0]), test_y_pred_lower_np, test_y_pred_upper_np, color='r', alpha=0.5, label='Confidence')
        ax.set_xlabel("Day number")
        ax.set_ylabel("Total cases")
        ax.legend()
        # plt.show()
        plt.savefig("./plots/Prediction_TC_DN.pdf")