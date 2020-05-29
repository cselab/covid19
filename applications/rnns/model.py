import torch
import torch.nn as nn
import numpy as np
from pytorch_models import RNN
import matplotlib.pyplot as plt
import os 

class Model():

    def __init__(self,params):

        self.country = params['country']

        # Data params
        self.n_features = params['n_features']
        self.seq_length = params['seq_length']

        # RNN params
        self.hidden_dim = params['hidden_dim']
        self.n_layers = params['n_layers']

        # Define hyperparameters
        self.n_epochs = params['n_epochs']
        self.lr = params['lr']

        self.leading_zeros = int(np.ceil(np.log10(self.n_epochs)))
        self.device = self.get_device()
        self.build_model()
        
        self.name = self.country + '/'
        self.name += 'seq'+str(self.seq_length)
        self.name += '_hdim'+str(self.hidden_dim)+'_l'+str(self.n_layers)
        self.name += '_e'+str(self.n_epochs)+'_lr'+str(self.lr) 
        self.path = './results/'+self.name + '/'

        self.create_folder(self.path)
        self.create_folder(self.path+'training_predictions/')

        self.training_loss = []

    def build_model(self):
        self.model = RNN(   input_size=self.n_features, 
                            output_size=self.n_features, 
                            hidden_dim=self.hidden_dim,     
                            n_layers=self.n_layers)


        # Define Loss, Optimizer
        self.criterion = nn.MSELoss()
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=self.lr)

# ---------------------------------------------------------------------------- #
# -------------------------- Training + Prediction --------------------------- #
# ---------------------------------------------------------------------------- #

    def train_model(self,x_data,y_data):
        ''' 
            input shape [batch,seq_length,n_features]
        '''

        x_seq, y_seq = self.build_sequences(x_data,y_data,self.seq_length)
        x_seq_tensor = torch.Tensor(x_seq)
        y_seq_tensor = torch.Tensor(y_seq)

        # Training Run
        if torch.cuda.is_available():
            self.model.send_to_gpu()

        self.model.init_hidden(x_seq_tensor.shape[0])

        x_seq_tensor = x_seq_tensor.to(self.device)
        y_seq_tensor = y_seq_tensor.to(self.device)

        for epoch in range(1, self.n_epochs + 1):
            
            x_seq_tensor, y_seq_tensor = self.shuffle(x_seq_tensor, y_seq_tensor)

            # Zero opt
            self.optimizer.zero_grad() # Clears existing gradients from previous epoch
            
            # Forward passs
            output = self.model(x_seq_tensor)
            loss = self.criterion(output, y_seq_tensor.permute(1, 0, 2))
            loss.backward() 
            self.optimizer.step()
            
            if epoch%1000 == 0:
                self.training_loss.append(loss.item()*1e9)
                print('Epoch: {}/{}.............'.format(epoch, self.n_epochs), end=' ')
                print("Loss: {:.8f}".format(loss.item()*1e9))
                path = self.path+'training_predictions/'+str(epoch).zfill(self.leading_zeros)+'.pdf'
                predictions = output.detach().cpu().numpy()[0]
                self.test_model(x_data,y_data,path)
                self.print_loss(self.training_loss)

        predictions = output.detach().cpu().numpy()[0]
        self.test_model(x_data,y_data,self.path+'one_step_pred.pdf')

    def test_model(self,x_data,y_data,path):
        x_seq, y_seq = self.build_sequences(x_data,y_data,self.seq_length)
        x_seq_tensor = torch.Tensor(x_seq).to(self.device)
        predictions = self.model(x_seq_tensor).detach().cpu().numpy()[0]

        fig = plt.figure()
        plt.plot(self.seq_length+np.arange(np.shape(predictions)[0]),predictions,'o',label='Predictions')
        plt.plot(np.arange(np.shape(y_data)[0]),y_data,'o',label='Data')
        plt.legend()
        plt.savefig(path)
        plt.close()

    def iterative_prediction(self,x_data,start_day,n_predictions):
        self.model.init_hidden(1)

        input_x = np.expand_dims(x_data[start_day-self.seq_length:start_day,:],axis=0)

        input_x = torch.Tensor(input_x)
        input_x = input_x.to(self.device)
        output = self.model(input_x)

        predictions = np.zeros([n_predictions,self.n_features])
        for i in range(n_predictions):
            input_x = input_x.roll(-1,1)
            input_x[0,-1,:] = input_x[0,-2,:] + output
            predictions[i,:] = input_x[0,-1,:].detach().cpu().numpy()
            output = self.model(input_x)

        fig = plt.figure()
        plt.plot(start_day+np.arange(np.shape(predictions)[0]),predictions,'o',label='Predictions')
        plt.plot(np.arange(np.shape(x_data)[0]),x_data,'o',label='Data')
        plt.legend()
        plt.savefig(self.path+'iterative_pred.pdf')
        plt.close()

        return predictions




# ---------------------------------------------------------------------------- #
# ----------------------------------- Utils----------------------------------- #
# ---------------------------------------------------------------------------- #
    
    def print_loss(self,training_loss):

        fig = plt.figure()
        plt.semilogy(np.arange(len(training_loss)),training_loss)
        plt.savefig(self.path+'loss.pdf')
        plt.close()

    def build_sequences(self,x,y,seq_length):

        n,n_features = x.shape 
        n = len(x)

        x_data = np.zeros([n-seq_length,seq_length,n_features])
        y_data = np.zeros([n-seq_length,1,n_features])

        for i in range(n-seq_length):
            x_data[i,::] = x[i:i+seq_length,:]
            y_data[i,::] = y[i+seq_length-1,:]

        return x_data, y_data

    def shuffle(self,x_data,y_data):
        n, _, _ = x_data.shape
        idx = np.arange(n)
        np.random.shuffle(idx)

        return x_data[idx,::], y_data[idx,::]

    def get_device(self,):
        is_cuda = torch.cuda.is_available()
        if is_cuda:
            device = torch.device("cuda")
            print("GPU is available")
        else:
            device = torch.device("cpu")
            print("GPU not available, CPU used")

        return device 


    def create_folder(self,name):
        if not os.path.exists(name):
            os.makedirs(name)


    def plot_data(self,confirmed):

        confirmed = np.array(confirmed)
        total  = confirmed[:-1]
        daily = confirmed[1:] - confirmed[:-1]

        f, ax = plt.subplots(2)

        ax[0].plot(np.arange(np.shape(daily)[0]),daily,'o')
        ax[1].plot(np.arange(np.shape(total)[0]),total,'--')
        ax[0].set_xlabel("Day number")
        ax[0].set_ylabel("Daily Incidences")
        ax[1].set_xlabel("Day number")
        ax[1].set_ylabel("Total cases")

        plt.savefig(self.path+'original_data.pdf')

