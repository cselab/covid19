import torch
import torch.nn as nn




class RNN(nn.Module):
    def __init__(self, input_size, output_size, hidden_dim, n_layers):
        super().__init__()
        # Defining some parameters
        self.hidden_dim = hidden_dim
        self.n_layers = n_layers

        # RNN Layer
        self.rnn = nn.RNN(input_size, hidden_dim, n_layers, batch_first=True)   
        # Fully connected layer
        self.fc = nn.Linear(hidden_dim, output_size)
    
    def forward(self, x):
        
        _, h_T = self.rnn(x, self.hidden_0)
        out = self.fc(h_T)
        
        return out
    
    def init_hidden(self, batch_size):
        # This method generates the first hidden state of zeros which we'll use in the forward pass
        # We'll send the tensor holding the hidden state to the device we specified earlier as well
        self.hidden_0 = torch.zeros(self.n_layers, batch_size, self.hidden_dim).cuda()
        print(self.hidden_0.device)

        return 

    def send_to_gpu(self):

        self.rnn.cuda()
        self.fc.cuda()
        return 0    