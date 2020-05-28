import torch
import torch.nn as nn


class Model():

    def __init__(params)

        # Data params
        self.n_features = params['n_features']
        self.seq_length = params['seq_length']

        # RNN params
        self.hidden_dim = params['hidden_dim']
        self.n_layers = params['n_layers']

        self.model = RNN(   input_size=self.n_features, 
                            output_size=self.n_features, 
                            hidden_dim=self.hidden_dim,     
                            n_layers=self.n_layers)

        self.device = self.get_device()

    def train_model(x_data,y_data):


def build_sequences(x,y,seq_length):

    n,n_features = x.shape 
    n = len(x)

    x_data = np.zeros([n-seq_length,seq_length,n_features])
    y_data = np.zeros([n-seq_length,1,n_features])

    for i in range(n-seq_length):
        x_data[i,::] = x[i:i+seq_length,:]
        y_data[i,::] = y[i+seq_length-1,:]

    return x_data, y_data

def shuffle(x_data,y_data):
    n, _, _ = x_data.shape
    idx = np.arange(n)
    np.random.shuffle(idx)

    return x_data[idx,::], y_data[idx,::]

def get_device():
    is_cuda = torch.cuda.is_available()
    if is_cuda:
        device = torch.device("cuda")
        print("GPU is available")
    else:
        device = torch.device("cpu")
        print("GPU not available, CPU used")

    return device 

def train_model(n_features,seq_length,train_x,train_y,val_x,val_y):
    '''
        
    '''

    device = get_device()

    train_x_seq, train_y_seq = build_sequences(train_x,train_y,seq_length)
    val_x_seq, val_y_seq = build_sequences(val_x,val_y,seq_length)

    # Instantiate the model with hyperparameters
    model = Model(input_size=n_features, output_size=n_features, hidden_dim=100, n_layers=1)
    model.to(device)

    # Define hyperparameters
    n_epochs = 20000
    lr=0.01

    # Define Loss, Optimizer
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    # Send data to torch and device
    train_x_tensor = torch.Tensor(train_x_seq)
    train_y_tensor = torch.Tensor(train_y_seq)
    # train_y_tensor = train_y_tensor.permute(1, 0, 2) # To match network output dim

    train_x_tensor.to(device)
    train_y_tensor.to(device)

    # val_x_seq.to(device)
    # val_y_seq.to(device)

    # Training Run
    for epoch in range(1, n_epochs + 1):
        
        train_x_tensor, train_y_tensor = shuffle(train_x_tensor, train_y_tensor)

        # Zero opt
        optimizer.zero_grad() # Clears existing gradients from previous epoch
        
        # Forward passs
        output = model(train_x_tensor)
        loss = criterion(output, train_y_tensor.permute(1, 0, 2))
        loss.backward() 
        optimizer.step()
        
        if epoch%10 == 0:
            print('Epoch: {}/{}.............'.format(epoch, n_epochs), end=' ')
            print("Loss: {:.4f}".format(loss.item()))


def prediction()



