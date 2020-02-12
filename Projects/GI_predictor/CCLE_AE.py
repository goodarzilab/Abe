import re, sys, os
from time import time
from sklearn.preprocessing import LabelEncoder, normalize
from sklearn.model_selection import train_test_split
from keras import metrics
from keras import backend as K
from keras.layers import Input, Dense
from keras.models import Model
from keras.callbacks import History, EarlyStopping
from keras.utils import plot_model, model_to_dot
import numpy as np
import pandas as pd
import matplotlib. pyplot as plt
import argparse
import warnings
warnings.filterwarnings("ignore")


def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--user", help="Absolute path to the input directories", type=str)
    parser.add_argument("-e", "--epochs", help="number of epochs", type=int, default=5)
    parser.add_argument("-s", "--structure", help="Structure of AE seprated by ','\n order:\
                        1) Encoder number of dimenstions \
                        2) Number of layer after input before encoder (min = 1)\
                        3) Number of layer after encoder before output (min = 1)\
                        4) Multiplication of (1) for dimntions of layers before/after encoded layer\
                        ", type=str, default='32,1,1,2')
    parser.add_argument("-b", "--batch_size", help="Batch size", type=int, default=32)
    parser.add_argument("-t", "--test_size", help="Test size ratio", type=float, default=0.1)
    parser.add_argument("-D", "--dataset", help="Absolute path to CCLE dataset", type=str)
    parser.set_defaults(dataset='/rumi/shams/abe/Datasets/CCLE/',user='/rumi/shams/abe/')
    args = parser.parse_args()
    return args


def read_data(dtype='counts'):
    global m_rna_counts, m_rna_rpkm, cell_annotations
    m_rna_counts = data_folder + 'CCLE_RNAseq_genes_counts_20180929.gct.gz'
    m_rna_rpkm = data_folder + 'CCLE_RNAseq_genes_rpkm_20180929.gct.gz'
    cell_annotations = data_folder + 'Cell_lines_annotations_20181226.txt'
    # select input format
    if dtype == 'counts':
        m_rna = m_rna_counts
    elif dtype == 'rpkm':
        m_rna = m_rna_rpkm

    # read raw data
    raw_m_rna = pd.read_csv(m_rna, skiprows=2, sep='\t')

    # make meta data dictionary
    meta = {'m_rna': raw_m_rna[['Name','Description']],
            'cell_lines':raw_m_rna.columns.values.tolist()[2:],
            'annotations':pd.read_csv(cell_annotations, sep='\t')
           }
    # normalize
    df_m_rna = raw_m_rna.drop(columns=['Name','Description']).to_numpy()
    df_m_rna = normalize(X=df_m_rna, axis=0, norm="max")
    data = {'df':df_m_rna, 'meta': meta}
    return data


def r_square(y_true, y_pred):
    # custom R2-score metrics for keras backend
    SS_res =  K.sum(K.square(y_true - y_pred))
    SS_tot = K.sum(K.square(y_true - K.mean(y_true)))
    return ( 1 - SS_res/(SS_tot + K.epsilon()) )


def create_model(input_size):
    # this is our input placeholder
    myinput = Input(shape=(input_size,))
    ### hidden layers
    # encoding
    n = n_en
    encoded = Dense(dim*(mu**n_en), activation='relu')(myinput)
    while n > 1:
        n-=1
        encoded = Dense(dim*(mu**n), activation='relu')(encoded)
    encoded = Dense(dim, activation='relu')(encoded)
    decoded = Dense(dim*mu, activation='relu')(encoded)
    # decoding
    while n < n_de:
        n+=1
        decoded = Dense(dim*(mu**n), activation='relu')(decoded)
    decoded = Dense(input_size, activation='linear')(decoded)
    ### Separate encoder & autoencoder models
    autoencoder = Model(myinput, decoded)
    encoder = Model(myinput, encoded)
    autoencoder.compile(
        optimizer='adadelta',
        loss='mse',
        metrics=[r_square]
    )
    return autoencoder, encoder


def plot_by_epochs(model,epochs):
    epochs = range(epochs)
    # plot loss
    loss = model.history.history['loss']
    val_loss = model.history.history['val_loss']
    plt.figure()
    plt.plot(epochs, loss, 'bo', label='Training loss')
    plt.plot(epochs, val_loss, 'b', label='Validation loss')
    plt.title('Training and validation loss')
    plt.legend()
    plt.savefig(model_path+'loss.pdf')
    plt.savefig(model_path+'loss.png')
    # plot R^2
    acc = model.history.history['r_square']
    val_acc = model.history.history['val_r_square']
    plt.figure()
    plt.plot(epochs, acc, 'bo', label='Training r_square')
    plt.plot(epochs, val_acc, 'b', label='Validation r_square')
    plt.title('Coefficient of determination')
    plt.ylabel('R^2')
    plt.xlabel('Epoch')
    plt.legend(['Trian', 'Test'])
    plt.savefig(model_path+'r_square.pdf')
    plt.savefig(model_path+'r_square.png')


def run_autoencoder(df):
    # make X and y from df
    X = df.T
    y_true = X
    X_train, X_test,_,_ = train_test_split(X, X, test_size=test_size, random_state=42)
    # callbacks
    print (f'Input shape: {X.shape[0]} Cell lines X {X.shape[1]} RNAseq genes')
    print (f'Batch size {batch_size}, epochs {epochs}')
    history = History()
    early_stopping = EarlyStopping(monitor='val_loss', patience=20)
    # make, save plot and fit the model
    autoencoder, encoder = create_model(X.shape[1])
    print(autoencoder.summary())
    autoencoder.fit(X_train, X_train,
              epochs=epochs,
              batch_size=batch_size,
              shuffle=False, # changed to false to keep val cell equal for comparing models
              callbacks=[history,early_stopping],
              validation_data=(X_test, X_test)
             )
    print("fitting has just been finished")
    X_pred = autoencoder.predict(X_test, batch_size=batch_size, verbose=2)
    print("prediction process has just been finished")
    # save the model and encoded-layer output
    autoencoder.save(filepath=model_path+"autoencoder.h5")
    encoder.save(filepath=model_path+"encoder.h5")
    # save the result and prediction value
    np.savetxt(X=X_test, fname=model_path+"X_test.csv", delimiter=",")
    np.savetxt(X=X_pred[0], fname=model_path+"X_pred.csv", delimiter=",")
    if early_stopping.stopped_epoch == 0:
        plot_by_epochs(autoencoder,epochs)
    elif early_stopping.stopped_epoch > 0:
        print (f'model stopped training at epoch {early_stopping.stopped_epoch}')
        # plot_by_epochs(autoencoder,early_stopping.stopped_epoch)
    print("model objects and metrics plots saved")


def main():
    global args, user, data_folder
    args = handler()
    user = args.user
    data_folder = args.dataset
    # Parse model parameters
    global test_size, epochs, batch_size, structure, dim,n_en,n_de,mu
    test_size = args.test_size
    batch_size = args.batch_size
    epochs = args.epochs
    dim,n_en,n_de,mu = list(map(int,args.structure.split(',') ))
    #
    global results_folder, model_structure, model_path
    results_folder = user + 'Projects/GI_predictor/Models/'
    model_structure = 'encoding-'+str(dim)+'-dim-mu-'+str(mu)+'--'+str(n_en)+'-en-'+str(n_de)+'-de-layers'
    model_path = results_folder+model_structure +"_batch_size_"+str(batch_size)+"_epochs-"+str(epochs)+'/'
    if not os.path.exists(model_path):
        os.mkdir(model_path)
    t0 = time()
    # read CCLE data
    print('\n************************* Read CCLE Data *************************\n')
    data = read_data()
    # cell_lines = data['meta']['cell_lines']
    # annotations = data['meta']['annotations']
    df = data['df']
    print('\n************************* run autoencoder ************************\n')
    run_autoencoder(df)
    print("done in %0.3fs" % (time() - t0))


if __name__ == '__main__':
    main()
