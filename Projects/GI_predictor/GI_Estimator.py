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
from main_functions import get_SLdataset, get_gene_pairs, r_square, plot_by_epochs

warnings.filterwarnings("ignore")


def handler():
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--user", help="Absolute path to the input directories", type=str)
    parser.add_argument("-e", "--epochs", help="number of epochs", type=int, default=5)
    parser.add_argument("-s", "--structure", help="Structure ...", type=str, default='')
    parser.add_argument("-c", "--cellline", help="train model on screening the cell line ", type=str, default='K562')
    parser.add_argument("-b", "--batch_size", help="Batch size", type=int, default=32)
    parser.add_argument("-t", "--test_size", help="Test size ratio", type=float, default=0.1)
    parser.set_defaults(dataset='/rumi/shams/abe/Datasets/CCLE/',user='/rumi/shams/abe/')
    args = parser.parse_args()
    return args


def read_data(cell_line):
    SL_dataset = get_SLdataset()
    f_path = 'cBioPortal/aml_ohsu_2018/' + 'data_RNA_Seq_expression_cpm_Zscores.txt'
    data = pd.read_csv(f_path,sep='\t', index_col='Hugo_Symbol', na_values ='NA').drop(columns='Entrez_Gene_Id').astype(float)   
    data_G1, data_G2 = get_gene_pairs(SL_dataset[cell_line], data)
    X = np.concatenate((    np.array (data_G1.T),    np.array (data_G2.T))).T
    # this is our input:
    print(X.shape)
    print (f'rows: {X.shape[0]} (#gene pairs)\ncolumns: {X.shape[1]} (2*#patients)')
    return X


def create_model(input_size, dim):
    # this is our input placeholder
    myinput = Input(shape=(input_size,))
    ### hidden layers
    deep = Dense(dim, activation='relu')(myinput)
    ### add cell-line features from encoder model
    deep #??
    # last layer 
    last = Dense(1, activation='sigmoid')(deep)
    model = Model(myinput, last)   
    model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy', r_square])
    return model


def run_model(X,y):
    # make X and y from df
    y_true = y
    X_train, X_test,y_train,y_test = train_test_split(X, y, test_size=test_size, random_state=42)
    # callbacks
    history = History()
    early_stopping = EarlyStopping(monitor='val_loss', patience=20)
    # make, save plot and fit the model
    model = create_model(X.shape[1])
    print(model.summary())
    model.fit(X_train, y_train,
              epochs=epochs,
              batch_size=batch_size,
              shuffle=False, # changed to false to keep val cell equal for comparing models
              callbacks=[history,early_stopping],
              validation_data=(X_test, y_test)
             )
    print("fitting has just been finished")
    X_pred = model.predict(X_test, batch_size=batch_size, verbose=2)
    print("prediction process has just been finished")
    # save the model and encoded-layer output
    model.save(filepath=model_path+"autoencoder.h5")
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
    global args, user
    args = handler()
    user = args.user
    # Parse model parameters
    global test_size, epochs, batch_size, structure
    test_size = args.test_size
    batch_size = args.batch_size
    epochs = args.epochs
    #
    global results_folder, model_structure, model_path
    results_folder = user + 'Projects/GI_predictor/Models/'
    model_structure = str(mu)+'--'+str(n_en)+'-en-'+str(n_de)+'-de-layers'
    model_path = results_folder+model_structure +"_batch_size_"+str(batch_size)+"_epochs-"+str(epochs)+'/'
    if not os.path.exists(model_path):
        os.mkdir(model_path)
    t0 = time()
    # read cohort data
    print('\n************************* Read trainable cohort data ***************\n')
    X = read_data()
    print('\n************************* run model ********************************\n')
    run_model(X,y)
    print("done in %0.3fs" % (time() - t0))


if __name__ == '__main__':
    main()
