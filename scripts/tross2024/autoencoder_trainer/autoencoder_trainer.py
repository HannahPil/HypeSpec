#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 16:12:10 2022

@author: mtross
"""
from tensorflow import keras  
from keras.callbacks import ModelCheckpoint
import pandas as pd
import tensorflow as tf
from sklearn.utils import shuffle

class AutoencoderTrainer:
    """
    A class for training an autoencoder model using TensorFlow and Keras.

    Attributes:
    dataset_path (str): Path to the dataset file.
    codings_size (int): Size of the encoding layer in the autoencoder.
    epochs (int): Number of epochs for training.
    patience (int): Patience for early stopping.
    learning_rate (float): Learning rate for the optimizer.
    random_seed (int): Seed for random number generation.
    prop_split (float): Proportion of data to be used for training.
    ts (str): Timestamp or unique identifier for file naming.
    """

    def __init__(self, **kwargs):
        """
        The constructor for AutoencoderTrainer class.

        Parameters:
        **kwargs: Arbitrary keyword arguments.
        """
        # Initialize instance attributes
        self.dataset_path = kwargs['dataset_path']
        self.codings_size = kwargs['codings_size']
        self.epochs = kwargs['epochs']
        self.patience = kwargs['patience']
        self.learning_rate = kwargs['learning_rate']
        self.random_seed = kwargs['random_seed']
        self.prop_split = kwargs['prop_split']
        self.ts = kwargs['ts']
        self.output_loc = kwargs['output_loc']

        # Load and preprocess data
        self.ref_data = pd.read_csv(self.dataset_path, compression="gzip")
        self.df_train, self.df_validation = self.my_data_split()
        self.df, self.df_train, self.df_validation = self.preProcess()

        # Initialize the autoencoder model
        self.ae = self.model()

    def my_data_split(self):
        """
        Splits the data into training and validation sets based on unique genotypes.

        Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: DataFrames for training and validation data.
        """
        train_labels = []
        validation_labels = []
        train_indices = []
        validation_indices = []

        # Shuffle indices for random sampling
        indices = list(self.ref_data.index)
        shuffled_indices = shuffle(indices, random_state=self.random_seed)

        shuffled_df = self.ref_data.loc[shuffled_indices]

        # Determine split based on unique genotypes
        genotypes = list(shuffled_df["PlotID"])
        ratio = int(self.prop_split * len(genotypes))
        count = 0

        for geno, ID in zip(genotypes, list(shuffled_df.index)):
            if geno not in train_labels:
                train_labels.append(geno)
                train_indices.append(ID)
                count += 1
            if count == ratio:
                break

        # Remaining data for validation
        for geno, ID in zip(genotypes, list(shuffled_df.index)):
            if ID not in train_indices:
                validation_labels.append(geno)
                validation_indices.append(ID)

        df_train = shuffled_df.loc[train_indices]
        df_validation = shuffled_df.loc[validation_indices]

        return df_train, df_validation

    def preProcess(self):
        """
        Preprocesses the data by selecting relevant columns and filling missing values.

        Returns:
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: DataFrames of original, training, and validation data.
        """
        # Select wavelength columns
        wavelengths = [i for i in self.ref_data.columns if i.startswith("X")]

        df = self.ref_data[wavelengths]
        df_train = self.df_train[wavelengths]
        df_validation = self.df_validation[wavelengths]

        # Impute missing values with median
        for col in df_train.columns:
            df_train[col] = df_train[col].fillna(df_train[col].median())
        assert df_train.isnull().values.any() == False, 'Failed to remove all nan values from training set' 

        for col in df_validation.columns:
            df_validation[col] = df_validation[col].fillna(df_validation[col].median())
        assert df_validation.isnull().values.any() == False, 'Failed to remove all nan values from validation set' 

        return df, df_train, df_validation

    def model(self):
        """
        Builds the autoencoder model.

        Returns:
        keras.Model: The autoencoder model.
        """
        # Define the input layer
        inputs = keras.Input(shape=(self.df_train.shape[1],))

        # Encoder layers
        z = keras.layers.Flatten()(inputs)
        z = keras.layers.Dense(2200, activation="selu")(z)
        z = keras.layers.Dropout(0.3)(z)
        z = keras.layers.Dense(3000, activation="selu")(z)
        z = keras.layers.Dropout(0.3)(z)
        z = keras.layers.Dense(2024, activation="selu")(z)
        z = keras.layers.Dense(self.codings_size)(z)
        self.encoder = keras.Model(inputs=[inputs], outputs=[z])

        # Decoder layers
        decoder_inputs = keras.layers.Input(shape=[self.codings_size])
        x = keras.layers.Dense(1024, activation="selu")(decoder_inputs)
        x = keras.layers.Dense(1536, activation="selu")(x)
        x = keras.layers.Dense(2500, activation="selu")(x)
        x = keras.layers.Dense(2500, activation="selu")(x)
        x = keras.layers.Dense(self.df_train.shape[1], activation="tanh")(x)
        self.decoder = keras.Model(inputs=[decoder_inputs], outputs=[x])

        # Full autoencoder model
        codings = self.encoder(inputs)
        reconstructions = self.decoder(codings)
        self.ae = keras.Model(inputs=[inputs], outputs=[reconstructions])
    
        return self.ae

    def train(self):
        """
        Trains the autoencoder model and saves the training and validation losses, model weights, and latent variables.
        """
        # Compile the model
        self.ae.compile(optimizer=tf.keras.optimizers.SGD(learning_rate=self.learning_rate), loss="mae")
        
        # Callbacks for saving the model and early stopping
        checkpoint = ModelCheckpoint(
            f"{self.output_loc}/model_mae_{self.ts}.h5",
            verbose=1,
            monitor="val_loss",
            save_best_only=True,
            mode="auto")
        
        early_stopping = tf.keras.callbacks.EarlyStopping(
            monitor="val_loss",
            min_delta=0,
            patience=self.patience,
            verbose=0,
            mode="auto",
            baseline=None,
            restore_best_weights=True)
        
        # Fit the model on training data
        history = self.ae.fit(
            self.df_train,
            self.df_train,
            epochs=self.epochs,
            shuffle=True,
            validation_data=(self.df_validation, self.df_validation),
            callbacks=[checkpoint, early_stopping])
        
        # Save training and validation loss histories
        validation_loss = pd.DataFrame(history.history["val_loss"])
        validation_loss.to_csv(f"{self.output_loc}/validation_loss_{self.ts}.csv")
        training_loss = pd.DataFrame(history.history["loss"])
        training_loss.to_csv(f"{self.output_loc}/training_loss_{self.ts}.csv")
        
        # Generate and save latent variables
        latent_variables = self.encoder.predict(self.df)
        self.ae.save_weights(f"{self.output_loc}/weights_{self.ts}.h5")
        self.ae.save(f"{self.output_loc}/model_{self.ts}.h5")
        
        # Create and save a DataFrame with latent variables
        final_df = pd.DataFrame()
        for i in range(0, self.codings_size):
            name = "LV{}".format(i + 1)
            final_df[name] = latent_variables[:, i]
        final_df.insert(0, "PlotID", self.ref_data["PlotID"])
        final_df["trt"] = self.ref_data["trt"]
        final_df.to_csv(f"{self.output_loc}/LVs_Reflectance_Data_{self.ts}.csv", index=False)
