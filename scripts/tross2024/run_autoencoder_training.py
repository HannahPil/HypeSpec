#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 16:07:27 2022

@author: mtross
"""
import os
import sys
import argparse
import time
sys.path.insert(0, f'{os.getcwd()}/')
from autoencoder_trainer import autoencoder_trainer
ts = time.strftime("%m_%d_%Y__%H")

# Set up argument parsing for command line options
parser = argparse.ArgumentParser(description='Training autoencoder')

# Define command line arguments that the script accepts
parser.add_argument('--dataset', help='Hyperspectral leaf reflectance dataset. Should be a compressed "gzip" file')
parser.add_argument('--num_variables', help='Number of latent variables (reduced dimensions)', type=int, default=10)
parser.add_argument('--epochs', help='Total number of training epochs', type=int, default=1000)
parser.add_argument('--patience', help='Patience until early stopping', type=int, default=50)
parser.add_argument('--learning_rate', help='Rate at which model is adapted to the problem. Value ranges between 0 and 1.', type=float, default=0.1)
parser.add_argument('--random_seed', help='Random seed for model weight initialization', type=int, default=0)
parser.add_argument('--prop_split', help='Proportion split of data for training set. Remaining portion allocated to validation data', type=float, default=0.8)
parser.add_argument('--output_loc', help='Location of output files', type = str, default= os.getcwd())

# Parse the arguments provided by the user
args = parser.parse_args()

# Extract arguments and assign them to variables
dataset_path = args.dataset
codings_size = args.num_variables
epochs = args.epochs
patience = args.patience
learning_rate = args.learning_rate
random_seed = args.random_seed
prop_split = args.prop_split
output_loc = args.output_loc

# Initialize the AutoencoderTrainer with the provided arguments
trainer = autoencoder_trainer.AutoencoderTrainer(dataset_path=dataset_path,
                             codings_size=codings_size,
                             epochs=epochs,
                             patience=patience,
                             learning_rate=learning_rate,
                             random_seed=random_seed,
                             prop_split=prop_split,
                             output_loc=output_loc,
                             ts=ts)

# Main execution
if __name__ == "__main__":
    # Train the model using the trainer instance
    trainer.train()
