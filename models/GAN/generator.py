# Generator for a GAN using fingerprints
# This will be used to generate data to train a VAE
# Which will then be used to create another GAN
# that can generate molecules from features
import os
import random as r
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
from keras.layers import Dense, Flatten, Dropout, Conv2D, MaxPooling2D, Conv2DTranspose, LeakyReLU, BatchNormalization, Reshape
from keras import losses
from sklearn.model_selection import train_test_split
from keras.models import Sequential
import tensorflow as tf

generator = Sequential([
    Dense(256, input_shape=256),
    BatchNormalization(),
    LeakyReLU(),
    Dense(512),
    BatchNormalization(),
    LeakyReLU(),
    Dense(1024),
    BatchNormalization(),
    LeakyReLU(),
    Dense(2048, activation='tanh')
])



