#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 14:05:10 2022

@author: robertoking
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


t1 = np.array([ 1.14992049e+00,  2.07734007e+00,  2.75926624e+00,  2.25918705e+00,
        3.85944045e+00,  2.33192584e+00,  2.83200503e+00,  3.50483884e+00,
        3.11386784e+00,  1.72273847e+00,  2.80472798e+00,  2.38647994e+00,
        3.03203670e+00,  2.31374114e+00,  3.11386784e+00,  3.67759347e+00,
        3.53211589e+00,  2.81382033e+00,  1.34995217e+00,  3.80488635e+00,
        2.66834275e+00,  1.97732424e+00,  1.13173580e+00,  2.40466463e+00,
        2.74108154e+00,  2.42284933e+00,  3.02294435e+00,  1.63181498e+00,
        3.82307105e+00,  1.93186249e+00,  2.32283349e+00,  1.40450626e+00,
        3.28662247e+00,  2.16826356e+00,  2.15007887e+00,  3.12296019e+00,
        2.81382033e+00,  3.05022140e+00,  3.44119240e+00,  3.24116073e+00,
        2.48649577e+00,  2.58651161e+00,  2.72289684e+00,  3.55030059e+00,
        3.55030059e+00,  3.38663831e+00,  2.40466463e+00,  3.58666998e+00,
        3.66850112e+00,  3.71396287e+00,  3.49574649e+00,  2.32283349e+00,
        2.45012638e+00,  3.11386784e+00,  2.33192584e+00,  2.36829524e+00,
        1.15901284e+00,  2.58651161e+00,  3.03203670e+00,  2.91383617e+00,
        2.13189417e+00,  2.38647994e+00,  1.65909203e+00,  3.17751428e+00,
        2.82291268e+00,  1.17719754e+00,  4.68684419e+00,  1.98641659e+00,
        2.44103403e+00,  2.35011054e+00,  2.55923456e+00,  1.86821605e+00,
        3.19569898e+00,  2.14098652e+00,  2.83200503e+00,  2.34101819e+00,
        3.11386784e+00,  1.95004719e+00,  2.96839026e+00,  2.84109738e+00,
        2.29555645e+00,  3.25934542e+00,  4.58901981e-01,  1.67727673e+00,
        1.75001552e+00,  4.01401038e+00,  2.21372531e+00,  2.98657496e+00,
        2.22281766e+00,  2.35011054e+00,  2.65925040e+00,  2.44103403e+00,
        2.80472798e+00,  2.02278598e+00,  3.22297603e+00,  2.35011054e+00,
        2.63197335e+00,  3.25025308e+00,  1.23175163e+00,  1.44087566e+00,
        1.95913954e+00,  2.96839026e+00,  2.15007887e+00,  2.09552477e+00,
        3.48665415e+00,  2.17735591e+00,  2.71380449e+00,  3.32299187e+00,
        3.12296019e+00,  1.60453794e+00,  1.34995217e+00,  3.70487052e+00,
        2.98657496e+00,  3.35026891e+00,  2.09552477e+00,  3.77760931e+00,
        3.54120824e+00,  5.86194865e-01,  2.06824773e+00,  3.07749845e+00,
        2.68652745e+00,  1.95004719e+00,  2.90474382e+00,  3.00475966e+00,
        3.30480717e+00,  1.98641659e+00,  1.54089149e+00,  2.63197335e+00,
        4.45953547e+00,  2.07734007e+00,  2.47740342e+00,  1.98641659e+00,
        2.78654328e+00,  1.78638491e+00,  3.25934542e+00,  2.70471214e+00,
        1.67727673e+00,  2.95929791e+00,  2.61378866e+00,  3.22297603e+00,
        2.56832691e+00,  2.09552477e+00,  2.85928208e+00,  2.32283349e+00,
        1.47724505e+00,  2.79563563e+00,  4.25041145e+00,  2.72289684e+00,
        2.43194168e+00,  2.68652745e+00,  3.25025308e+00,  2.77745094e+00,
        1.37722921e+00,  2.84109738e+00,  2.28646410e+00,  2.20463296e+00,
        3.94127159e+00,  1.94095484e+00,  2.21372531e+00,  1.85912370e+00,
        2.15917121e+00,  1.85003135e+00,  2.74108154e+00,  1.36813687e+00,
        3.24116073e+00,  3.53211589e+00,  2.40466463e+00,  2.95020556e+00,
        3.70487052e+00,  2.50468047e+00,  3.45028475e+00,  2.37738759e+00,
        3.43210005e+00,  2.28646410e+00,  4.25950380e+00,  3.04112905e+00,
        3.11386784e+00,  3.31389952e+00,  3.11386784e+00,  2.99566731e+00,
        1.38632156e+00,  2.44103403e+00,  2.75017389e+00,  3.18660663e+00,
        2.10461712e+00,  3.38663831e+00,  3.50483884e+00,  1.78638491e+00,
        1.78638491e+00,  2.43194168e+00,  3.29571482e+00,  2.01369363e+00,
        3.37754596e+00,  2.02278598e+00,  2.95020556e+00,  2.85928208e+00,
        9.86258214e-01,  9.40796470e-01,  3.32299187e+00,  2.02278598e+00,
        3.15932959e+00,  2.22281766e+00,  3.15023724e+00,  1.40450626e+00,
        2.85928208e+00,  3.61394703e+00,  2.45921873e+00,  3.85034810e+00,
        2.00460128e+00,  3.60485468e+00,  3.04112905e+00,  1.96823189e+00,
        4.23222675e+00,  3.51393119e+00,  3.66850112e+00,  3.20479133e+00,
        2.05915538e+00,  2.91383617e+00,  3.02294435e+00,  2.18644826e+00,
        2.49558812e+00,  1.98641659e+00,  2.84109738e+00,  3.11386784e+00,
        1.83184666e+00,  2.98657496e+00,  2.56832691e+00,  1.44087566e+00,
        2.65925040e+00,  3.98673333e+00,  2.57741926e+00,  3.21388368e+00,
        2.11370947e+00,  1.28630573e+00,  1.43178331e+00,  3.30480717e+00,
        2.75926624e+00,  2.09552477e+00,  1.42269096e+00,  1.98641659e+00,
        2.45012638e+00,  2.80472798e+00,  3.10477549e+00,  3.45028475e+00,
        2.95929791e+00,  4.81413708e+00,  1.93186249e+00,  7.77134190e-01,
        1.29539807e+00,  2.34101819e+00,  3.00475966e+00,  3.29571482e+00,
        2.93202087e+00,  2.43194168e+00,  1.46815270e+00,  2.68652745e+00,
        3.47756180e+00,  3.35936126e+00,  1.22265928e+00,  4.26859615e+00,
        1.95913954e+00,  4.23222675e+00,  2.65925040e+00,  1.39541391e+00,
        1.60453794e+00,  1.14082814e+00,  1.79547726e+00,  3.25025308e+00,
        2.82291268e+00,  3.30480717e+00,  3.45937710e+00,  2.37738759e+00,
        3.75942461e+00,  2.35011054e+00,  2.38647994e+00,  2.95020556e+00,
        3.66850112e+00,  2.08643242e+00,  3.34117656e+00,  3.10477549e+00,
        2.31374114e+00,  3.95036394e+00,  2.84109738e+00,  2.20463296e+00,
        2.80472798e+00,  2.70471214e+00,  2.89565147e+00,  3.69577817e+00,
        2.40466463e+00,  2.94111321e+00,  1.80456961e+00,  2.61378866e+00,
        3.62303938e+00,  3.58666998e+00,  3.54120824e+00,  1.27721338e+00,
        2.69561980e+00,  3.30480717e+00,  2.97748261e+00,  2.21372531e+00,
        2.96839026e+00,  3.11386784e+00,  3.10477549e+00,  2.71380449e+00,
        4.18676501e+00,  3.95036394e+00,  3.53211589e+00,  2.34101819e+00,
        2.85928208e+00,  3.31389952e+00,  1.85912370e+00,  2.25009470e+00,
        2.33192584e+00,  2.68652745e+00,  2.55923456e+00,  3.56848528e+00,
        2.75926624e+00,  3.83216340e+00,  2.20463296e+00,  2.32283349e+00,
        4.07765682e+00,  3.77760931e+00,  3.75033226e+00,  3.02294435e+00,
        3.10477549e+00,  3.05022140e+00,  2.15007887e+00,  2.09552477e+00,
        2.88655912e+00,  2.67743510e+00,  2.41375698e+00,  3.37754596e+00,
        2.25009470e+00,  2.13189417e+00,  2.81382033e+00,  2.16826356e+00,
        1.23175163e+00,  3.29571482e+00,  3.47756180e+00,  2.35920289e+00,
        2.68652745e+00,  3.47756180e+00,  1.20447459e+00,  2.06824773e+00,
        1.88640075e+00,  2.65925040e+00,  1.29539807e+00,  3.41391535e+00,
        3.63213173e+00,  2.77745094e+00,  2.71380449e+00,  3.07749845e+00,
        2.33192584e+00,  2.55014221e+00,  2.75926624e+00,  1.95004719e+00,
        2.01369363e+00,  2.93202087e+00,  2.50468047e+00,  1.55907619e+00,
        2.76835859e+00,  1.02262761e+00,  2.23191001e+00,  1.74092317e+00,
        3.10477549e+00,  2.28646410e+00,  2.01369363e+00,  2.26827940e+00,
        2.58651161e+00,  3.94127159e+00,  2.83200503e+00,  3.91399454e+00,
        3.41391535e+00,  2.80472798e+00,  1.45906035e+00,  2.42284933e+00,
        2.36829524e+00,  2.43194168e+00,  3.01385201e+00,  2.84109738e+00,
        2.69561980e+00,  2.91383617e+00,  4.02310273e+00,  2.38647994e+00,
        2.95020556e+00,  2.66834275e+00,  2.83200503e+00,  2.17735591e+00,
        2.59560396e+00,  2.15007887e+00,  2.71380449e+00,  3.15932959e+00,
        4.22532586e-01,  2.83200503e+00,  3.09568315e+00,  1.82275431e+00,
        2.89565147e+00,  2.15007887e+00,  2.72289684e+00,  3.86853280e+00,
        4.00491803e+00,  1.54998384e+00,  2.65015805e+00,  3.16842194e+00,
        2.50468047e+00,  1.70455377e+00,  3.33208421e+00,  1.45906035e+00,
        1.38632156e+00,  2.91383617e+00,  2.75926624e+00,  1.72273847e+00,
        2.06824773e+00,  3.98673333e+00,  3.31389952e+00,  1.88640075e+00,
        2.95929791e+00,  2.17735591e+00,  2.49558812e+00,  3.34117656e+00,
        2.65015805e+00,  2.75017389e+00,  1.53179914e+00,  1.59544559e+00,
        1.54998384e+00,  2.08643242e+00,  1.95004719e+00,  1.98641659e+00,
        3.65031642e+00,  2.06824773e+00,  1.68636907e+00,  4.34133494e+00,
        3.59576233e+00,  2.18644826e+00,  7.95318888e-01,  2.47740342e+00,
        3.98673333e+00,  4.08674917e+00,  1.99550894e+00,  2.83200503e+00,
        1.25902868e+00,  2.05006303e+00,  3.04112905e+00,  3.70487052e+00,
        2.01369363e+00,  2.97748261e+00,  3.32299187e+00,  2.85928208e+00,
        9.22611772e-01,  2.80472798e+00,  1.73183082e+00,  3.09568315e+00,
        1.78638491e+00,  3.04112905e+00,  3.62303938e+00, -4.80780988e-03,
        2.37738759e+00,  3.40482301e+00,  2.83200503e+00,  2.35920289e+00,
        2.79563563e+00,  2.67743510e+00,  3.95036394e+00,  2.29555645e+00,
        3.32299187e+00,  2.12280182e+00,  2.71380449e+00,  1.70455377e+00,
        1.62272263e+00,  3.23206838e+00,  2.45012638e+00,  2.91383617e+00,
        3.45028475e+00,  3.06840610e+00,  1.50452210e+00,  2.54104987e+00,
        1.45906035e+00,  2.38647994e+00,  2.75017389e+00,  2.34101819e+00,
        2.70471214e+00,  1.41359861e+00,  2.49558812e+00,  3.32299187e+00,
        2.28646410e+00,  2.95929791e+00,  1.74092317e+00,  2.88655912e+00,
        2.86837442e+00,  3.26843777e+00,  1.65909203e+00,  3.76851696e+00,
        2.35920289e+00,  2.37738759e+00,  2.29555645e+00,  2.98657496e+00,
        2.87746677e+00,  2.62288101e+00,  2.14098652e+00,  3.08659080e+00])

t2 = np.array([2.60469631, 3.33208421, 2.75926624, 1.68636907, 2.06824773,
       3.19569898, 3.51393119, 1.29539807, 1.88640075, 2.32283349,
       2.3046488 , 3.11386784, 1.83184666, 2.99566731, 1.28630573,
       1.90458545, 3.15023724, 2.6956198 , 1.53179914, 2.08643242,
       2.2864641 , 1.840939  , 3.78670166, 2.49558812, 3.25025308,
       3.36845361, 2.56832691, 1.79547726, 2.52286517, 2.2864641 ,
       2.00460128, 2.2864641 , 2.05915538, 1.54089149, 1.78638491,
       1.54089149, 3.32299187, 2.83200503, 3.33208421, 1.67727673,
       2.18644826, 2.56832691, 2.51377282, 3.37754596, 1.16810519,
       0.76804184, 2.08643242, 2.45012638, 2.33192584, 3.48665415,
       1.99550894, 2.85928208, 2.59560396, 2.6956198 , 2.35011054,
       2.53195752, 2.53195752, 2.12280182, 2.23191001, 3.05931375,
       2.08643242, 2.6410657 , 2.95020556, 3.26843777, 3.52302354,
       3.66850112, 1.04081231, 2.99566731, 1.61363028, 2.35920289,
       3.23206838, 2.43194168, 2.42284933, 1.53179914, 1.19538224,
       2.08643242, 1.93186249, 2.18644826, 2.46831107, 1.64090733,
       2.40466463, 3.49574649, 3.90490219, 2.59560396, 2.24100235,
       1.83184666, 2.05006303, 3.19569898, 2.34101819, 2.99566731,
       2.97748261, 1.71364612, 0.09520803, 2.48649577, 1.840939  ,
       3.30480717, 2.73198919, 1.99550894, 3.95036394, 2.32283349,
       2.00460128, 3.4230077 , 2.2500947 , 3.25934542, 3.82307105,
       3.0865908 , 2.47740342, 2.21372531, 1.22265928, 3.0684061 ,
       3.64122408, 3.50483884, 2.48649577, 2.20463296, 2.88655912,
       1.64090733, 2.55923456, 3.75942461, 2.80472798, 1.86821605,
       4.05947212, 2.65015805, 2.97748261, 1.66818438, 3.71396287,
       4.04128742, 1.38632156, 2.48649577, 1.39541391, 2.16826356,
       3.45028475, 3.61394703, 2.21372531, 2.88655912, 1.32267512,
       2.20463296, 2.25918705, 3.57757763, 3.45028475, 2.76835859,
       1.8954931 , 2.50468047, 3.60485468, 2.17735591, 2.17735591,
       2.90474382, 3.41391535, 2.2682794 , 2.74108154, 2.49558812,
       2.63197335, 1.86821605, 3.05931375, 3.00475966, 4.46862782,
       2.21372531, 3.59576233, 2.21372531, 2.38647994, 2.19554061,
       2.59560396, 1.61363028, 2.02278598, 3.35936126, 3.17751428,
       2.79563563, 4.39588903, 3.46846945, 2.55923456, 2.72289684,
       2.15917121, 4.05037977, 1.81366196, 3.52302354, 1.21356693,
       2.60469631, 1.62272263, 2.77745094, 2.27737175, 2.97748261,
       2.2682794 , 3.32299187, 2.93202087, 2.54104987, 3.04112905,
       2.16826356, 2.65015805, 3.00475966, 0.97716586, 1.96823189,
       2.49558812, 1.840939  , 1.04990466, 3.93217924, 2.44103403,
       2.76835859, 2.78654328, 2.95929791, 2.45921873, 3.75942461,
       2.73198919, 1.75910787, 2.15007887, 1.41359861, 1.01353526,
       4.47772017, 2.99566731, 3.0320367 , 3.15932959, 1.80456961,
       2.6592504 , 1.59544559, 1.4681527 , 2.53195752, 2.48649577,
       3.15023724, 2.73198919, 3.52302354, 2.03187833, 2.05006303,
       3.33208421, 3.31389952, 3.4411924 , 2.44103403, 2.95020556,
       1.56816854, 3.09568315, 3.25934542, 2.35920289, 2.62288101,
       2.46831107, 2.6774351 , 1.60453794, 2.84109738, 2.02278598,
       2.45012638, 3.4230077 , 1.51361445, 3.21388368, 3.02294435,
       3.35936126, 2.17735591, 1.90458545, 1.96823189, 3.60485468,
       3.12296019, 2.15917121, 0.57710252, 1.47724505, 2.49558812,
       2.83200503, 2.37738759, 1.62272263, 3.25934542, 3.4593771 ,
       1.8591237 , 2.20463296, 1.97732424, 4.2049497 , 2.65015805,
       1.95004719, 1.66818438, 2.58651161, 2.57741926, 2.45921873,
       2.33192584, 1.99550894, 3.4411924 , 2.91383617, 1.80456961,
       4.07765682, 1.20447459, 2.05915538, 2.6592504 , 2.70471214,
       2.82291268, 2.50468047, 2.13189417, 1.90458545, 1.86821605,
       2.58651161, 1.54089149, 2.80472798, 3.8321634 , 1.76820021,
       3.57757763, 1.38632156, 3.48665415, 2.84109738, 1.02262761,
       2.91383617, 3.12296019, 2.07734007, 3.25025308, 2.02278598,
       2.19554061, 1.94095484, 1.76820021, 4.07765682, 3.32299187,
       1.54089149, 2.6410657 , 3.22297603, 2.45921873, 2.88655912,
       2.34101819, 2.99566731, 3.65940877, 2.01369363, 0.53164077,
       2.37738759, 2.6410657 , 2.11370947, 2.82291268, 2.82291268,
       2.6410657 , 4.44135077, 2.46831107, 2.22281766, 0.668026  ,
       3.04112905, 2.45012638, 3.49574649, 1.51361445, 2.65015805,
       1.95004719, 0.78622654, 3.04112905, 3.35936126, 2.88655912,
       3.13205254, 1.34995217, 0.67711835, 2.97748261, 1.61363028,
       2.29555645, 2.68652745, 2.16826356, 2.98657496, 1.99550894,
       2.47740342, 1.29539807, 1.51361445, 1.70455377, 3.35026891,
       1.77729256, 3.91399454, 1.88640075, 2.81382033, 3.0865908 ,
       2.39557228, 2.95020556, 3.38663831, 3.48665415, 2.86837442,
       3.11386784, 2.2864641 , 3.13205254, 4.04128742, 2.61378866,
       2.15007887, 3.58666998, 3.13205254, 2.75017389, 2.52286517,
       2.84109738, 2.59560396, 1.79547726, 1.83184666, 2.57741926,
       0.82259593, 3.26843777, 2.62288101, 3.35026891, 0.84078063,
       3.31389952, 2.15917121, 2.27737175, 3.13205254, 1.04081231,
       2.57741926, 1.80456961, 1.5227068 , 1.69546142, 2.05006303,
       1.63181498, 2.24100235, 4.17767266, 3.0502214 , 3.4411924 ,
       4.21404205, 3.38663831, 1.68636907, 3.22297603, 1.67727673,
       2.87746677, 2.90474382, 1.98641659, 2.16826356, 3.35026891,
       2.20463296, 1.64090733, 1.00444291, 2.83200503, 1.88640075,
       2.42284933, 2.21372531, 2.84109738, 2.91383617, 2.2682794 ,
       2.42284933, 3.24116073, 3.0320367 , 2.35920289, 0.88624238,
       1.34085982, 1.1317358 , 2.66834275, 1.68636907, 2.22281766,
       1.08627405, 3.4775618 , 2.53195752, 2.55923456, 3.18660663,
       3.34117656, 3.00475966, 1.94095484, 2.06824773, 1.4681527 ,
       2.75017389, 2.58651161, 2.6410657 , 0.53164077, 1.14992049,
       2.58651161, 2.76835859, 3.43210005, 2.20463296, 2.2864641 ,
       3.67759347, 2.56832691, 1.80456961, 1.79547726, 1.58635324,
       2.75926624, 2.6410657 , 3.8503481 , 2.80472798, 3.04112905,
       2.35920289, 3.10477549, 2.39557228, 2.43194168, 2.80472798,
       2.00460128, 3.71396287, 1.79547726, 2.87746677, 3.78670166,
       2.92292852, 1.75001552, 2.54104987, 2.81382033, 2.03187833,
       2.11370947, 1.97732424, 3.84125575, 3.54120824, 2.95929791,
       3.55939294, 2.87746677, 2.61378866, 2.29555645, 2.84109738,
       4.75049063, 3.4411924 , 1.78638491, 1.75910787, 3.63213173,
       1.96823189, 2.87746677, 3.89580984, 1.28630573, 2.09552477,
       3.4775618 , 2.6956198 , 1.53179914, 2.04097068, 2.95929791,
       2.59560396, 2.58651161, 2.19554061, 2.3046488 , 2.71380449,
       2.92292852, 2.6774351 , 2.41375698, 3.10477549, 2.91383617])

t3 = np.array([2.92292852, 1.42269096, 1.42269096, 3.77760931, 3.19569898,
       2.79563563, 2.62288101, 2.43194168, 2.45012638, 3.75033226,
       2.43194168, 1.71364612, 3.12296019, 3.55030059, 0.67711835,
       2.52286517, 2.16826356, 3.4775618 , 1.14992049, 4.04128742,
       2.15917121, 2.52286517, 3.82307105, 2.15917121, 0.73167245,
       3.15932959, 2.97748261, 2.00460128, 1.9136778 , 2.61378866,
       3.80488635, 2.95020556, 2.21372531, 2.27737175, 2.2864641 ,
       1.37722921, 1.30449042, 0.668026  , 2.60469631, 2.75017389,
       2.61378866, 1.60453794, 3.09568315, 2.81382033, 1.99550894,
       3.43210005, 3.50483884, 3.8139787 , 3.95945628, 3.98673333,
       2.97748261, 1.5045221 , 2.85018973, 4.04128742, 3.28662247,
       2.89565147, 2.47740342, 4.25041145, 1.14992049, 3.18660663,
       4.40498138, 3.07749845, 2.41375698, 1.35904452, 2.15007887,
       2.13189417, 2.15007887, 4.32315024, 1.57726089, 1.34085982,
       2.57741926, 2.90474382, 2.35920289, 3.27753012, 3.77760931,
       2.82291268, 1.95913954, 3.84125575, 1.98641659, 2.10461712,
       3.0684061 , 2.00460128, 5.10509224, 1.38632156, 2.45921873,
       2.81382033, 2.84109738, 2.35011054, 3.04112905, 2.10461712,
       2.38647994, 1.9136778 , 2.61378866, 1.4863374 , 3.32299187,
       1.840939  , 1.81366196, 3.28662247, 1.98641659, 1.78638491,
       3.49574649, 2.35920289, 3.94127159, 3.93217924, 3.4775618 ,
       3.78670166, 2.84109738, 1.28630573, 2.6410657 , 2.35920289,
       2.08643242, 4.45953547, 3.0502214 , 1.27721338, 2.86837442,
       2.20463296, 3.22297603, 2.73198919, 3.55030059, 3.74123991,
       2.56832691, 2.50468047, 1.840939  , 2.62288101, 1.34995217,
       1.57726089, 3.0502214 , 1.21356693, 2.2864641 , 2.38647994,
       2.89565147, 2.34101819, 1.94095484, 3.55030059, 2.52286517,
       2.13189417, 4.31405789, 3.93217924, 2.75017389, 0.49527138,
       3.11386784, 3.0320367 , 2.18644826, 2.95020556, 4.43225842,
       2.50468047, 3.25934542, 1.449968  , 2.39557228, 2.82291268,
       3.11386784, 1.63181498, 2.06824773, 2.58651161, 3.05931375,
       2.95020556, 3.30480717, 2.76835859, 2.80472798, 2.65015805,
       2.2500947 , 1.72273847, 3.0684061 , 2.42284933, 2.62288101,
       3.25934542, 2.55923456, 2.65015805, 2.29555645, 2.58651161,
       1.78638491, 1.94095484, 3.8139787 , 2.65015805, 3.07749845,
       1.53179914, 2.00460128, 2.75926624, 1.92277014, 2.87746677,
       3.15932959, 2.22281766, 1.08627405, 3.25025308, 3.62303938,
       1.97732424, 3.63213173, 2.37738759, 3.65940877, 1.8591237 ,
       2.38647994, 3.61394703, 3.25934542, 3.55939294, 2.33192584,
       3.02294435, 2.08643242, 1.63181498, 2.78654328, 3.46846945,
       2.34101819, 2.93202087, 2.38647994, 3.67759347, 1.79547726,
       2.59560396, 1.8954931 , 1.72273847, 2.73198919, 3.76851696,
       3.12296019, 3.90490219, 2.13189417, 1.00444291, 1.40450626,
       3.55939294, 2.35920289, 3.38663831, 3.12296019, 2.05915538,
       3.14114489, 2.10461712, 2.16826356, 2.07734007, 2.15007887,
       2.56832691, 2.22281766, 2.85018973, 1.27721338, 1.82275431,
       1.26812103, 3.0684061 , 1.85003135, 3.0320367 , 2.55014221,
       2.75017389, 2.20463296, 2.77745094, 2.61378866, 1.95004719,
       2.51377282, 1.82275431, 2.35011054, 0.56801017, 2.72289684,
       3.02294435, 2.21372531, 1.42269096, 0.73167245, 1.94095484,
       1.62272263, 2.97748261, 3.62303938, 2.20463296, 3.97764098,
       2.37738759, 3.36845361, 1.64090733, 2.95929791, 2.52286517,
       2.6592504 , 2.93202087, 1.29539807, 2.6956198 , 2.54104987,
       2.15007887, 3.12296019, 2.39557228, 1.81366196, 3.15932959,
       3.41391535, 3.96854863, 2.73198919, 3.45028475, 2.38647994,
       1.15901284, 2.45921873, 2.6774351 , 3.95036394, 2.75926624,
       3.15023724, 2.61378866, 2.63197335, 1.8954931 , 1.5227068 ,
       1.73183082, 1.54998384, 2.95020556, 2.45012638, 2.15007887,
       2.81382033, 2.14098652, 1.44087566, 2.76835859, 3.67759347,
       2.47740342, 2.00460128, 1.88640075, 2.12280182, 2.95020556,
       3.11386784, 2.92292852, 2.58651161, 2.95020556, 2.92292852,
       3.0320367 , 2.41375698, 2.91383617, 1.78638491, 2.39557228,
       2.76835859, 2.78654328, 1.12264345, 2.20463296, 1.95004719,
       2.96839026, 3.75942461, 3.77760931, 2.38647994, 3.52302354,
       3.22297603, 1.36813687, 1.78638491, 2.44103403, 4.00491803,
       1.63181498, 2.09552477, 2.12280182, 2.20463296, 2.23191001,
       3.63213173, 1.90458545, 3.0320367 , 4.25041145, 1.10445875,
       2.66834275, 2.89565147, 2.07734007, 1.840939  , 2.55014221,
       0.60437956, 2.96839026, 2.44103403, 2.39557228, 2.40466463,
       2.88655912, 3.84125575, 1.73183082, 2.60469631, 2.88655912,
       1.53179914, 3.28662247, 2.48649577, 3.69577817, 2.6592504 ,
       2.45921873, 2.33192584, 3.99582568, 3.46846945, 3.15023724,
       2.45012638, 2.85928208, 2.04097068, 2.91383617, 2.95929791,
       3.41391535, 2.11370947, 1.75001552, 2.44103403, 1.55907619,
       3.84125575, 2.12280182, 3.25025308, 2.6410657 , 3.31389952,
       2.80472798, 2.45921873, 3.62303938, 3.23206838, 3.14114489,
       2.08643242, 1.8591237 , 2.56832691, 3.29571482, 2.34101819,
       2.15917121, 3.0320367 , 1.79547726, 2.63197335, 1.69546142,
       4.56864366, 2.51377282, 2.35011054, 1.10445875, 3.09568315,
       2.77745094, 0.98625821, 3.95945628, 3.36845361, 2.57741926,
       1.25902868, 3.48665415, 3.09568315, 3.37754596, 2.52286517,
       3.63213173, 4.26859615, 1.99550894, 2.05006303, 3.27753012,
       2.40466463, 2.39557228, 3.0865908 , 4.23222675, 3.20479133,
       1.53179914, 1.75910787, 3.10477549, 3.60485468, 2.63197335,
       3.27753012, 2.97748261, 2.59560396, 1.70455377, 3.61394703,
       3.50483884, 2.59560396, 1.98641659, 3.43210005, 1.62272263,
       3.24116073, 2.11370947, 2.85928208, 3.13205254, 3.21388368,
       2.2682794 , 3.74123991, 2.07734007, 2.55923456, 2.37738759,
       2.74108154, 1.8954931 , 2.37738759, 2.17735591, 4.42316608,
       3.11386784, 2.66834275, 2.21372531, 3.11386784, 1.75001552,
       3.30480717, 4.62319775, 1.78638491, 3.88671749, 2.87746677,
       3.71396287, 2.75926624, 2.2682794 , 2.57741926, 2.32283349,
       2.16826356, 3.25025308, 2.97748261, 0.60437956, 2.95020556,
       3.01385201, 2.00460128, 2.84109738, 3.56848528, 1.08627405,
       1.18628989, 2.55923456, 2.02278598, 1.8591237 , 2.09552477,
       2.54104987, 2.6592504 , 0.78622654, 2.13189417, 2.97748261,
       1.94095484, 4.70502889, 2.41375698, 2.62288101, 2.54104987,
       3.38663831, 2.25918705, 3.49574649, 2.47740342, 1.80456961,
       1.88640075, 1.96823189, 3.25934542, 1.35904452, 3.14114489,
       3.11386784, 2.50468047, 2.65015805, 2.76835859, 3.25025308])

t4 = np.array([2.46831107, 3.21388368, 3.22297603, 3.13205254, 2.43194168,
       2.72289684, 2.60469631, 3.4411924 , 1.22265928, 1.93186249,
       3.18660663, 2.59560396, 1.81366196, 1.51361445, 3.22297603,
       2.88655912, 2.21372531, 1.39541391, 1.30449042, 2.2864641 ,
       3.63213173, 3.29571482, 2.98657496, 2.56832691, 2.44103403,
       3.01385201, 1.37722921, 2.61378866, 1.47724505, 3.09568315,
       1.88640075, 1.75910787, 2.83200503, 2.95020556, 2.56832691,
       1.86821605, 1.95913954, 2.54104987, 1.56816854, 2.05006303,
       2.05006303, 2.59560396, 3.53211589, 2.25918705, 3.34117656,
       2.37738759, 3.01385201, 3.74123991, 3.58666998, 2.88655912,
       3.98673333, 2.71380449, 4.45953547, 2.45012638, 1.04081231,
       3.31389952, 2.24100235, 1.93186249, 2.12280182, 4.25041145,
       2.58651161, 2.90474382, 1.81366196, 2.6956198 , 1.20447459,
       2.84109738, 2.36829524, 3.30480717, 3.24116073, 3.72305522,
       1.86821605, 3.56848528, 2.27737175, 1.75910787, 2.74108154,
       2.97748261, 3.10477549, 2.46831107, 3.15932959, 3.20479133,
       2.47740342, 2.21372531, 3.16842194, 2.00460128, 2.84109738,
       3.70487052, 2.79563563, 2.41375698, 3.22297603, 3.32299187,
       2.06824773, 3.15932959, 2.53195752, 3.55939294, 2.81382033,
       2.63197335, 1.22265928, 2.17735591, 2.93202087, 1.02262761,
       2.51377282, 2.05006303, 3.13205254, 2.53195752, 2.57741926,
       3.24116073, 1.97732424, 1.17719754, 2.15917121, 3.73214756,
       2.53195752, 2.23191001, 3.74123991, 1.33176747, 2.77745094,
       2.12280182, 2.06824773, 2.85018973, 3.20479133, 2.44103403,
       2.27737175, 1.72273847, 2.76835859, 3.32299187, 2.71380449,
       3.14114489, 2.63197335, 4.53227426, 3.88671749, 3.82307105,
       3.68668582, 2.59560396, 2.49558812, 2.62288101, 2.3046488 ,
       1.85003135, 3.8503481 , 4.29587319, 3.4775618 , 2.2864641 ,
       2.82291268, 3.4230077 , 3.35936126, 2.45921873, 2.80472798,
       1.8591237 , 0.62256426, 3.8139787 , 2.14098652, 2.60469631,
       3.37754596, 2.47740342, 2.05915538, 0.45890198, 3.05931375,
       0.57710252, 3.62303938, 3.0865908 , 4.00491803, 2.85018973,
       2.57741926, 3.35026891, 1.62272263, 4.28678084, 3.15023724,
       2.75017389, 0.99535056, 2.89565147, 1.54998384, 2.79563563,
       2.2500947 , 2.78654328, 1.4681527 , 2.83200503, 2.87746677,
       2.6956198 , 1.88640075, 1.8954931 , 1.449968  , 3.45028475,
       2.24100235, 2.6592504 , 2.68652745, 2.37738759, 2.91383617,
       3.8503481 , 2.40466463, 3.36845361, 4.12311856, 1.90458545,
       4.16858031, 1.64090733, 2.52286517, 2.55014221, 4.11402622,
       1.79547726, 2.75926624, 3.55939294, 2.97748261, 2.85928208,
       2.42284933, 2.77745094, 3.4775618 , 2.33192584, 2.78654328,
       1.24993633, 1.74092317, 3.0684061 , 3.33208421, 2.92292852,
       3.00475966, 2.50468047, 0.36797849, 2.79563563, 2.37738759,
       0.94079647, 2.37738759, 3.20479133, 3.4775618 , 3.15023724,
       3.29571482, 4.31405789, 4.53227426, 2.41375698, 2.94111321,
       0.55891782, 4.06856447, 2.76835859, 4.04128742, 3.43210005,
       2.6774351 , 1.95913954, 1.53179914, 2.10461712, 2.92292852,
       2.27737175, 2.6774351 , 2.66834275, 2.33192584, 3.56848528,
       3.29571482, 2.35011054, 4.47772017, 1.64999968, 2.27737175,
       1.75001552, 2.15917121, 2.01369363, 2.54104987, 1.94095484,
       2.15007887, 3.12296019, 2.29555645, 2.15007887, 3.60485468,
       1.5045221 , 3.02294435, 2.02278598, 3.82307105, 2.45012638,
       2.81382033, 2.99566731, 2.15917121, 1.64090733, 2.35011054,
       1.67727673, 3.04112905, 1.32267512, 1.34995217, 1.08627405,
       3.16842194, 3.64122408, 1.73183082, 1.20447459, 2.10461712,
       2.36829524, 2.77745094, 3.66850112, 4.10493387, 1.77729256,
       2.45921873, 3.95036394, 3.35936126, 2.60469631, 3.21388368,
       4.01401038, 2.6592504 , 3.48665415, 4.12311856, 4.27768849,
       3.09568315, 2.78654328, 3.25934542, 3.4411924 , 2.48649577,
       3.4230077 , 2.95929791, 2.21372531, 3.71396287, 1.57726089,
       2.70471214, 2.91383617, 3.15023724, 2.22281766, 1.95913954,
       2.29555645, 1.24993633, 2.50468047, 3.26843777, 3.32299187,
       2.50468047, 3.20479133, 1.64090733, 2.57741926, 2.98657496,
       1.8773084 , 3.51393119, 3.05931375, 2.6592504 , 3.40482301,
       2.17735591, 2.39557228, 1.69546142, 1.95004719, 3.31389952,
       0.89533473, 2.60469631, 2.35920289, 3.12296019, 1.19538224,
       2.15007887, 2.66834275, 1.5227068 , 2.09552477, 2.55014221,
       3.68668582, 2.36829524, 4.02310273, 2.73198919, 1.92277014,
       3.66850112, 2.05915538, 2.01369363, 3.0865908 , 2.55923456,
       1.8954931 , 1.82275431, 3.11386784, 1.38632156, 2.00460128,
       1.49542975, 2.95929791, 3.61394703, 2.57741926, 2.91383617,
       2.70471214, 2.12280182, 3.01385201, 2.17735591, 2.36829524,
       2.16826356, 3.12296019, 1.21356693, 2.00460128, 2.11370947,
       3.34117656, 2.82291268, 2.84109738, 1.78638491, 3.25934542,
       2.34101819, 1.75001552, 2.19554061, 3.16842194, 3.4411924 ,
       2.93202087, 3.10477549, 2.99566731, 4.51408956, 1.95913954,
       2.77745094, 2.99566731, 3.37754596, 2.6410657 , 2.96839026,
       2.58651161, 1.67727673, 1.93186249, 2.10461712, 1.64090733,
       3.04112905, 2.50468047, 2.97748261, 1.34995217, 1.93186249,
       2.48649577, 2.23191001, 1.66818438, 1.8954931 , 2.38647994,
       2.46831107, 0.61347191, 3.77760931, 2.72289684, 3.4775618 ,
       2.75926624, 1.74092317, 2.66834275, 2.6774351 , 3.91399454,
       3.17751428, 2.08643242, 3.0320367 , 2.75926624, 2.45921873,
       2.88655912, 2.61378866, 0.97716586, 1.95913954, 0.81350359,
       1.02262761, 3.68668582, 2.04097068, 2.74108154, 1.98641659,
       3.35026891, 3.35026891, 1.39541391, 3.66850112, 2.39557228,
       2.86837442, 2.50468047, 3.75942461, 3.18660663, 3.30480717,
       1.5045221 , 1.99550894, 3.15932959, 2.84109738, 3.25934542,
       3.49574649, 2.31374114, 3.04112905, 2.82291268, 2.36829524,
       3.19569898, 3.46846945, 2.61378866, 1.449968  , 2.6410657 ,
       3.68668582, 1.93186249, 3.13205254, 3.12296019, 2.46831107,
       2.18644826, 0.67711835, 1.45906035, 2.08643242, 1.44087566,
       2.33192584, 2.27737175, 2.35920289, 3.58666998, 2.45012638,
       2.79563563, 2.52286517, 1.93186249, 2.2864641 , 1.25902868,
       2.88655912, 3.19569898, 2.6956198 , 2.88655912, 2.84109738,
       1.449968  , 4.02310273, 3.29571482, 3.4411924 , 2.92292852,
       3.38663831, 2.93202087, 2.33192584, 3.21388368, 2.81382033,
       2.2864641 , 1.20447459, 3.4775618 , 2.08643242, 2.16826356,
       2.65015805, 2.49558812, 2.96839026, 0.99535056, 2.88655912,
       3.59576233, 1.72273847, 2.74108154, 1.70455377, 3.29571482])



#%%
Nd  = 1000          #Data sample size
Nmc = 10000          #MCarlo sample size
myK = 10        #No of nearest neighbours checked
NoofTs = 2000
T_ValsN1000 = np.concatenate((t1,t2,t3,t4))
T_ValsN100 = np.array(pd.read_csv('MSM100'))
T_ValsN200 = np.array(pd.read_csv('MSM200'))
T_ValsN500 = np.array(pd.read_csv('MSM500'))

T_ValsN100 = [i[1] for i in T_ValsN100]
T_ValsN200 = [i[1] for i in T_ValsN200]
T_ValsN500 = [i[1] for i in T_ValsN500]

Tabove = [i for i in T_ValsN100 if i > 1.28]
rejectionpower = len(Tabove)/len(T_ValsN1000)
print('Rejection Power =',rejectionpower)

maxheight = 190


plt.hist(T_ValsN200,bins='auto',density = False,histtype = 'step',lw=2,label='$N_{d} = 100$')
plt.hist(T_ValsN100,bins='auto',density = False,histtype = 'step',lw=2,label='$N_{d} = 200$')
plt.hist(T_ValsN500,bins='auto',density = False,histtype = 'step',lw=2,label='$N_{d} = 500$')
plt.hist(T_ValsN1000,bins='auto',density = False,histtype = 'step',lw=2,label='$N_{d} = 1000$')

plt.plot([1.28,1.28],[0,maxheight],color='black',ls='dashed',label = '10% significance line')
plt.plot([1.64,1.64],[0,maxheight],color='black',ls='solid',label = '5% significance line')
plt.xlabel('Pull-T')
plt.ylabel('Frequency')
plt.legend()
plt.ylim(0,maxheight)
# plt.title('Mixed Sampled Method Pull-T Distribution, Neut v Nuwro,\n Hist Events = %s, $ N_{d} = 10N_{mc}, N_{k} = %s $' %(NoofTs,myK))
# plt.savefig('MixedSamp1')













