import numpy as np
import pandas as pd
from scipy.stats import qmc


def get_displacement_gradient_from_principal_stretches(sample_space):
    l1 = sample_space[:, 0]
    l2 = sample_space[:, 1]
    theta = sample_space[:, 2]

    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta), np.cos(theta)]]).transpose(2, 0, 1)

    U_eig = np.array([[l1, np.zeros_like(l1)],
                      [np.zeros_like(l2), l2]]).transpose(2, 0, 1)
    
    Identity = np.array([[np.ones_like(l1), np.zeros_like(l1)],
                         [np.zeros_like(l2), np.ones_like(l2)]]).transpose(2, 0, 1)

    F = np.einsum('mij,mjk->mik', R, U_eig)

    dU = F - Identity

    return pd.DataFrame(dU.reshape(-1, 4), columns=['dU11', 'dU12', 'dU21', 'dU22'])

def get_data_samples(input_dim = 3, lower_bounds = [0.2, 0.2, 0], upper_bounds = [4, 4, np.pi], number_of_samples_exponent=11):

    number_of_samples = 2**number_of_samples_exponent
    m = number_of_samples_exponent
    sampler = qmc.Sobol(d=input_dim, scramble=False)
    uniform_samples = sampler.random_base2(m=m)
    assert uniform_samples.shape[0] >= number_of_samples
    data_samples = get_displacement_gradient_from_principal_stretches(qmc.scale(uniform_samples, lower_bounds, upper_bounds))
    return data_samples

if __name__ == '__main__':
    input_data_space = get_data_samples(input_dim=3,
                                    lower_bounds=[0.2, 0.2, 0],
                                    upper_bounds=[4, 4, np.pi],
                                    number_of_samples_exponent=12)
    input_data_space.to_csv('doe.csv', index=False)