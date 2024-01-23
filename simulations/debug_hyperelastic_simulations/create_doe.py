import numpy as np
import pandas as pd
from scipy.stats import qmc
import jax
jax.config.update('jax_enable_x64', True)
import jax.numpy as jnp


def get_displacement_gradient_from_principal_stretches(sample):
    l1 = sample[0]
    l2 = sample[1]
    theta = sample[2]

    R = jnp.array([[jnp.cos(theta), -jnp.sin(theta)],
                  [jnp.sin(theta), jnp.cos(theta)]])

    U_eig = jnp.diag(jnp.array([l1, l2]))

    Identity = jnp.eye(2)

    F = jnp.matmul(R, jnp.matmul(U_eig, R.T))

    dU = F - Identity

    return dU

def get_data_samples(input_dim=3,
                     lower_bounds=[0.2, 0.2, 0],
                     upper_bounds=[4, 4, np.pi],
                     number_of_samples_exponent=11):

    number_of_samples = 2**number_of_samples_exponent
    m = number_of_samples_exponent
    sampler = qmc.Sobol(d=input_dim, scramble=False)
    uniform_samples = sampler.random_base2(m=m)
    assert uniform_samples.shape[0] >= number_of_samples
    data_samples = jax.vmap(
        get_displacement_gradient_from_principal_stretches)(
            qmc.scale(uniform_samples, lower_bounds, upper_bounds)
            ).reshape(-1, 4)

    return data_samples

if __name__ == '__main__':

    input_data_space = get_data_samples(input_dim=3,
                                        lower_bounds=[2, 2, 0],
                                        upper_bounds=[4, 4, np.pi/2],
                                        number_of_samples_exponent=4)

    input_dataframe = pd.DataFrame(input_data_space,
                                   columns=['dU11', 'dU12', 'dU21', 'dU22']).round(4)
    # input_dataframe = input_dataframe.tail(-16).reset_index(drop=True)
    input_dataframe.to_csv('doe_1.csv', index=True)
