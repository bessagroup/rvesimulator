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


def uniaxial_test(uniaxial_strain):
    gl_strain = jnp.array([uniaxial_strain, 0, 0, 0]).reshape(2, 2)
    C = 2.0 * gl_strain + jnp.eye(2)
    S, V, D = jnp.linalg.svd(C)
    F = S @ jnp.diag(jnp.sqrt(V)) @ D
    U = F - jnp.eye(2)
    return U.reshape(-1)

def shear_test(shear_strain):
    gl_strain = jnp.array([0, shear_strain, shear_strain, 0]).reshape(2, 2)
    C = 2.0 * gl_strain + jnp.eye(2)
    S, V, D = jnp.linalg.svd(C)
    F = S @ jnp.diag(jnp.sqrt(V)) @ D
    U = F - jnp.eye(2)
    return U.reshape(-1)

def get_uniaxial_test_data():
    uniaxial_strains = jnp.linspace(-0.3, 2, 50)
    data_samples = jax.vmap(uniaxial_test)(uniaxial_strains)
    return data_samples

def get_shear_test_data():
    shear_strains = jnp.linspace(-0.42, 42, 50)
    data_samples = jax.vmap(shear_test)(shear_strains)
    return data_samples




if __name__ == '__main__':

    # input_data_space = get_data_samples(input_dim=3,
    #                                     lower_bounds=[0.5, 0.5, 0],
    #                                     upper_bounds=[2, 2, np.pi/2],
    #                                     number_of_samples_exponent=11)
    
    input_data_space = get_shear_test_data()

    input_dataframe = pd.DataFrame(input_data_space,
                                   columns=['dU11', 'dU12', 'dU21', 'dU22']).round(4)
    # input_dataframe = input_dataframe.tail(-128).reset_index(drop=True)
    input_dataframe.to_csv('doe_shear_test.csv', index=True)
