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
                     number_of_samples_exponent=11,
                     seed=0):

    number_of_samples = 2**number_of_samples_exponent
    m = number_of_samples_exponent
    sampler = qmc.Sobol(d=input_dim, scramble=True, seed=seed)
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
    uniaxial_strains = jnp.linspace(-0.3, 0, 50)
    data_samples = jax.vmap(uniaxial_test)(uniaxial_strains)
    return data_samples

def get_shear_test_data():
    shear_strains = jnp.linspace(-0.42, 0.42, 50)
    data_samples = jax.vmap(shear_test)(shear_strains)
    return data_samples




if __name__ == '__main__':

    # single rve per volume fraction
    volume_fractions = [0.1, 0.2, 0.3, 0.4, 0.5]
    seed_for_volume_fractions = [10, 20, 30, 40, 50]

    full_input_dataframe = pd.DataFrame(columns=['vol_req', 'seed', 'dU11', 'dU12', 'dU21', 'dU22'])

    for vf, seed in zip(volume_fractions, seed_for_volume_fractions):
        input_data_space = get_data_samples(input_dim=3,
                                            lower_bounds=[0.75, 0.75, 0],
                                            upper_bounds=[1.75, 1.75, np.pi/2],
                                            number_of_samples_exponent=10,
                                            seed=seed)
        
        # add two columns for volume fraction and seed
        input_data_space = np.hstack((np.full((input_data_space.shape[0], 1), vf),
                                        np.full((input_data_space.shape[0], 1), seed),
                                        input_data_space))
        
        input_dataframe = pd.DataFrame(input_data_space,
                                    columns=['vol_req', 'seed', 'dU11', 'dU12', 'dU21', 'dU22']).round(4)
        full_input_dataframe = pd.concat([full_input_dataframe, input_dataframe], axis=0)
    
    # reset index and keep it in the csv file
    full_input_dataframe = full_input_dataframe.reset_index(drop=True)
    
    full_input_dataframe.to_csv('doe_rubber_cb_all_vf.csv', index=True)
