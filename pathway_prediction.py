import numpy as np
import pymc3 as pm
import theano.tensor as tt
import os
from util import write_data
from scipy.stats import logistic

def pathway_prediction_from_modelConfig(model_config, record_samples=True):
    a_init = np.asarray(model_config["a_init"])
    a_init = np.reshape(a_init, [len(a_init)])
    return pathway_prediction(model_config["landa"],
                              a_init,
                              model_config["mu"],
                              model_config["gamma"],
                              np.asarray(model_config["eta"]),
                              np.asarray(model_config["tau"]),
                              np.asarray(model_config["observed_weight_vector"]),
                              model_config["pathway_dict"],
                              record_samples)

# for printing debugging traces
# https://docs.pymc.io/notebooks/howto_debugging.html

def pathway_prediction(landa, a_init, mu, gamma, eta, tau, observed_weight_vector, pathway_dict,
                       record_samples=True):
    number_of_pathways = np.size(eta, 0)
    number_of_metabolites = np.size(eta, 1)
    myModel = pm.Model()
    with myModel:

        landa_value = pm.Beta('landa_value', alpha=1, beta=1)
        # define prior
        a = pm.Bernoulli('a', p=landa_value, shape=number_of_pathways)  # 1 x p
        # define posterior:  p (w|a)
        l = pm.math.dot(a, eta)  # 1xf: number of pathways that can generate each metabolite f
        phi = 1 - tt.exp(tt.log(1 - mu) * l)  # 1xf: p(m_j = 1| a)
        psi = 1 - tt.exp(tt.dot(tt.log(1 - (gamma * phi)), tau))  # 1xk: p(w_k=1 | a)
        w = pm.Bernoulli('w', p=psi, observed=observed_weight_vector, shape=observed_weight_vector.shape)

        start_point = {'landa_value': landa, 'a': a_init.astype(np.int32)}
        step1 = pm.Metropolis([landa_value])
        step2 = pm.BinaryGibbsMetropolis([a])
        trace = pm.sample(draws=1000, step=[step1, step2], start=start_point, random_seed=42)

    landa_value_samples_logodds = trace.get_values(trace.varnames[0], burn=100)
    landa_value_samples = logistic.pdf(landa_value_samples_logodds)
    pathways_samples = trace.get_values(trace.varnames[1],  burn=100)

    mean_pathways_activity = np.mean(pathways_samples, axis=0)
    if record_samples:
        outdata_dir = os.environ['PUMA_OUTPUT_DATA']
        pathway_prediction_output = os.path.join(outdata_dir, 'pathway_prediction_output.xlsx')
        mean_pathways_activity_in_samples = np.squeeze(mean_pathways_activity).reshape(1, -1)
        write_data(mean_pathways_activity_in_samples, pathway_prediction_output, sheetname="samples",
                   header=pathway_dict["pathway"])

    print("mean_pathways_activity_PUMA_detected:", list(mean_pathways_activity))
    n_active_pathways = len(
        [pathway_activity for pathway_activity in np.mean(pathways_samples, axis=0) if pathway_activity >= 0.5])
    print("number_active_pathways [PUMA detected]:", n_active_pathways)
    active_pathways_indices = np.nonzero(mean_pathways_activity >= 0.5)[0]
    active_pathways_ID = [pathway_dict["pathway"][index] for index in active_pathways_indices]
    print("active_pathways_PUMA_detected:", active_pathways_ID)
    not_active_pathways_indices = np.nonzero(mean_pathways_activity < 0.5)[0]
    not_active_pathways_ID = [pathway_dict["pathway"][index] for index in not_active_pathways_indices]
    print("not_active_pathways_PUMA_detected:", not_active_pathways_ID)
    return pathways_samples


def main():
    pass


if __name__ == "__main__":
    main()
