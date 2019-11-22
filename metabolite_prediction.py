import numpy as np
import os
from util import write_data

def predict_metabolite_from_config(model_config, samples, record_mets=True):
    weightObservations = model_config["observed_weight_vector"]
    return predict_metabolties(samples,
                        np.asarray(model_config["eta"]),
                        np.asarray(model_config["tau"]),
                        model_config["mu"],
                        model_config["gamma"],
                        np.asarray(weightObservations),
                        model_config["compounds"],
                        record_mets)


def predict_metabolties(samples, eta, tau, mu, gamma, weight_observations, compounds, record_mets=True):
    if (np.sum(samples) == 0):
        print("PUMA ERROR - Metabolite identification: samples are all zero")
        exit()

    c = np.dot(samples, eta)
    phi = 1 - np.array(np.exp(np.log(1 - mu) * c))
    b = np.log(1 - gamma * phi)
    psi = np.dot(b, np.dot(tau, tau.transpose())) - b
    v = np.dot(weight_observations, tau.transpose())
    m_one = v * (phi * (1 - ((1 - gamma) * np.exp(psi)))) + (1 - v) * (phi * (1 - gamma))
    m_zero = v * ((1 - phi) * (1 - np.exp(psi))) + (1 - v) * (1 - phi)

    mean_mets_activity = np.mean(m_one, axis=0) / (np.mean(m_one, axis=0) + np.mean(m_zero, axis=0))
    print("mean_mets_activity_PUMA_detected:", list(mean_mets_activity))

    # write out m_one, m_zero, and mean_mets_activity files:
    if record_mets:
        outdata_dir = os.environ['PUMA_OUTPUT_DATA']
        metabolite_prediction_output = os.path.join(outdata_dir, 'metabolite_prediction_output.xlsx')
        mean_mets_activity = np.squeeze(mean_mets_activity).reshape(1, -1)
        write_data(mean_mets_activity, metabolite_prediction_output, sheetname="results", header=compounds)

    n_active_metabolites = len([metabolite for metabolite in mean_mets_activity[0] if metabolite >= 0.5])
    print("number_active_metabolites_PUMA_detected:", n_active_metabolites)
    active_mets_indices = np.nonzero(mean_mets_activity[0] >= 0.5)[0]
    active_mets_ID = [compounds[index] for index in active_mets_indices]
    print("active_metabolites_PUMA_detected:", active_mets_ID)
    return





