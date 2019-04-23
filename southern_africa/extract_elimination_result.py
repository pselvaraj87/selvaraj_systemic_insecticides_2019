import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from simtools.Analysis.AnalyzeManager import AnalyzeManager

from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.Analysis.BaseAnalyzers.SimulationDirectoryMapAnalyzer import SimulationDirectoryMapAnalyzer
from simtools.Utilities.Experiments import retrieve_experiment


class ExtractInfectionResult(BaseAnalyzer):
    def __init__(self):
        filenames = ['output/ReportMalariaFilteredCatchment.json']
        super().__init__(filenames=filenames)

    # FOR TESTING PURPOSES ONLY
    # def filter(self, simulation):
    #     return simulation.tags['Run_Number'] == 1


    def select_simulation_data(self, data, simulation):
        eoy1 = np.array(data[self.filenames[0]]["Channels"]["Infected"]["Data"][3*365-1])
        eoy2 = np.array(data[self.filenames[0]]["Channels"]["Infected"]["Data"][4*365-1])

        return_df = pd.DataFrame({"time": ["eoy1", "eoy2"],
                                  "infec": [eoy1, eoy2]})
        return_df["sim_id"] = simulation.id
        return return_df


    def combine(self, all_data):
        data_list = list(all_data.values())
        # clean = [x for x in data_list if x != None]
        clean = list(filter(None.__ne__, data_list))
        df = pd.concat(clean, ignore_index=True)
        return df

    def finalize(self, all_data):
        sim_data_full = self.combine(all_data)
        sim_data_full.to_csv("raw_infec_data.csv", index=False)
        return sim_data_full



def convert_infection_csv_to_elim():
    df = pd.read_csv("raw_infec_data.csv")
    sim_map = pd.read_csv("sim_map.csv")

    full = df.merge(sim_map, how="left", left_on="sim_id", right_on="id")

    # Set flag for whether elimination happened
    full["elim_eoy1"] = 0
    full["elim_eoy2"] = 0
    full["elim_eoy1"][np.logical_and(full["time"]=="eoy1", full["infec"]==0)] = 1
    full["elim_eoy2"][np.logical_and(full["time"]=="eoy2", full["infec"]==0)] = 1

    # Group by different intervention variants, summing over elimination flag
    grp_by_package = full.groupby(["vc_pack", "drug_pack", "mda_coverage", "ivm_duration","vc_coverage"]).agg({"elim_eoy1": "sum",
                                                                                                 "elim_eoy2": "sum"})

    # Convert to fractions out of 100 which eliminated
    grp_by_package["elim_eoy1_frac"] = grp_by_package["elim_eoy1"]/100
    grp_by_package["elim_eoy2_frac"] = grp_by_package["elim_eoy2"]/100

    grp_by_package.reset_index(inplace=True)
    grp_by_package.to_csv("elim_data.csv", index=False)


def plot_elim_heatmap():
    pass

def plot_elim_curves(y=1):
    if y == 1:
        chan = "elim_eoy1_frac"
    elif y == 2:
        chan = "elim_eoy2_frac"




    df = pd.read_csv("elim_data.csv")

    plt.figure(figsize=(5,15))

    # No vector control:
    vc_packs = ["none", "ITN only", "ITN and IRS"]
    for i in range(3):
        plt.subplot(3,1,i+1)
        sub_df = df[df["vc_pack"]==vc_packs[i]]
        zero_cov = np.array(sub_df[sub_df["mda_coverage"]==-1][chan])[0]
        lines_dict = {}
        lines_dict["mda_wo_ivm"] = np.array(sub_df[sub_df["drug_pack"]=="MDA without IVM"].sort_values(by=['mda_coverage'], ascending=True)[chan])
        for ivm in [14,30,60,90]:
            hold = sub_df[np.logical_and(sub_df["drug_pack"] == "MDA with IVM", sub_df["ivm_duration"] == ivm)]
            lines_dict["mda_w_ivm{}".format(ivm)] = np.array(hold.sort_values(by=['mda_coverage'], ascending=True)[chan])

        # Add zero-coverage point to all curves:
        for key in lines_dict.keys():
            lines_dict[key] = np.append(zero_cov, lines_dict[key])

        x = np.array([0.,0.2,0.4,0.6,0.8,1.0])
        plt.plot(x,lines_dict["mda_wo_ivm"], c="C0",label="mda_wo_ivm")
        plt.plot(x,lines_dict["mda_w_ivm14"], c="C1",label="mda_w_ivm14")
        plt.plot(x,lines_dict["mda_w_ivm30"], c="C2",label="mda_w_ivm30")
        plt.plot(x,lines_dict["mda_w_ivm60"], c="C3",label="mda_w_ivm60")
        plt.plot(x,lines_dict["mda_w_ivm90"], c="C4",label="mda_w_ivm90")
        plt.xlabel("MDA Coverage")
        plt.ylabel("Fraction of sims eliminating")
        plt.title(vc_packs[i])
        plt.legend()

    plt.show()



if __name__=="__main__":
    if True:
        # analyzer_list = [ExtractInfectionResult()]
        analyzer_list = [ExtractInfectionResult(),
                         SimulationDirectoryMapAnalyzer(save_file="sim_map.csv")]
        exp_list = ["520818ca-ae3b-e911-a2c5-c4346bcb7273"]

        am = AnalyzeManager(force_analyze=True)

        for exp_name in exp_list:
            am.add_experiment(retrieve_experiment(exp_name))
        for a in analyzer_list:
            am.add_analyzer(a)

        am.analyze()

    if True:
        convert_infection_csv_to_elim()


    if False:
        plot_elim_curves(y=2)
