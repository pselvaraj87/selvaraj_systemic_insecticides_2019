import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

def plot_elim_heatmap():
    df = pd.read_csv("elim_data.csv")

    df["mda_coverage"][df["mda_coverage"]==-1]=0
    df["ivm_duration"][df["ivm_duration"]==-1]=0

    df.rename(columns={"mda_coverage": "MDA Coverage",
                      "ivm_duration": "Endectocide Duration"},inplace=True)

    plt.figure(figsize=(20,5))

    plt.subplot(131)
    sub_df = df[df["vc_pack"]=="none"].reset_index()
    foo = sub_df.pivot("Endectocide Duration", "MDA Coverage", "elim_eoy2_frac")
    sns.heatmap(foo, annot=True)
    plt.title("No ITN or IRS")

    plt.subplot(132)
    sub_df = df[df["vc_pack"]=="ITN only"].reset_index()
    foo = sub_df.pivot("Endectocide Duration", "MDA Coverage", "elim_eoy2_frac")
    sns.heatmap(foo, annot=True)
    plt.title("ITN only")

    plt.subplot(133)
    sub_df = df[df["vc_pack"]=="ITN and IRS"].reset_index()
    foo = sub_df.pivot("Endectocide Duration", "MDA Coverage", "elim_eoy2_frac")
    sns.heatmap(foo, annot=True)
    plt.title("ITN and IRS")

    plt.show()



def plot_elim_curves():
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
    plot_elim_curves()
    plot_elim_heatmap()