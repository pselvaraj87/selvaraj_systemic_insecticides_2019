import os
import sys
sys.path.append(os.path.dirname(__file__))

from simtools.Analysis.AnalyzeManager import AnalyzeManager
from summary_report_analyzer import SummaryAnalyzer
from inset_chart_analyzer import InsetAnalyzer


if __name__ == "__main__":

    experiments = {
                    # "IVM_final": "dd4ce507-c72f-e911-a2c5-c4346bcb7273",
                    # "No_interventions": "e9e19499-c52f-e911-a2c5-c4346bcb7273"
                    # "eSMC_IVM_everyone": "9609582b-7b3f-e911-a2c5-c4346bcb7273",
                    # "eSMC10_IVM": "a47b2382-7b3f-e911-a2c5-c4346bcb7273",
                    # "Ivermectin_excludingchildbearing": "1afccb61-8846-e911-a2c0-c4346bcb1554",
                    'SMC_10': '154263e1-c031-e911-a2c5-c4346bcb7273'
                   }

    for expt_name, exp_id in experiments.items():
        if 'Reference' in expt_name:
            am = AnalyzeManager(exp_list=exp_id,
                                analyzers=[
                                    SummaryAnalyzer(expt_name=expt_name,
                                                    report_names=['Daily_Report'],
                                                    sweep_variables=["Run_Number"]),
                                    InsetAnalyzer(expt_name=expt_name,
                                                  report_names=['InsetChart'],
                                                  sweep_variables=["Run_Number"])

                                ],
                                force_analyze=False)

        else:
            am = AnalyzeManager(exp_list=exp_id,
                                analyzers=[
                                    SummaryAnalyzer(expt_name=expt_name,
                                                    report_names=['Daily_Report'],
                                                    sweep_variables=["Run_Number", "Coverage", "Intervention_type"]),
                                    InsetAnalyzer(expt_name=expt_name,
                                                  report_names=['InsetChart'],
                                                  sweep_variables=["Run_Number", "Coverage", "Intervention_type"])

                                ],
                                force_analyze=False)

        print(am.experiments)
        am.analyze()

