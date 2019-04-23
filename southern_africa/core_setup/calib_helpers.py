from dtk.vector.species import set_larval_habitat
from gridded_sims.run.site import catch_3_yr_spline


def set_ento(cb, a_sc, f_sc, arab_times, arab_spline, funest_times, funest_spline):
    hab = {
        'arabiensis': {
            "LINEAR_SPLINE": {
                "Capacity_Distribution_Number_Of_Years": int(len(arab_times)/12),
                "Capacity_Distribution_Over_Time": {
                    "Times": arab_times,
                    "Values": arab_spline
                },
                "Max_Larval_Capacity": pow(10, a_sc)
            }
        },
        'funestus': {
            "LINEAR_SPLINE": {
                "Capacity_Distribution_Number_Of_Years": int(len(funest_times)/12),
                "Capacity_Distribution_Over_Time": {
                    "Times": funest_times,
                    "Values": funest_spline
                },
                "Max_Larval_Capacity": pow(10, f_sc)
            }
        }
    }

    set_larval_habitat(cb, hab)

def map_sample_to_model_input(cb, catch, sample):
    a_sc = sample['arabiensis_scale']
    f_sc = sample['funestus_scale']

    arab_times, arab_spline = catch_3_yr_spline(catch, "gambiae")
    funest_times, funest_spline = catch_3_yr_spline(catch, "funestus")

    hab = {
        'arabiensis': {
            "LINEAR_SPLINE": {
                "Capacity_Distribution_Number_Of_Years": 3,
                "Capacity_Distribution_Over_Time": {
                    # "Capacity_Distribution_Per_Year": {
                    "Times": arab_times,
                    "Values": arab_spline
                },
                "Max_Larval_Capacity": pow(10, a_sc)
            }
        },
        'funestus': {
            "LINEAR_SPLINE": {
                "Capacity_Distribution_Number_Of_Years": 3,
                "Capacity_Distribution_Over_Time": {
                    # "Capacity_Distribution_Per_Year": {
                    "Times": funest_times,
                    "Values": funest_spline
                },
                "Max_Larval_Capacity": pow(10, f_sc)
                # "Max_Larval_Capacity": pow(10,a_sc)/arab_funest_ratio
            }
        }
    }

    set_larval_habitat(cb, hab)

    return {"arab": a_sc,
            "funest": f_sc}