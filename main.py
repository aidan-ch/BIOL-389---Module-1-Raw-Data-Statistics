import statistics
import scipy.stats as stats

import pandas as pd

control_data = pd.read_excel("BIOL 389 - CONTROL ERG data.xlsx", sheet_name=None)
test_data = pd.read_excel("BIOL 389 - ERG data.xlsx", sheet_name=None)

class Sweep :
    def __init__(self, sweep_n):
        self.sweep = sweep_n.reset_index(drop=True)
        self.stim1 = Stimulus(self.sweep.loc[0, "Time 1(ms)"], self.sweep.loc[0, "Time 2(ms)"],
                              self.sweep.loc[0, "Y1"], self.sweep.loc[0, "Y2"])
        self.stim2 = Stimulus(self.sweep.loc[0, "Time 3(ms)"], self.sweep.loc[0, "Time 4(ms)"],
                              self.sweep.loc[0, "Y3"], self.sweep.loc[0, "Y4"])
        self.stim3 = Stimulus(self.sweep.loc[1, "Time 1(ms)"], self.sweep.loc[1, "Time 2(ms)"],
                              self.sweep.loc[1, "Y1"], self.sweep.loc[1, "Y2"])
        self.stim4 = Stimulus(self.sweep.loc[1, "Time 3(ms)"], self.sweep.loc[1, "Time 4(ms)"],
                              self.sweep.loc[1, "Y3"], self.sweep.loc[1, "Y4"])

class Stimulus :
    def __init__(self, t_init, t_final, v_init, v_final):
        self.t_init = t_init
        self.t_final = t_final
        self.v_init = v_init
        self.v_final = v_final
        self.delta_v = v_final-v_init
        self.delta_t = t_final-t_init

class Fly :
    # across all 4 sweeps, first 3 stims
    # key = condition
    def __init__(self):
        self.avg_recoveries = dict()
        self.avg_depolarizations = dict()
        self.is_control = False

control_flies = [Fly(), Fly(), Fly(), Fly(), Fly()]
test_flies = [Fly(), Fly(), Fly(), Fly()]



control_out = {
    "Condition": [], "Avg_depol" : [], "SD avg depol" : [], "U value depol" : [], "p value depol" : [], "Avg_depol stim1": [], "Avg_depol stim2": [], "Avg_depol stim3": [], "Avg_depol stim4": [],
    "%stims with transients": [], "U value avg recovery" : [], "p value avg recovery" : [], "Avg % recovery" : [], "Avg % recovery stim1" : [], "Avg % recovery stim2" : [],
    "Avg % recovery stim3" : [], "Avg % recovery stim4" : [], "Stim 1 recovery time" : [], "Stim 2 recovery time" : [],
    "Stim 3 recovery time" : [], "Stim 4 recovery time" : []
}
test_out = {
    "Condition": [], "Avg_depol" : [] , "SD avg depol" : [], "U value depol" : [], "p value depol" : [], "Avg_depol stim1": [], "Avg_depol stim2": [], "Avg_depol stim3": [], "Avg_depol stim4": [],
    "%stims with transients": [], "U value avg recovery" : [], "p value avg recovery" :[], "Avg % recovery" : [], "Avg % recovery stim1" : [], "Avg % recovery stim2" : [],
    "Avg % recovery stim3" : [], "Avg % recovery stim4" : [], "Stim 1 recovery time" : [], "Stim 2 recovery time" : [],
    "Stim 3 recovery time" : [], "Stim 4 recovery time" : []
}

# the value for the values for the flies in the dictionary will be a list
mannwhit_control_depol = {
    "Condition" : [], "Depol. average for each fly" : []
}
mannwhit_control_recovery = {
    "Condition" : [], "Recovery average for each fly" : []
}
mannwhit_test_depol={
    "Condition" : [], "Depol. average for each fly" : []
}
mannwhit_test_recovery={
    "Condition" : [], "Recovery average for each fly" : []
}

for sheet_name, df in control_data.items():
    fly_dataframes = [df[ df["Fly #"] == 1], df[ df["Fly #"] == 2], df[ df["Fly #"] == 3], df[ df["Fly #"] == 4], df[ df["Fly #"] == 5]]
    control_out["Condition"].append(sheet_name)


    # average receptor depolarization for the ith stimulus, across all sweeps
    avg_depol1 = 0
    avg_depol2 = 0
    avg_depol3 = 0
    avg_depol4 = 0

    percent_stim_transients = 0

    avg_percent_recovery_1 = 0
    avg_percent_recovery_2 = 0
    avg_percent_recovery_3 = 0
    avg_percent_recovery_4 = 0


    # 64 stims total (5 flies * 4 sweeps * 4 stims) --> each row in excel sheet is 2 stims
    control_out["%stims with transients"].append( (df["On/Off transients?"].sum(axis=0) * 2)/80 * 100)

    #list of all average depols of each fly (across all stims) to get SD
    list_avg_depols = list()
    list_avg_recovery = list()


    mannwhit_control_depol["Condition"].append(sheet_name)
    mannwhit_control_recovery["Condition"].append(sheet_name)


    # iterate over each fly
    for i in range(0,5):
        sweep_1 = Sweep(fly_dataframes[i][fly_dataframes[i]["Sweep"] == 1])
        sweep_2 = Sweep(fly_dataframes[i][fly_dataframes[i]["Sweep"] == 2])
        sweep_3 = Sweep(fly_dataframes[i][fly_dataframes[i]["Sweep"] == 3])
        sweep_4 = Sweep(fly_dataframes[i][fly_dataframes[i]["Sweep"] == 4])

        control_flies[i].avg_recoveries[sheet_name] = []
        control_flies[i].avg_depolarizations[sheet_name] = []

        fly_avg_depol=0
        fly_avg_recovery=0
        for sweep in [sweep_1, sweep_2, sweep_3, sweep_4]:

            avg_depol1 += sweep.stim1.delta_v/20
            avg_depol2 += sweep.stim2.delta_v/20
            avg_depol3 += sweep.stim3.delta_v /20
            avg_depol4 += sweep.stim4.delta_v /20

            # to calculate average depolarization of this individual fly --> used for SDev
            fly_avg_depol += (sweep.stim1.delta_v + sweep.stim2.delta_v + sweep.stim3.delta_v + sweep.stim4.delta_v)

            avg_percent_recovery_1 += (sweep.stim2.v_init - sweep.stim1.v_final) / (sweep.stim1.v_init - sweep.stim1.v_final) * 100/20
            avg_percent_recovery_2 += (sweep.stim3.v_init - sweep.stim2.v_final) / (sweep.stim2.v_init - sweep.stim2.v_final) * 100/20
            avg_percent_recovery_3 += (sweep.stim4.v_init - sweep.stim3.v_final) / (sweep.stim3.v_init - sweep.stim3.v_final) * 100/20

            fly_avg_recovery += (((sweep.stim2.v_init - sweep.stim1.v_final) / (sweep.stim1.v_init - sweep.stim1.v_final) * 100)
            + ((sweep.stim3.v_init - sweep.stim2.v_final) / (sweep.stim2.v_init - sweep.stim2.v_final) * 100)
            + ((sweep.stim4.v_init - sweep.stim3.v_final) / (sweep.stim3.v_init - sweep.stim3.v_final) * 100))

        # average over 4 sweeps, 4 stimuli each sweep
        list_avg_depols.append(fly_avg_depol/(4*4))

        # average over 4 sweeps, looking at first 3 stims per sweep
        list_avg_recovery.append(fly_avg_recovery/(4*3))

        control_flies[i].avg_depolarizations[sheet_name] = fly_avg_depol/16
        control_flies[i].avg_recoveries[sheet_name]=fly_avg_recovery/12
        control_flies[i].is_control = True

        avg_percent_recovery_4 += ( (sweep_2.stim1.v_init - sweep_1.stim4.v_final) / ( sweep_1.stim4.v_init - sweep_1.stim4.v_final)) * 100/15
        avg_percent_recovery_4 += ( (sweep_3.stim1.v_init - sweep_2.stim4.v_final) / ( sweep_2.stim4.v_init - sweep_2.stim4.v_final)) * 100/15
        avg_percent_recovery_4 += ( (sweep_4.stim1.v_init - sweep_3.stim4.v_final) / ( sweep_3.stim4.v_init - sweep_3.stim4.v_final)) * 100/15

    mannwhit_control_depol["Depol. average for each fly"].append(list_avg_depols)
    mannwhit_control_recovery["Recovery average for each fly"].append(list_avg_recovery)



    control_out["Avg_depol"].append((avg_depol1 + avg_depol2 + avg_depol3 + avg_depol4) / 4)

    control_out["Avg % recovery"].append((avg_percent_recovery_1 + avg_percent_recovery_2 + avg_percent_recovery_3
                                       +avg_percent_recovery_4) / 4)

    control_out["Avg_depol stim1"].append(avg_depol1)
    control_out["Avg_depol stim2"].append(avg_depol2)
    control_out["Avg_depol stim3"].append(avg_depol3)
    control_out["Avg_depol stim4"].append(avg_depol4)

    control_out["Avg % recovery stim1"].append(avg_percent_recovery_1)
    control_out["Avg % recovery stim2"].append(avg_percent_recovery_2)
    control_out["Avg % recovery stim3"].append(avg_percent_recovery_3)
    control_out["Avg % recovery stim4"].append(avg_percent_recovery_4)

    control_out["Stim 1 recovery time"].append(sweep_1.stim2.t_init - sweep_1.stim1.t_final)
    control_out["Stim 2 recovery time"].append(sweep_1.stim3.t_init - sweep_1.stim2.t_final)
    control_out["Stim 3 recovery time"].append(sweep_1.stim4.t_init - sweep_1.stim3.t_final)

    time_between_sweeps = fly_dataframes[0].loc[1, "Trace start"] - fly_dataframes[0].loc[0, "Trace start"]
    control_out["Stim 4 recovery time"].append(sweep_2.stim1.t_init - sweep_1.stim4.t_final + time_between_sweeps)

    control_out["SD avg depol"].append(statistics.stdev(list_avg_depols))

for sheet_name, df in test_data.items():
    fly_dataframes = [df[ df["Fly #"] == 1], df[ df["Fly #"] == 2], df[ df["Fly #"] == 3], df[ df["Fly #"] == 4]]
    test_out["Condition"].append(sheet_name)


    # average receptor depolarization for the ith stimulus, across all sweeps
    avg_depol1 = 0
    avg_depol2 = 0
    avg_depol3 = 0
    avg_depol4 = 0

    percent_stim_transients = 0

    avg_percent_recovery_1 = 0
    avg_percent_recovery_2 = 0
    avg_percent_recovery_3 = 0
    avg_percent_recovery_4 = 0


    # 64 stims total (4 flies * 4 sweeps * 4 stims) --> each row in excel sheet is 2 stims
    test_out["%stims with transients"].append( (df["On/Off transients?"].sum(axis=0) * 2)/64 * 100)

    list_avg_depols = list()

    # NOT INCLUDE stim 4 as recovery time is diff
    list_avg_recovery = list()

    mannwhit_test_depol["Condition"].append(sheet_name)
    mannwhit_test_recovery["Condition"].append(sheet_name)
    # iterate over each fly
    for i in range(0,4):
        sweep_1 = Sweep(fly_dataframes[i][fly_dataframes[i]["Sweep"] == 1])
        sweep_2 = Sweep(fly_dataframes[i][fly_dataframes[i]["Sweep"] == 2])
        sweep_3 = Sweep(fly_dataframes[i][fly_dataframes[i]["Sweep"] == 3])
        sweep_4 = Sweep(fly_dataframes[i][fly_dataframes[i]["Sweep"] == 4])

        test_flies[i].avg_recoveries[sheet_name] = []
        test_flies[i].avg_depolarizations[sheet_name] = []

        fly_avg_depol = 0
        fly_avg_recovery = 0
        for sweep in [sweep_1, sweep_2, sweep_3, sweep_4]:

            avg_depol1 += sweep.stim1.delta_v/16
            avg_depol2 += sweep.stim2.delta_v/16
            avg_depol3 += sweep.stim3.delta_v /16
            avg_depol4 += sweep.stim4.delta_v /16

            # to calculate average depolarization of this individual fly --> used for SDev
            fly_avg_depol += (sweep.stim1.delta_v + sweep.stim2.delta_v + sweep.stim3.delta_v + sweep.stim4.delta_v)

            fly_avg_recovery += (((sweep.stim2.v_init - sweep.stim1.v_final) / (sweep.stim1.v_init - sweep.stim1.v_final) * 100)
            + ((sweep.stim3.v_init - sweep.stim2.v_final) / (sweep.stim2.v_init - sweep.stim2.v_final) * 100)
            + ((sweep.stim4.v_init - sweep.stim3.v_final) / (sweep.stim3.v_init - sweep.stim3.v_final) * 100))

            avg_percent_recovery_1 += (sweep.stim2.v_init - sweep.stim1.v_final) / (sweep.stim1.v_init - sweep.stim1.v_final) * 100/16
            avg_percent_recovery_2 += (sweep.stim3.v_init - sweep.stim2.v_final) / (sweep.stim2.v_init - sweep.stim2.v_final) * 100/16
            avg_percent_recovery_3 += (sweep.stim4.v_init - sweep.stim3.v_final) / (sweep.stim3.v_init - sweep.stim3.v_final) * 100/16


        avg_percent_recovery_4 += ( (sweep_2.stim1.v_init - sweep_1.stim4.v_final) / ( sweep_1.stim4.v_init - sweep_1.stim4.v_final)) * 100/12
        avg_percent_recovery_4 += ( (sweep_3.stim1.v_init - sweep_2.stim4.v_final) / ( sweep_2.stim4.v_init - sweep_2.stim4.v_final)) * 100/12
        avg_percent_recovery_4 += ( (sweep_4.stim1.v_init - sweep_3.stim4.v_final) / ( sweep_3.stim4.v_init - sweep_3.stim4.v_final)) * 100/12

        # average across 4 sweeps and 4 stims
        list_avg_depols.append(fly_avg_depol/(4*4))

        # average across 4 sweeps and the first 3 stims
        list_avg_recovery.append(fly_avg_recovery/(4*3))

        test_flies[i].avg_depolarizations[sheet_name]=fly_avg_depol/16
        test_flies[i].avg_recoveries[sheet_name]=fly_avg_recovery/12

    mannwhit_test_depol["Depol. average for each fly"].append(list_avg_depols)
    mannwhit_test_recovery["Recovery average for each fly"].append(list_avg_recovery)


    test_out["Avg_depol"].append(statistics.mean(list_avg_depols))

    test_out["Avg % recovery"].append((avg_percent_recovery_1 + avg_percent_recovery_2 + avg_percent_recovery_3
                                       +avg_percent_recovery_4) / 4)

    test_out["Avg_depol stim1"].append(avg_depol1)
    test_out["Avg_depol stim2"].append(avg_depol2)
    test_out["Avg_depol stim3"].append(avg_depol3)
    test_out["Avg_depol stim4"].append(avg_depol4)

    test_out["Avg % recovery stim1"].append(avg_percent_recovery_1)
    test_out["Avg % recovery stim2"].append(avg_percent_recovery_2)
    test_out["Avg % recovery stim3"].append(avg_percent_recovery_3)
    test_out["Avg % recovery stim4"].append(avg_percent_recovery_4)

    test_out["Stim 1 recovery time"].append(sweep_1.stim2.t_init - sweep_1.stim1.t_final)
    test_out["Stim 2 recovery time"].append(sweep_1.stim3.t_init - sweep_1.stim2.t_final)
    test_out["Stim 3 recovery time"].append(sweep_1.stim4.t_init - sweep_1.stim3.t_final)

    time_between_sweeps = fly_dataframes[0].loc[1, "Trace start"] - fly_dataframes[0].loc[0, "Trace start"]

    test_out["Stim 4 recovery time"].append(sweep_2.stim1.t_init - sweep_1.stim4.t_final + time_between_sweeps)

    test_out["SD avg depol"].append(statistics.stdev(list_avg_depols))

for idx, condition_name in enumerate(mannwhit_test_depol["Condition"]):
    test_vals_depol = mannwhit_test_depol["Depol. average for each fly"][idx]
    control_vals_depol = mannwhit_control_depol["Depol. average for each fly"][idx]
    depol_mw = stats.mannwhitneyu(test_vals_depol, control_vals_depol, alternative="greater")

    # u value depol and p value depol for test and control flies are the same for each condition
    # just putting it on each excel sheet bc weird to have only on one
    test_out["U value depol"].append(depol_mw.statistic)
    test_out["p value depol"].append(depol_mw.pvalue)

    control_out["U value depol"].append(depol_mw.statistic)
    control_out["p value depol"].append(depol_mw.pvalue)

    test_vals_recovery = mannwhit_test_recovery["Recovery average for each fly"][idx]
    control_vals_recovery = mannwhit_control_recovery["Recovery average for each fly"][idx]
    recovery_mw = stats.mannwhitneyu(test_vals_recovery, control_vals_recovery, alternative="two-sided")

    control_out["U value avg recovery"].append(recovery_mw.statistic)
    control_out["p value avg recovery"].append(recovery_mw.pvalue)
    test_out["p value avg recovery"].append(recovery_mw.pvalue)
    test_out["U value avg recovery"].append(recovery_mw.statistic)


test_out_df = pd.DataFrame(test_out)
control_out_df = pd.DataFrame(control_out)

with pd.ExcelWriter("test_out_data_v4.xlsx", engine="xlsxwriter") as w:
    test_out_df.to_excel(w, index=False, sheet_name="Summary")

with pd.ExcelWriter("control_out_data_v4.xlsx", engine="xlsxwriter") as w:
    control_out_df.to_excel(w, index=False, sheet_name="Summary")


# ONLY USE STIMS 1-3 FOR THESE
# OFF time data analysis
    # does off time change receptor responses for a given stim length
    # compare the difference in response magnitude WITHIN a fly a fixed on time
        # calculating a delta V for one off time
        # e.g fix ON = 1s, compare A (2s off) and C (1s off), delta off = A - C voltage



# Map protocols into matched ON-time pairs that differ only by OFF time
# OFF=2s protocols: A, B, E
# OFF=1s protocols: C, D, F
OFF_PAIRS = [
    {"on_time_s": 1.0, "off2": "A", "off1": "C"},
    {"on_time_s": 0.5, "off2": "B", "off1": "D"},
    {"on_time_s": 0.1, "off2": "E", "off1": "F"},
]

def _sheet_lookup_by_protocol_and_color(sheet_dict):
    """
    Returns dict: (protocol_letter, color_lower) -> sheet_name
    Works with sheet names like "A blue", "Protocol A Blue", etc.
    """
    out = {}
    for sheet_name in sheet_dict.keys():
        s = sheet_name.strip().lower()
        # infer protocol as first A-F that appears as a standalone token or start
        proto = None
        for p in ["a", "b", "c", "d", "e", "f"]:
            if s.startswith(p) or f" {p} " in f" {s} ":
                proto = p.upper()
                break
        if proto is None:
            # fallback: first character
            proto = sheet_name.strip()[0].upper()

        color = "blue" if "blue" in s else ("white" if "white" in s else None)
        if color is None:
            continue
        if (proto, color) in out: raise ValueError("Duplicate mapping...")
        out[(proto, color)] = sheet_name
    return out


# so I have 2 protocols with different off times
# I want to find the change in RECOVERY PERCENTAGE between the different off times
# I want to find the change in AVERAGE DEPOLARIZATION between the different off times

# So I need to go through each fly, and go through each sweep, calculate average depolarization and then average recovery % (first 3 stims)

def delta_off_potential():
    ctrl_lookup = _sheet_lookup_by_protocol_and_color(control_data)
    test_lookup = _sheet_lookup_by_protocol_and_color(test_data)

    data = {"Fly #" : [], "Fly type" : [], "Color" : [], "On time(s)" : [], "Off2 protocol" : [], "Off1 protocol" : [],
            "Off2 depolarization" : [], "Off1 depolarization" : [], "Delta Off2 Off1 voltage" : []}

    switched_to_test_flies = False
    # fly number
    n = 1
    for fly in control_flies + test_flies:
        if fly.is_control == False and not switched_to_test_flies:
            n = 1
            switched_to_test_flies = True

        lookup = ctrl_lookup if not switched_to_test_flies else test_lookup

        for color in ["blue", "white"]:

            for pair in OFF_PAIRS:
                data["Fly #"].append(n)
                data["Fly type"].append("Control" if fly.is_control else "Test")

                data["On time(s)"].append(pair["on_time_s"])
                data["Off2 protocol"].append(pair["off2"])
                data["Off1 protocol"].append(pair["off1"])
                data["Color"].append(color)

                off2_depol = fly.avg_depolarizations[lookup[(pair["off2"], color)]]
                off1_depol = fly.avg_depolarizations[lookup[(pair["off1"], color)]]
                data["Off2 depolarization"].append(off2_depol)
                data["Off1 depolarization"].append(off1_depol)

                data["Delta Off2 Off1 voltage"].append(off2_depol - off1_depol)

        n+=1

    return data

def delta_off_recovery():
    ctrl_lookup = _sheet_lookup_by_protocol_and_color(control_data)
    test_lookup = _sheet_lookup_by_protocol_and_color(test_data)

    data = {"Fly #" : [], "Fly type" : [], "Color" : [], "On time(s)" : [], "Off2 protocol" : [], "Off1 protocol" : [],
            "Off2 recovery" : [], "Off1 recovery" : [], "Delta Off2 Off1 recovery" : []}

    switched_to_test_flies = False
    # fly number
    n = 1
    for fly in control_flies + test_flies:
        if fly.is_control == False and not switched_to_test_flies:
            n = 1
            switched_to_test_flies = True

        lookup = ctrl_lookup if not switched_to_test_flies else test_lookup

        for color in ["blue", "white"]:

            for pair in OFF_PAIRS:
                data["Fly #"].append(n)
                data["Fly type"].append("Control" if fly.is_control else "Test")

                data["On time(s)"].append(pair["on_time_s"])
                data["Off2 protocol"].append(pair["off2"])
                data["Off1 protocol"].append(pair["off1"])
                data["Color"].append(color)

                off2_recovery = fly.avg_recoveries[lookup[(pair["off2"], color)]]
                off1_recovery= fly.avg_recoveries[lookup[(pair["off1"], color)]]
                data["Off2 recovery"].append(off2_recovery)
                data["Off1 recovery"].append(off1_recovery)

                data["Delta Off2 Off1 recovery"].append(off2_recovery - off1_recovery)

        n+=1

    return data

data_off_potential = delta_off_potential()
data_off_recovery = delta_off_recovery()




potential_df = pd.DataFrame(data_off_potential)
recoveries_df = pd.DataFrame(data_off_recovery)


def average_off_across_conditions(potential_df, recoveries_df):

    data = {"Fly Type" : [], "Color" : [], "On time(s)" : [], "Off2 protocol" : [], "Off1 protocol" : [],
            "Mean Off2_potential - Off1_potential" : [], "Mean Off2_recovery - Off1_recovery" : [], "U value depol" : [],
            "p value depol" : [], "U value recovery" : [], "p value recovery" : []}

    ctrl_lookup = _sheet_lookup_by_protocol_and_color(control_data)
    test_lookup = _sheet_lookup_by_protocol_and_color(test_data)

    potential_on_time_1 = potential_df[potential_df["On time(s)"] == 1]
    potential_on_time_point5 = potential_df[potential_df["On time(s)"] == 0.5]
    potential_on_time_point1 = potential_df[potential_df["On time(s)"] == 0.1]

    recovery_on_time_1 = recoveries_df[recoveries_df["On time(s)"] == 1]
    recovery_on_time_point5 = recoveries_df[recoveries_df["On time(s)"] == 0.5]
    recovery_on_time_point1 = recoveries_df[recoveries_df["On time(s)"] == 0.1]


    for on_time, recovery_on_time in [(potential_on_time_1, recovery_on_time_1),(potential_on_time_point1, recovery_on_time_point1), (potential_on_time_point5, recovery_on_time_point5)]:

        for filtered, recovery_filtered in [(on_time[on_time["Color"] == "blue"].reset_index(drop=True), recovery_on_time[recovery_on_time["Color"] == "blue"].reset_index(drop=True)),
                                            (on_time[on_time["Color"] == "white"].reset_index(drop=True), recovery_on_time[recovery_on_time["Color"] == "white"].reset_index(drop=True))]:




            mean_delta_v_control = filtered[filtered["Fly type"] == "Control"]["Delta Off2 Off1 voltage"].mean()
            mean_delta_v_test = filtered[filtered["Fly type"] == "Test"]["Delta Off2 Off1 voltage"].mean()


            mean_delta_recovery_control = recovery_filtered[recovery_filtered["Fly type"] == "Control"]["Delta Off2 Off1 recovery"].mean()
            mean_delta_recovery_test = recovery_filtered[recovery_filtered["Fly type"] == "Test"]["Delta Off2 Off1 recovery"].mean()

            ctrl_vals_voltage = filtered[filtered["Fly type"]=="Control"]["Delta Off2 Off1 voltage"].dropna().to_numpy()
            test_vals_voltage = filtered[filtered["Fly type"]=="Test"]["Delta Off2 Off1 voltage"].dropna().to_numpy()

            delta_v_stat = stats.mannwhitneyu(test_vals_voltage, ctrl_vals_voltage, alternative="two-sided")

            ctrl_vals_recovery = recovery_filtered[recovery_filtered["Fly type"]=="Control"]["Delta Off2 Off1 recovery"].dropna().to_numpy()
            test_vals_recovery = recovery_filtered[recovery_filtered["Fly type"]=="Test"]["Delta Off2 Off1 recovery"].dropna().to_numpy()

            delta_recovery_stat = stats.mannwhitneyu(ctrl_vals_recovery, test_vals_recovery, alternative="two-sided")


            data["Fly Type"].append("Control")
            data["Color"].append(filtered.loc[0, "Color"])
            data["On time(s)"].append(filtered.loc[0, "On time(s)"])
            data["Off2 protocol"].append(filtered.loc[0, "Off2 protocol"])
            data["Off1 protocol"].append(filtered.loc[0, "Off1 protocol"])

            data["Mean Off2_potential - Off1_potential"].append(mean_delta_v_control)
            data["Mean Off2_recovery - Off1_recovery"].append(mean_delta_recovery_control)

            data["U value depol"].append(delta_v_stat.statistic)
            data["p value depol"].append(delta_v_stat.pvalue)
            data["U value recovery"].append(delta_recovery_stat.statistic)
            data["p value recovery"].append(delta_recovery_stat.pvalue)

            data["Fly Type"].append("Test")
            data["Color"].append(filtered.loc[0, "Color"])
            data["On time(s)"].append(filtered.loc[0, "On time(s)"])
            data["Off2 protocol"].append(filtered.loc[0, "Off2 protocol"])
            data["Off1 protocol"].append(filtered.loc[0, "Off1 protocol"])

            data["Mean Off2_potential - Off1_potential"].append(mean_delta_v_test)
            data["Mean Off2_recovery - Off1_recovery"].append(mean_delta_recovery_test)

            data["U value depol"].append(delta_v_stat.statistic)
            data["p value depol"].append(delta_v_stat.pvalue)
            data["U value recovery"].append(delta_recovery_stat.statistic)
            data["p value recovery"].append(delta_recovery_stat.pvalue)

    return data

off_averaged_data = average_off_across_conditions(potential_df, recoveries_df)
off_averaged_df = pd.DataFrame(off_averaged_data)


# Write workbook
with pd.ExcelWriter("off_time_data.xlsx", engine="xlsxwriter") as w:
    potential_df.to_excel(w, index=False, sheet_name="OFF_Delta_potentials_PerFly")
    recoveries_df.to_excel(w, index=False, sheet_name="OFF_Delta_recoveries_PerFly")
    off_averaged_df.to_excel(w, index=False, sheet_name="OFF_averaged_data")
