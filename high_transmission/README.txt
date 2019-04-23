High transmission burden reduction use scripts
Prashanth Selvaraj, March 2019

-------------
Requirements:
dtk-tools package
dtk-tools-malaria package
malaria-toolbox package
These are available upon request from support@idmod.org

COMPS system for HPC job management
Email support@idmod.org

Input files, executable, and DLL's are elsewhere in the Additional File.

-------------
Scripts:
fig1_serialization_file.py - creates a population with immune profile matching a high transmission setting with no internvetions through which different interventions scenarios can be explored

fig1_main_run_file.py - creates and runs endectocide scenarios

fig1_plotting.py - generates plots for fig1

addendum_file1.py - generates plots for addendum file 1

analyzers/inset_chart_analyzer.py - analyze output files to generate EIR reduction

analyzers/run_analysis.py - collates simulation output for multiple scenarios and saves data in csv format

analyzers/summary_chart_analyzer.py - analyze output files to generate clinical case reduction by age or in the total population
