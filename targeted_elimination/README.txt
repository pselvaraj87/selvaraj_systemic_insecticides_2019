Targeted elimination use case scripts
Jaline Gerardin, March 2019

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
run_scenarios.py - creates and runs endectocide scenarios.

configure_forest_system.py - vector param updates, other config param updates, set up seasonal migration, set
input/exe/dll assets for HPC (helper to run_scenarios.py)

analyze_scenarios.py - grabs simulation output for multiple scenarios and packages relevant data into csv

plot_scenarios.py - generates plots of output from analyze_scenarios.py
