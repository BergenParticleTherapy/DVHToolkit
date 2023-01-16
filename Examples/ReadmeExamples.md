# Instructions for example files

## ExampleCohort1:
In this cohort the toxicity is defined in the ``ExampleCohort1_tox.csv`` file, followed by the setting "Load toxicity from: CSV_file":
This ``csv`` file contains rows
```
filename, toxlevel
```
where the ``filename``s reflect one patient each.

Load the file by configuring ``DVHToolkit`` with "DVH File type: simple"; "decimal/column separator: autodetect" or "./," and "dose unit: autodetect" or "cGy". Skip 1 row.

## ExampleCohort2:
The data here is the same, but instead of providing an ``ExampleCohort2_tox.csv`` file, the status is coded directly into the filenames themselves:
```
patient1.csv, patient2tox.csv
```
Configure the program similarly, but remember to set "Load toxicity from: 'tox' in filename".